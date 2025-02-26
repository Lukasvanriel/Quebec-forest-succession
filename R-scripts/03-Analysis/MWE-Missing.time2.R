# The first 155 lines contains code to generate an artificial dataset, and transformations
# to a dataset suitable for INLAjoint (data_inla). Events are coded as interval censored (value 3)
# lines 156 - 178 convert data_inla to a list where each element contains the specific transitions (event.list)
# 179 onwards contains code that creates the list of survival objects Surv.list (setting time2 as the end of the intervals)
# 195 onwards attempts to run the model using 'joint'
if(F){
install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE) 
install.packages("R.rsp") # (only for the vignette)

devtools::install_github('DenisRustand/INLAjoint', build_vignettes = TRUE)
}
### 0: Load packages ####
library(tidyverse)
library(INLA)
library(INLAjoint)
library(here)

### 1: Create artificial observation dataset ####
# Define the number of subjects and states
set.seed(123) 
n_subjects <- 1000
n_states <- 3

# I will be using 2 covariates (1 cat, 1 cont) to test
# I define the effects of the covariates through the transition matrices:
# The continuous variable has no effect, the categorical does

# Define the transition matrices. 
transition_matrix1 <- matrix(c( #Baseline
  0.94, 0.03, 0.03,  # Probabilities from state 1 to state 1, 2, 3
  0.03, 0.91, 0.06,  # Probabilities from state 2 to state 1, 2, 3
  0.06, 0.04, 0.9   # Probabilities from state 3 to state 1, 2, 3
), nrow = n_states, byrow = TRUE)

# For cat 2, p(1,2) and p(3,2) are increased, p(2,3) is decreased
transition_matrix2 <- matrix(c( #Baseline
  0.91, 0.06, 0.03,  # Probabilities from state 1 to state 1, 2, 3
  0.03, 0.94, 0.03,  # Probabilities from state 2 to state 1, 2, 3
  0.06, 0.08, 0.86   # Probabilities from state 3 to state 1, 2, 3
), nrow = n_states, byrow = TRUE)

# For cat 3, p(1,2) and p(1,3) are decreased, p(2,3) is decreased 
transition_matrix3 <- matrix(c( #Baseline
  0.98, 0.01, 0.01,  # Probabilities from state 1 to state 1, 2, 3
  0.03, 0.93, 0.04,  # Probabilities from state 2 to state 1, 2, 3
  0.06, 0.04, 0.9   # Probabilities from state 3 to state 1, 2, 3
), nrow = n_states, byrow = TRUE)

# Function to simulate Markov chain for one subject
simulate_markov_chain <- function(n_steps, transition_matrix) {
  states <- numeric(n_steps)
  states[1] <- sample(1:n_states, 1)  # Initial state chosen randomly
  
  for (i in 2:n_steps) {
    current_state <- states[i - 1]
    states[i] <- sample(1:n_states, 1, prob = transition_matrix[current_state, ])
  }
  
  return(states)
}

# Simulate Markov chains for all subjects
n_steps <- 50  # Number of steps/observations per subject
data1 <- matrix(0, nrow = n_subjects%/%3, ncol = n_steps)
data2 <- matrix(0, nrow = n_subjects%/%3, ncol = n_steps)
data3 <- matrix(0, nrow = n_subjects%/%3, ncol = n_steps)

for (subject in 1:n_subjects%/%3) {
  data1[subject, ] <- simulate_markov_chain(n_steps, transition_matrix1)
  data2[subject, ] <- simulate_markov_chain(n_steps, transition_matrix2)
  data3[subject, ] <- simulate_markov_chain(n_steps, transition_matrix3)
}

# Convert to dataframe
df <- as.data.frame(rbind(data1, data2, data3))
colnames(df) <- paste("Step", 1:n_steps, sep = "_")

# Generate covariates
df$ContinuousCovariate <- rnorm(n_subjects-1, mean = 0, sd = 1) # No significant effect
df$CategoricalCovariate <- as.factor(rep(c("Level1", "Level2", "Level3"), each=n_subjects%/%3))

# Select only a few observations:
df_select <- df[,c(1, 14, 23, 35, 49, (ncol(df)-1):ncol(df))]

### 2: Convert this to a dataset I can use for INLA: ####
## I first transform it to a format the msm package can use (this is where I'm coming from before trying INLAjoint)
data_msm <- df_select %>% 
  pivot_longer(cols = starts_with("Step"),
               names_to = "Time",
               names_prefix = "Step_",
               values_to = "State") %>%
  mutate(Time=(as.numeric(Time)-1)) %>% 
  mutate(ID=rep(1:(n_subjects-1), each=5)) %>% 
  arrange(ID, as.numeric(Time)) %>% 
  select(ID, Time, State, ContinuousCovariate, CategoricalCovariate)

## Now transform it to a format INLA can use 
# Some helper functions
expand.trans <- function(line.fr, line.to, transitions = 1:3){
  
  expansions <- lapply(transitions[-line.fr$State], FUN = function(x) {
    lines <- data.frame(ID = line.fr$ID,
                        Tstart = line.fr$Time, 
                        Tstop = line.to$Time,
                        from = line.fr$State,
                        to = x,
                        Status = ifelse(x == line.to$State, 3, 0),
                        ContCov = line.fr$ContinuousCovariate,
                        CatCov = line.fr$CategoricalCovariate ) 
  } )
  bind_rows(expansions)
}

expand.stay <- function(line.fr, line.to, transitions = 1:3){
  inclu <- list(transitions[-line.fr$State])[[1]]
  expansions <- lapply(inclu, FUN = function(x) {
    lines <- data.frame(ID = line.fr$ID,
                        Tstart = line.fr$Time, 
                        Tstop = line.to$Time,
                        from = line.fr$State,
                        to = x,
                        Status = 0,
                        ContCov = line.fr$ContinuousCovariate,
                        CatCov = line.fr$CategoricalCovariate ) 
  } )
  bind_rows(expansions)
}

expand.all <- function(block){
  expansion <- lapply(2:nrow(block), FUN = function(x) {
    if(block$State[x] == block$State[x - 1]) {
      no.fut.change <- all(block$State[x:nrow(block)] == block$State[x])
      if(no.fut.change) {
        if(x == nrow(block)){
          line.cor <- block[x - 1,]
          r = 2
          while((x - r) > 0 && block$State[x - 1] == block$State[x - r]) {
            line.cor <- block[x - r,]
            r <- r + 1
          }
          expand.stay(line.cor, block[x,])}
      } else {
        first.fut <- block$State[x:nrow(block)][which(block$State[x:nrow(block)] != block$State[x])[1]]
        expand.stay(block[x - 1,], block[x,]) # Only for future state because rest will be taken care of then as well
      }
    } else {
      line.cor <- block[x - 1,]
      
      expand.trans(line.cor, block[x,])
    }
  })
  bind_rows(expansion)
}

# Run the functions to convert to msm format to inla format
data_inla <- bind_rows(lapply(unique(data_msm$ID), FUN = function(x) {
  expand.all(data_msm %>% filter(ID == x))
} ))
head(data_inla, 10)

### 3: Convert to survival objects ####
Nb.states <- 3
Surv.list <- vector("list", Nb.states * (Nb.states - 1))

# Helper function to extract the 'from' and 'to' states in case of transition number k and amount of states N
get_transition_states <- function(k, N) {
  from_state <- (k - 1) %/% (N - 1) + 1
  position <- (k - 1) %% (N - 1)
  
  to_states <- setdiff(1:N, from_state)  # Exclude from-state
  to_state <- to_states[position + 1]   # Pick the corresponding to-state
  
  return(c(from_state, to_state))
}

# Create a list, where each element contains the specific transitions (e.g. [[1]] contains all 1->2)
event.list <- lapply(1:length(Surv.list), FUN = function(x){
  states.involved <- get_transition_states(x, Nb.states)
  data_inla %>% 
    filter(from == states.involved[1],
           to == states.involved[2])
})

### 4: Define the interval times, including time2 ####
# Surv.list is a list of survival objects, where [[1]] is for transition 1->2, [[2]] for 1->3 etc.
for(i in 1:length(Surv.list)) { 
  Surv.list[[i]] <- inla.surv(time = event.list[[i]]$Tstart,
                              event = event.list[[i]]$Status,
                              time2 = event.list[[i]]$Tstop) # Here I define time2, not providing it here throws an
}

# The joint function only seems to work if I store the survival objects 
s12 <- Surv.list[[1]]
s13 <- Surv.list[[2]]
s21 <- Surv.list[[3]]
s23 <- Surv.list[[4]]
s31 <- Surv.list[[5]]
s32 <- Surv.list[[6]]

### 5: Run the joint function ####
weib.surv <- joint(formSurv = list(
                        s12 ~ ContCov + CatCov,
                        s13 ~ ContCov + CatCov,
                        s21 ~ ContCov + CatCov,
                        s23 ~ ContCov + CatCov,
                        s31 ~ ContCov + CatCov,
                        s32 ~ ContCov + CatCov),
  basRisk = rep("weibullsurv", 6), 
  dataSurv = event.list,
  control = list(config=TRUE)) #, verbose=TRUE
summary(weib.surv)

