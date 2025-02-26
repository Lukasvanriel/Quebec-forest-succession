#### Load packages ####
library(tidyverse)
library(msm)
library(INLA)
library(INLAjoint)
library(here)
library(funrar)
library(p3state.msm)

#### Create artificial dataset ####

# Load necessary library
set.seed(123) # For reproducibility

# Define the number of subjects and states
n_subjects <- 1000
n_states <- 3

# Define the transition matrix. We will be using 2 covariates (1 cat, 1 cont)
# I define the effects of the covariates on through the transition matrices:
# For now the continuous has no effect
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

(transition_matrix1 + transition_matrix2 + transition_matrix3) / 3

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

# Convert to data frame for easier handling
df <- as.data.frame(rbind(data1, data2, data3))
colnames(df) <- paste("Step", 1:n_steps, sep = "_")

# Generate covariates
df$ContinuousCovariate1 <- rnorm(n_subjects-1, mean = 0, sd = 1) # No significant effect
df$CategoricalCovariate <- as.factor(rep(c("Level1", "Level2", "Level3"), each=n_subjects%/%3))

# Select only a few observations:
df_select <- df[,c(1, 14, 23, 35, 49, (ncol(df)-1):ncol(df))]

# Print first few rows of the dataset
head(df_select)

#### Wrangle ####
### msm
data_msm <- df_select %>% 
  pivot_longer(cols = starts_with("Step"),
               names_to = "Time",
               names_prefix = "Step_",
               values_to = "State") %>%
  mutate(Time=(as.numeric(Time)-1)) %>% 
  mutate(ID=rep(1:(n_subjects-1), each=5)) %>% 
  arrange(ID, as.numeric(Time)) %>% 
  select(ID, Time, State, ContinuousCovariate1, CategoricalCovariate)

### INLA
## Functions
expand.trans <- function(line.fr, line.to, transitions = 1:3){
  
  expansions <- lapply(transitions[-line.fr$State], FUN = function(x) {
    lines <- data.frame(ID = line.fr$ID,
                        Tstart = line.fr$Time, 
                        Tstop = line.to$Time,
                        from = line.fr$State,
                        to = x,
                        Status = ifelse(x == line.to$State, 3, 0),
                        ContCov = line.fr$ContinuousCovariate1,
                        CatCov = line.fr$CategoricalCovariate ) 
    } )
  bind_rows(expansions)
}

expand.stay <- function(line.fr, line.to, future = -1, transitions = 1:3){
  inclu <- list(transitions[-line.fr$State])[[1]]
  #inclu <- ifelse(future == -1, list(transitions[-line.fr$State]), list(future))[[1]]
  expansions <- lapply(inclu, FUN = function(x) {
    lines <- data.frame(ID = line.fr$ID,
                        Tstart = line.fr$Time, 
                        Tstop = line.to$Time,
                        from = line.fr$State,
                        to = x,
                        Status = 0,
                        ContCov = line.fr$ContinuousCovariate1,
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
          expand.stay(line.cor, block[x,], future = -1)}
      } else {
        first.fut <- block$State[x:nrow(block)][which(block$State[x:nrow(block)] != block$State[x])[1]]
        expand.stay(block[x - 1,], block[x,], first.fut) # Only for future state because rest will be taken care of then as well
      }
    } else {
      line.cor <- block[x - 1,]

      expand.trans(line.cor, block[x,])
    }
  })
  bind_rows(expansion)
}

if(F){
  
  expand.stay <- function(line.fr, line.to, future = -1, transitions = 1:3){
    inclu <- ifelse(future == -1, list(transitions[-line.fr$State]), list(future))[[1]]
    expansions <- lapply(inclu, FUN = function(x) {
      lines <- data.frame(ID = line.fr$ID,
                          Tstart = line.fr$Time, 
                          Tstop = line.to$Time,
                          from = line.fr$State,
                          to = x,
                          Status = 0,
                          ContCov = line.fr$ContinuousCovariate1,
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
            expand.stay(line.cor, block[x,], future = -1)}
        } else{
          first.fut <- block$State[x:nrow(block)][which(block$State[x:nrow(block)] != block$State[x])[1]]
          expand.stay(block[x - 1,], block[x,], first.fut) # Only for future state because rest will be taken care of then as well
        }
      } else {
        line.cor <- block[x - 1,]
        r = 2
        while((x - r) > 0 && block$State[x - 1] == block$State[x - r]) {
          line.cor <- block[x - r,]
          r <- r + 1
        }
        expand.trans(line.cor, block[x,])
      }
    })
    bind_rows(expansion)
  }
} # Older versions


## The data
data_inla <- bind_rows(lapply(unique(data_msm$ID), FUN = function(x) {
  expand.all(data_msm %>% filter(ID == x))
} ))

## Convert to survival objects:
get_transition_states <- function(k, N) {
  from_state <- (k - 1) %/% (N - 1) + 1
  position <- (k - 1) %% (N - 1)
  
  to_states <- setdiff(1:N, from_state)  # Exclude from-state
  to_state <- to_states[position + 1]   # Pick the corresponding to-state
  
  return(c(from_state, to_state))
}

Nb.states <- 3
Surv.list <- vector("list", Nb.states * (Nb.states - 1))

event.list <- lapply(1:length(Surv.list), FUN = function(x){
  data_inla %>% 
    filter(from == get_transition_states(x, Nb.states)[1],
           to == get_transition_states(x, Nb.states)[2])
})

for(i in 1:length(Surv.list)) { 
  Surv.list[[i]] <- inla.surv(time = event.list[[i]]$Tstart,
                              event = event.list[[i]]$Status,
                              time2 = event.list[[i]]$Tstop)
}

?inla.surv

s12 <- Surv.list[[1]]
s13 <- Surv.list[[2]]
s21 <- Surv.list[[3]]
s23 <- Surv.list[[4]]
s31 <- Surv.list[[5]]
s32 <- Surv.list[[6]]

#### Run models ####
### INLA

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

weib.surv <- joint(formSurv = list(
  inla.surv(time = event.list[[1]]$Tstart, time2 = event.list[[1]]$Tstop, event = event.list[[1]]$Status) ~ ContCov + CatCov,
  inla.surv(time = event.list[[2]]$Tstart, time2 = event.list[[2]]$Tstop, event = event.list[[2]]$Status) ~ ContCov + CatCov,
  inla.surv(time = event.list[[3]]$Tstart, time2 = event.list[[3]]$Tstop, event = event.list[[3]]$Status) ~ ContCov + CatCov,
  inla.surv(time = event.list[[4]]$Tstart, time2 = event.list[[4]]$Tstop, event = event.list[[4]]$Status) ~ ContCov + CatCov,
  inla.surv(time = event.list[[5]]$Tstart, time2 = event.list[[5]]$Tstop, event = event.list[[5]]$Status) ~ ContCov + CatCov,
  inla.surv(time = event.list[[6]]$Tstart, time2 = event.list[[6]]$Tstop, event = event.list[[6]]$Status) ~ ContCov + CatCov),
  basRisk = rep("weibullsurv", 6), 
  dataSurv = event.list,
  control = list(config=TRUE))
summary(weib.surv)

Surv.list[[i]] <- inla.surv(time = event.list[[i]]$Tstart,
                            time2 = event.list[[i]]$Tstop,
                            event = event.list[[i]]$Status)

# No covariate:
weib.surv.base <- joint(formSurv = list(
  s12 ~ 1,
  s13 ~ 1,
  s21 ~ 1,
  s23 ~ 1,
  s31 ~ 1,
  s32 ~ 1),
  basRisk = rep("weibullsurv", 6), dataSurv = event.list,
  control = list(config=TRUE))
summary(weib.surv.base)

exp.surv <- joint(formSurv = list(
  s12 ~ CatCov,
  s13 ~ CatCov,
  s21 ~ CatCov,
  s23 ~ CatCov,
  s31 ~ CatCov,
  s32 ~ CatCov),
  basRisk = rep("exponentialsurv", 6), dataSurv = event.list,
  control = list(config=TRUE))
summary(exp.surv) 

exp.surv.base <- joint(formSurv = list(
  s12 ~ 1,
  s13 ~ 1,
  s21 ~ 1,
  s23 ~ 1,
  s31 ~ 1,
  s32 ~ 1),
  basRisk = rep("exponentialsurv", 6), dataSurv = event.list,
  control = list(config=TRUE))
summary(exp.surv.base)

exp.surv <- joint(formSurv = list(
  s12 ~ CatCov,
  s13 ~ CatCov,
  s21 ~ CatCov,
  s23 ~ CatCov,
  s31 ~ CatCov,
  s32 ~ CatCov),
  basRisk = rep("exponentialsurv", 6), dataSurv = event.list,
  control = list(config=TRUE))
summary(exp.surv) 

exp.surv.base <- joint(formSurv = list(
  s12 ~ 1,
  s13 ~ 1,
  s21 ~ 1,
  s23 ~ 1,
  s31 ~ 1,
  s32 ~ 1),
  basRisk = rep("exponentialsurv", 6), dataSurv = event.list,
  control = list(config=TRUE))
summary(exp.surv.base)

# Set the time
t <- seq(0,100, by=1)

# Create statetable
msm_state <- statetable.msm(State, ID, data=data_msm)
round(funrar::make_relative(msm_state), 3)

### msm
# Define the model (all transitions are possible here)
Q.model <- rep(0.3, 9)
dim(Q.model) <- c(3,3)

# Get initial estimate for Q

Q.init <- crudeinits.msm(State ~ Time, ID, data = data_msm, qmatrix = Q.model)

# Run base msm
output.msm.base.CG <-  msm(State ~ Time, ID, data = data_msm, qmatrix = Q.init, method= "CG",
                           control = list(fnscale = 100000, maxit=1000))
output.msm.base.CG

output.msm.cov.CG <-  msm(State ~ Time, ID, data=data_msm, qmatrix = Q.init, method= "CG",
                          control = list(fnscale = 100000,  maxit=100),
                          covariates = ~ CategoricalCovariate + ContinuousCovariate)
output.msm.cov.CG


#### Create output objects ####

# msm: Get a list of p-matrices
covariates <- list(c(CategoricalCovariate="Level1"), c(CategoricalCovariate="Level2"), c(CategoricalCovariate="Level3"))

# For normal
pmatrix.msm(output.msm.cov.CG, max(t), covariates = as.list(covariates[[1]]))

p <- lapply(covariates, FUN=function(c) lapply(1:max(t), FUN = function(t) pmatrix.msm(output.msm.cov.CG, t, covariates = as.list(c))))

# Convert into dataframe
p.vector <- lapply(p, FUN = function(x) lapply(x, FUN = function(z) as.vector(t(z))))
p.msm.df <- lapply(1:3, function(x) data.frame(matrix(ncol = (length(p.vector[[1]][[1]]) + 1), nrow = max(t))))
p.msm.df <- lapply(1:3, function(y) setNames(as.data.frame(cbind(1:max(t), sapply(1:(ncol(p.msm.df[[y]])-1), FUN=function(x) p.msm.df[[y]][,x+1] <- sapply(p.vector[[y]], FUN=function(z) z[x])))),
                                             c("time", "t11", "t12", "t13", "t21", "t22", "t23", "t31", "t32", "t33")))

# For exponential
p <- lapply(1:max(t), FUN = function(t) pmatrix.msm(output.msm.base.CG, t))
# Convert into dataframe
p.vectorE <- lapply(p, FUN = function(z) as.vector(t(z)))
p.msm.dfE <- data.frame(matrix(ncol = (length(p.vectorE[[1]]) + 1), nrow = max(t)))
p.msm.dfE <- setNames(as.data.frame(cbind(1:max(t), sapply(1:(ncol(p.msm.dfE)-1), FUN=function(x) p.msm.dfE[,x+1] <- sapply(p.vectorE, FUN=function(z) z[[x]])))),
                     c("time", "t11", "t12", "t13", "t21", "t22", "t23", "t31", "t32", "t33"))


# INLA: Q and P matrices
# Functions
dP_dt_full <- function(t, P0, params) {
  # PREP: define initial conditions and get all get all hij(t)
  P_A <- P0[1]
  P_B <- P0[2]
  P_C <- P0[3]
  
  all_transitions <- paste0(rep(LETTERS[1:3], each=3), rep(LETTERS[1:3], 3))
  
  h_all <- sapply(all_transitions, FUN = function(x) h(x, t, params))
  dim(h_all) <- c(3, 3) 
  h_all[1,1] <- h_all[2,2] <- h_all[3,3] <- 0
  
  #EQUATIONS
  
  dPA <- - sum(h_all[1,], na.rm = T) * P_A + 
    h_all[2,1] * P_B + h_all[3,1] * P_C 
  
  dPB <- - sum(h_all[2,], na.rm = T) * P_B + 
    h_all[1,2] * P_A + h_all[3,2] * P_C
  dPC <- - sum(h_all[3,], na.rm = T) * P_C + 
    h_all[1,3] * P_A + h_all[2,3] * P_B
  
  
  
  #OUTPUT AS LIST
  list(c(dPA, dPB, dPC))
}

Wb <- function(t, X, l, a, b) {
  l * a * t^(a - 1) * exp(sum(b * X))
}

Exp <- function(t, X, l, b) {
  l * exp(sum(b * X))
}

h <- function(string, time, params){
  #convert the letters to a transition number between 1 and 6
  f <- match(substr(string, 1,1), LETTERS)
  t <- match(substr(string, 2,2), LETTERS)
  
  if(t == f) {return(-1)}
  
  helper <- seq(1,3, by=1)
  transitions <- c()
  for(i in helper) {
    for(j in helper[-i]) {
      transitions <- c(transitions, 10*i + j%%10)
    }
  }
  info.transitions <- data.frame(inla_trans=seq(1:6), from = transitions %/% 10, to = transitions %% 10)
  
  transition_number <- info.transitions %>% 
    filter(from==f, to==t) %>% 
    pull(inla_trans)
  
  # Based on the transition number, lets now get the necessary parameters to feed to the Weibull hazard
  
  X <- params$X
  l <- exp(params$beta[names(params$beta) %in% paste0(c("Interc"), transition_number)])
  a <- params$alpha[transition_number]
  b <- params$beta[names(params$beta) %in% paste0(c("ContCov", "CatCov2", "CatCov3"), transition_number)]
  
  Wb(time, X, l, a, b)
}

generate_Q <- function(t, params){
  n_states <- 3
  all_transitions <- paste0(rep(LETTERS[1:n_states], each=n_states), rep(LETTERS[1:n_states], n_states))
  h_all <- lapply(all_transitions, FUN = function(x) h(x, t, params))
  h_all <- matrix(h_all, ncol=3, nrow=3, byrow=TRUE)
  
  sapply(1:n_states, function(i) {h_all[i,i] <<- list(-mapply(`+`, h_all[i,seq(1,n_states)[-i][1]][[1]],
                                                              h_all[i,seq(1,n_states)[-i][2]][[1]]))})
  
  h_all
}

h.exp <- function(string, time, params){
  #convert the letters to a transition number between 1 and 6
  f <- match(substr(string, 1,1), LETTERS)
  t <- match(substr(string, 2,2), LETTERS)
  
  if(t == f) {return(-1)}
  
  helper <- seq(1,3, by=1)
  transitions <- c()
  for(i in helper) {
    for(j in helper[-i]) {
      transitions <- c(transitions, 10*i + j%%10)
    }
  }
  info.transitions <- data.frame(inla_trans=seq(1:6), from = transitions %/% 10, to = transitions %% 10)
  
  transition_number <- info.transitions %>% 
    filter(from==f, to==t) %>% 
    pull(inla_trans)
  
  # Based on the transition number, lets now get the necessary parameters to feed to the Weibull hazard
  
  X <- params$X
  l <- exp(params$beta[names(params$beta) %in% paste0(c("Interc"), transition_number)])
  b <- params$beta[names(params$beta) %in% paste0(c("CatCov2", "CatCov3"), transition_number)]
  
  Exp(time, X, l, b)
}

generate_QExp <- function(t, params){
  n_states <- 3
  all_transitions <- paste0(rep(LETTERS[1:n_states], each=n_states), rep(LETTERS[1:n_states], n_states))
  h_all <- lapply(all_transitions, FUN = function(x) h.exp(x, t, params))
  h_all <- matrix(unlist(h_all), ncol=3, nrow=3, byrow=TRUE) #for now like this to lose names
  
  sapply(1:n_states, function(i) {h_all[i,i] <<- -mapply(`+`, h_all[i,seq(1,n_states)[-i][1]],
                                                         h_all[i,seq(1,n_states)[-i][2]])})
  
  h_all
}

# Gather the output parameters
t <- 1:60
parameters <- list()

parameters$alpha <- weib.surv$summary.hyperpar$mean
names(parameters$alpha) <- as.character(1:6)

parameters$beta <- weib.surv$summary.fixed$mean
names(parameters$beta) <- paste0(rep(c("Interc", "ContCov", "CatCov2", "CatCov3"), each=1, 6), rep(1:6, each=4))

parameters$X <- c(0,0,0)

# Create Q and P matrices
a <- generate_Q(0:max(t), parameters)

Q <- lapply(1:(max(t)+1), function(i) {matrix(sapply(a, function(m) m[[i]]), nrow = nrow(a), ncol = ncol(a))})

P <- vector("list", length = length(t))
p <- diag(3)
sapply(1:(max(t)+1), function(x) {
  P[[x]] <<- p %*% expm(Q[[x]])
  p <<- P[[x]]
} )

p.vector <- lapply(P, function(p) {
  p <- as.vector(t(p))
  names(p) <- c("p11", "p12", "p13", "p21", "p22", "p23", "p31", "p32", "p33")
  p
} )
p.df <- bind_rows(p.vector)

parameters$X <- c(0,1,0)
a2 <- generate_Q(0:max(t), parameters)

Q2 <- lapply(1:(max(t)+1), function(i) {matrix(sapply(a2, function(m) m[[i]]), nrow = nrow(a2), ncol = ncol(a2))})

P2 <- vector("list", length = length(t))
p2 <- diag(3)
sapply(1:(max(t)+1), function(x) {
  P2[[x]] <<- p2 %*% expm(Q2[[x]])
  p2 <<- P2[[x]]
} )

p.vector2 <- lapply(P2, function(p) {
  p <- as.vector(t(p))
  names(p) <- c("p11", "p12", "p13", "p21", "p22", "p23", "p31", "p32", "p33")
  p
} )
p.df2 <- bind_rows(p.vector2)

parameters$X <- c(0,0,1)
a3 <- generate_Q(0:max(t), parameters)

Q3 <- lapply(1:(max(t)+1), function(i) {matrix(sapply(a3, function(m) m[[i]]), nrow = nrow(a3), ncol = ncol(a3))})

P3 <- vector("list", length = length(t))
p3 <- diag(3)
sapply(1:(max(t)+1), function(x) {
  P3[[x]] <<- p3 %*% expm(Q3[[x]])
  p3 <<- P3[[x]]
} )

p.vector3 <- lapply(P3, function(p) {
  p <- as.vector(t(p))
  names(p) <- c("p11", "p12", "p13", "p21", "p22", "p23", "p31", "p32", "p33")
  p
} )
p.df3 <- bind_rows(p.vector3)

## without covariates:

parameters <- list()

parameters$alpha <- weib.surv.base$summary.hyperpar$mean
names(parameters$alpha) <- as.character(1:6)

parameters$beta <- weib.surv.base$summary.fixed$mean
names(parameters$beta) <- paste0(rep(c("Interc"), each=1, 6), rep(1:6, each=1))

parameters$X <- c(0,0,0)

#
an <- generate_Q(0:max(t), parameters)

Qn <- lapply(1:(max(t)+1), function(i) {matrix(sapply(an, function(m) m[[i]]), nrow = nrow(an), ncol = ncol(an))})

Pn <- vector("list", length = length(t))
pn <- diag(3)
sapply(1:(max(t)+1), function(x) {
  Pn[[x]] <<- pn %*% expm(Qn[[x]])
  pn <<- Pn[[x]]
} )

p.vectorn <- lapply(Pn, function(p) {
  p <- as.vector(t(p))
  names(p) <- c("p11", "p12", "p13", "p21", "p22", "p23", "p31", "p32", "p33")
  p
} )
p.dfn <- bind_rows(p.vectorn)

#Exponentials

parameters$beta <- exp.surv$summary.fixed$mean
names(parameters$beta) <- paste0(rep(c("Interc", "CatCov2", "CatCov3"), each=1, 6), rep(1:6, each=3))

parameters$X <- c(0,0)

QE1 <- generate_QExp(t, parameters)
expm(QE1)

PE1 <- vector("list", length = length(t))
pE1 <- diag(3)

sapply(1:(max(t)+1), function(x) {
  PE1[[x]] <<- pE1 %*% expm(QE1)
  pE1 <<- PE1[[x]]
} )

parameters$X <- c(1,0)

QE2 <- generate_QExp(t, parameters)
expm(QE2)

PE2 <- vector("list", length = length(t))
pE2 <- diag(3)

sapply(1:(max(t)+1), function(x) {
  PE2[[x]] <<- pE2 %*% expm(QE2)
  pE2 <<- PE2[[x]]
} )

parameters$X <- c(0,1)

QE3 <- generate_QExp(t, parameters)
expm(QE3)

PE3 <- vector("list", length = length(t))
pE3 <- diag(3)

sapply(1:(max(t)+1), function(x) {
  PE3[[x]] <<- pE3 %*% expm(QE3)
  pE3 <<- PE3[[x]]
} )

p.vectorE1 <- lapply(PE1, function(p) {
  p <- as.vector(t(p))
  names(p) <- c("p11", "p12", "p13", "p21", "p22", "p23", "p31", "p32", "p33")
  p
} )
p.dfE1 <- bind_rows(p.vectorE1)

p.vectorE2 <- lapply(PE2, function(p) {
  p <- as.vector(t(p))
  names(p) <- c("p11", "p12", "p13", "p21", "p22", "p23", "p31", "p32", "p33")
  p
} )
p.dfE2 <- bind_rows(p.vectorE2)

p.vectorE3 <- lapply(PE3, function(p) {
  p <- as.vector(t(p))
  names(p) <- c("p11", "p12", "p13", "p21", "p22", "p23", "p31", "p32", "p33")
  p
} )
p.dfE3 <- bind_rows(p.vectorE3)


#### Check out base results / Analysis ####
## Look at some results:
# msm
qmatrix.msm(output.msm.base.CG)

pmatrix.msm(output.msm.base.CG, 1)
pmatrix.msm(output.msm.base.CG, 10)

pnext.msm(output.msm.base.CG)

# with covariates
qmatrix.msm(output.msm.cov.CG)
pmatrix.msm(output.msm.cov.CG, 10, ci="normal")
pmatrix.msm(output.msm.cov.CG, 1, covariates = list(CategoricalCovariate="Level1")) %>% round(2)
pmatrix.msm(output.msm.cov.CG, 1, covariates = list(CategoricalCovariate="Level2")) %>% round(2)
pmatrix.msm(output.msm.cov.CG, 1, covariates = list(CategoricalCovariate="Level3")) %>% round(2)

HR.msm.cov.CG <- hazard.msm(output.msm.cov.CG)
HR.msm.cov.CG

#INLA

log(HR.msm.cov.CG$CategoricalCovariateLevel2)
weib.surv$summary.fixed$mean

plot(weib.surv)$Baseline
plot(weib.surv)$Outcomes

# What after 10 iterations?
p.df[11,]
p.msm.df[[1]][10,2:10]

p.df2[11,]
p.msm.df[[2]][10,2:10]

p.df3[11,]
p.msm.df[[3]][10,2:10]

#### Plots #####
# First let's compare msm with full weibull survival for starting in state1:
plot(0:max(t), p.df$p11, pch=20, ylim=c(0,1))
points(0:max(t), p.df$p12, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.df$p13, pch=20, ylim=c(0,1), col="blue")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t11, pch=20, col="darkgreen")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t12, pch=20, col="orange")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t13, pch=20, col="lightblue")
abline(v = 10)

plot(0:max(t), p.df$p22, pch=20, ylim=c(0,1))
points(0:max(t), p.df$p21, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.df$p23, pch=20, ylim=c(0,1), col="blue")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t22, pch=20, col="darkgreen")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t21, pch=20, col="orange")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t23, pch=20, col="lightblue")

plot(0:max(t), p.df$p33, pch=20, ylim=c(0,1))
points(0:max(t), p.df$p31, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.df$p32, pch=20, ylim=c(0,1), col="blue")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t33, pch=20, col="darkgreen")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t31, pch=20, col="orange")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t32, pch=20, col="lightblue")

# Compare weibull with exponential:
plot(0:max(t), p.df$p11, pch=20, ylim=c(0,1))
points(0:max(t), p.df$p12, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.df$p13, pch=20, ylim=c(0,1), col="blue")
points(0:max(t), p.dfE1$p11, pch=20, col="darkgreen")
points(0:max(t), p.dfE1$p12, pch=20, col="orange")
points(0:max(t), p.dfE1$p13, pch=20, col="lightblue")
abline(v = 10)

plot(0:max(t), p.df$p22, pch=20, ylim=c(0,1))
points(0:max(t), p.df$p21, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.df$p23, pch=20, ylim=c(0,1), col="blue")
points(0:max(t), p.dfE1$p22, pch=20, col="darkgreen")
points(0:max(t), p.dfE1$p21, pch=20, col="orange")
points(0:max(t), p.dfE1$p23, pch=20, col="lightblue")

plot(0:max(t), p.df$p33, pch=20, ylim=c(0,1))
points(0:max(t), p.df$p31, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.df$p32, pch=20, ylim=c(0,1), col="blue")
points(0:max(t), p.dfE1$p33, pch=20, col="darkgreen")
points(0:max(t), p.dfE1$p31, pch=20, col="orange")
points(0:max(t), p.dfE1$p32, pch=20, col="lightblue")

# Compare msm with exponential:
plot(0:max(t), p.dfE1$p11, pch=20, ylim=c(0,1))
points(0:max(t), p.dfE1$p12, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.dfE1$p13, pch=20, ylim=c(0,1), col="blue")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t11, pch=20, col="darkgreen")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t12, pch=20, col="orange")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t13, pch=20, col="lightblue")
abline(v = 10)

plot(0:max(t), p.dfE1$p22, pch=20, ylim=c(0,1))
points(0:max(t), p.dfE1$p21, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.dfE1$p23, pch=20, ylim=c(0,1), col="blue")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t22, pch=20, col="darkgreen")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t21, pch=20, col="orange")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t23, pch=20, col="lightblue")

plot(0:max(t), p.dfE1$p33, pch=20, ylim=c(0,1))
points(0:max(t), p.dfE1$p31, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.dfE1$p32, pch=20, ylim=c(0,1), col="blue")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t33, pch=20, col="darkgreen")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t31, pch=20, col="orange")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t32, pch=20, col="lightblue")

# Covariates:
plot(0:max(t), p.df$p11, pch=20, ylim=c(0,1))
points(0:max(t), p.df2$p11, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.df3$p11, pch=20, ylim=c(0,1), col="blue")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t11, pch=20, col="darkgreen")
points(p.msm.df[[2]]$time, p.msm.df[[2]]$t11, pch=20, col="orange")
points(p.msm.df[[3]]$time, p.msm.df[[3]]$t11, pch=20, col="lightblue")
abline(v = 10)


plot(0:max(t), p.df$p11, pch=20, ylim=c(0,1))
points(0:max(t), p.df2$p11, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.df3$p11, pch=20, ylim=c(0,1), col="blue")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t11, pch=20, col="darkgreen")
points(p.msm.df[[2]]$time, p.msm.df[[2]]$t11, pch=20, col="orange")
points(p.msm.df[[3]]$time, p.msm.df[[3]]$t11, pch=20, col="lightblue")
abline(v = 10)


 ## Create some plots p(t)

plot(p.msm.df[[1]]$time, p.msm.df[[1]]$t12, pch=20, ylim = c(0, 1))
points(p.msm.df[[2]]$time, p.msm.df[[2]]$t12, pch=20, col='red')
points(p.msm.df[[3]]$time, p.msm.df[[3]]$t12, pch=20, col='blue')

# Make the plots

plots.msm.p <- function(p.df, tr.name){
  cbind(p.df[[1]]$time, p.df[[1]][tr.name], p.df[[2]][tr.name], p.df[[3]][tr.name]) %>% 
    as.data.frame() %>% 
    `colnames<-`(c("time", "level1", "level2", "level3")) %>% 
    pivot_longer(cols = c(level1, level2, level3), names_to = "Covariate") %>% 
    ggplot +
    geom_point(aes(x=time, y=value, color=Covariate))
}

plots.msm.p(p.msm.df, "t12")
plots.msm.p(p.msm.df, "t11")
plots.msm.p(p.msm.df, "t31")


# For cat 2, p(1,2) and p(3,2) are increased, p(2,3) is decreased
# For cat 3, p(1,2) and p(1,3) are decreased, p(2,3) is decreased


p.df$p11[11]
p.df$p12[11]
p.df$p13[11]

plot(0:max(t), p.df$p11, pch=20, ylim=c(0,1))
points(0:max(t), p.df$p12, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.df$p13, pch=20, ylim=c(0,1), col="blue")

points(p.msm.df[[1]]$time, p.msm.df[[1]]$t11, pch=20, col="darkgreen")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t12, pch=20, col="orange")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t13, pch=20, col="lightblue")


plot(0:max(t), p.df$p11, pch=20, ylim=c(0,1))
points(start1.X1[,1], start1.X1[,2], pch=20, col="red")

plot(0:max(t), p.df$p12, pch=20, ylim=c(0,1))
points(start1.X1[,1], start1.X1[,3], pch=20, col="red")

plot(0:max(t), p.df$p13, pch=20, ylim=c(0,1))
points(start1.X1[,1], start1.X1[,4], pch=20, col="red")


######### OLD ##########
# Some backups

expand.stay <- function(line.fr, line.to, future = -1, transitions = 1:3){
  inclu <- ifelse(future == -1, list(transitions[-line.fr$State]), list(future))[[1]]
  expansions <- lapply(inclu, FUN = function(x) {
    lines <- data.frame(ID = line.fr$ID,
                        Tstart = line.fr$Time, 
                        Tstop = line.to$Time,
                        from = line.fr$State,
                        to = x,
                        Status = 0,
                        ContCov = line.fr$ContinuousCovariate1,
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
          expand.stay(line.cor, block[x,], future = -1)}
      } else{
        first.fut <- block$State[x:nrow(block)][which(block$State[x:nrow(block)] != block$State[x])[1]]
        expand.stay(block[x - 1,], block[x,], first.fut) # Only for future state because rest will be taken care of then as well
      }
    } else {
      line.cor <- block[x - 1,]
      r = 2
      while((x - r) > 0 && block$State[x - 1] == block$State[x - r]) {
        line.cor <- block[x - r,]
        r <- r + 1
      }
      expand.trans(line.cor, block[x,])
    }
  })
  bind_rows(expansion)
}





#### Analysis ####
# Let's try the differential equations
library(deSolve)

# Set the initial conditions:

P0A <- c(P_A = 1, P_B = 0, P_C = 0)
P0B <- c(P_A = 0, P_B = 1, P_C = 0)
P0C <- c(P_A = 0, P_B = 0, P_C = 1)

# Show figures
#df

plot(start1.X1[,1], start1.X1[,2], pch=20, ylim=c(0,1), xlab="time [yr]", ylab="Pr")
points(start1.X2[,1], start1.X2[,2], pch=20, ylim=c(0,1), col="red")
points(start1.X3[,1], start1.X3[,2], pch=20, ylim=c(0,1), col="darkgreen")


# Compare with msm
plot(0:max(t), p.df$p12, pch=20, col="black", ylim=c(0,1), xlab="time [yr]", ylab="Pr")
points(0:max(t), p.df2$p12, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.df3$p12, pch=20, ylim=c(0,1), col="darkblue")

points(p.msm.df[[1]]$time, p.msm.df[[1]]$t12, pch=20, col="grey")
points(p.msm.df[[2]]$time, p.msm.df[[2]]$t12, pch=20, ylim=c(0,1), col="orange")
points(p.msm.df[[3]]$time, p.msm.df[[3]]$t12, pch=20, ylim=c(0,1), col="cyan")

# Compare with exponential
plot(0:max(t), p.df$p12, pch=20, col="black", ylim=c(0,1), xlab="time [yr]", ylab="Pr")
points(0:max(t), p.df2$p12, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.df3$p12, pch=20, ylim=c(0,1), col="darkblue")

points(0:max(t), p.dfE1$p12, pch=20, col="grey")
points(0:max(t), p.dfE2$p12, pch=20, ylim=c(0,1), col="orange")
points(0:max(t), p.dfE3$p12, pch=20, ylim=c(0,1), col="cyan")


plot(0:max(t), p.df$p11, pch=20, col="black", ylim=c(0,1), xlab="time [yr]", ylab="Pr")
points(0:max(t), p.df2$p11, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.df3$p11, pch=20, ylim=c(0,1), col="darkblue")

points(0:max(t), p.dfE1$p11, pch=20, col="grey")
points(0:max(t), p.dfE2$p11, pch=20, ylim=c(0,1), col="orange")
points(0:max(t), p.dfE3$p11, pch=20, ylim=c(0,1), col="cyan")

####### Different models ########


plot(0:max(t), p.df$p11, pch=20, ylim=c(0,1))
points(0:max(t), p.df$p12, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.df$p13, pch=20, ylim=c(0,1), col="red")



plots.inla.p <- function(p.df, tr.name){
    p.df %>% 
    as.data.frame() %>% 
    `colnames<-`(c("time", "level1", "level2", "level3")) %>% 
    pivot_longer(cols = c(level1, level2, level3), names_to = "Covariate") %>% 
    ggplot +
    geom_point(aes(x=time, y=value, color=Covariate))
}

plots.inla.p(cbind(time=0:max(t), level1=p.df["p13"], level2=p.df2["p13"], level3=p.df3["p13"]), "p13")
plots.msm.p(p.msm.df, "t13")


plots.inla.p(cbind(time=0:max(t), level1=p.df["p13"], level2=p.df2["p13"], level3=p.df3["p13"]), "p13")

plots.msm.p(p.msm.df, "t12")

#########
p <- lapply(1:max(t), FUN = function(t) pmatrix.msm(output.msm.base.CG, t))

# Convert into dataframe
p.vector <- lapply(p, FUN = function(z) as.vector(t(z)))
p.msm.df <- data.frame(matrix(ncol = (length(p.vector[[1]]) + 1), nrow = max(t)))
p.msm.df <- setNames(as.data.frame(cbind(1:max(t), sapply(1:(ncol(p.msm.df)-1), FUN=function(x) p.msm.df[,x+1] <- sapply(p.vector, FUN=function(z) z[[x]])))),
                                             c("time", "t11", "t12", "t13", "t21", "t22", "t23", "t31", "t32", "t33"))



plot(0:max(t), p.df$p11, pch=20, ylim=c(0,1))
points(0:max(t), p.df$p12, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.df$p13, pch=20, ylim=c(0,1), col="red")
points(p.msm.df$time, p.msm.df$t11, pch=20, col="darkblue")
points(p.msm.df$time, p.msm.df$t12, pch=20, ylim=c(0,1), col="orange")
points(p.msm.df$time, p.msm.df$t13, pch=20, ylim=c(0,1), col="orange")


## Just exponential, no covariates





## Different rows:
# from1
plot(0:max(t), p.dfE$p11, pch=20, ylim=c(0,1))
points(0:max(t), p.dfE$p12, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.dfE$p13, pch=20, ylim=c(0,1), col="red")
abline(h = c(0.2575, 0.3968, 0.3457))
points(p.msm.df$time, p.msm.df$t11, pch=20, col="darkblue")
points(p.msm.df$time, p.msm.df$t12, pch=20, ylim=c(0,1), col="orange")
points(p.msm.df$time, p.msm.df$t13, pch=20, ylim=c(0,1), col="orange")
abline(h = c(0.4545, 0.337, 0.2086))

# from2
plot(0:max(t), p.dfE$p22, pch=20, ylim=c(0,1))
points(0:max(t), p.dfE$p21, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.dfE$p23, pch=20, ylim=c(0,1), col="red")
points(p.msm.df$time, p.msm.df$t22, pch=20, col="darkblue")
points(p.msm.df$time, p.msm.df$t21, pch=20, ylim=c(0,1), col="orange")
points(p.msm.df$time, p.msm.df$t23, pch=20, ylim=c(0,1), col="orange")

# from3
plot(0:max(t), p.dfE$p33, pch=20, ylim=c(0,1))
points(0:max(t), p.dfE$p31, pch=20, ylim=c(0,1), col="red")
points(0:max(t), p.dfE$p32, pch=20, ylim=c(0,1), col="red")
points(p.msm.df$time, p.msm.df$t33, pch=20, col="darkblue")
points(p.msm.df$time, p.msm.df$t31, pch=20, ylim=c(0,1), col="orange")
points(p.msm.df$time, p.msm.df$t32, pch=20, ylim=c(0,1), col="orange")

######### Compare exponential output:
covariates <- list(c(CategoricalCovariate="Level1"), c(CategoricalCovariate="Level2"), c(CategoricalCovariate="Level3"))
p <- lapply(covariates, FUN=function(c) lapply(1:max(t), FUN = function(t) pmatrix.msm(output.msm.cov.CG, t, covariates = as.list(c))))

# Convert into dataframe
p.vector <- lapply(p, FUN = function(x) lapply(x, FUN = function(z) as.vector(t(z))))
p.msm.df <- lapply(1:3, function(x) data.frame(matrix(ncol = (length(p.vector[[1]][[1]]) + 1), nrow = max(t))))
p.msm.df <- lapply(1:3, function(y) setNames(as.data.frame(cbind(1:max(t), sapply(1:(ncol(p.msm.df[[y]])-1), FUN=function(x) p.msm.df[[y]][,x+1] <- sapply(p.vector[[y]], FUN=function(z) z[x])))),
                                             c("time", "t11", "t12", "t13", "t21", "t22", "t23", "t31", "t32", "t33")))

# For cat 2, p(1,2) and p(3,2) are increased, p(2,3) is decreased
# For cat 3, p(1,2) and p(1,3) are decreased, p(2,3) is decreased
# Effect of covariates

plot(0:max(t), p.dfE1$p12, pch=20, ylim=c(0,1))
points(0:max(t), p.dfE2$p12, pch=20, col = "maroon")
points(0:max(t), p.dfE3$p12, pch=20, col = "darkgreen")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t12, pch=20, col = "grey")
points(p.msm.df[[2]]$time, p.msm.df[[2]]$t12, pch=20, col = "red")
points(p.msm.df[[3]]$time, p.msm.df[[3]]$t12, pch=20, col = "green")

plot(0:max(t), p.dfE1$p11, pch=20, ylim=c(0,1))
points(0:max(t), p.dfE2$p11, pch=20, col = "maroon")
points(0:max(t), p.dfE3$p11, pch=20, col = "darkgreen")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t11, pch=20, col = "grey")
points(p.msm.df[[2]]$time, p.msm.df[[2]]$t11, pch=20, col = "red")
points(p.msm.df[[3]]$time, p.msm.df[[3]]$t11, pch=20, col = "green")

plot(0:max(t), p.dfE1$p23, pch=20, ylim=c(0,1))
points(0:max(t), p.dfE2$p23, pch=20, col = "maroon")
points(0:max(t), p.dfE3$p23, pch=20, col = "darkgreen")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t23, pch=20, col = "grey")
points(p.msm.df[[2]]$time, p.msm.df[[2]]$t23, pch=20, col = "red")
points(p.msm.df[[3]]$time, p.msm.df[[3]]$t23, pch=20, col = "green")





pmatrix.msm(output.msm.base.CG, 1)
pmatrix.msm(output.msm.cov.CG, 1, covariates = list(CategoricalCovariate="Level1")) %>% round(2)
pmatrix.msm(output.msm.cov.CG, 1, covariates = list(CategoricalCovariate="Level2")) %>% round(2)
pmatrix.msm(output.msm.cov.CG, 1, covariates = list(CategoricalCovariate="Level3")) %>% round(2)






## Run a differential equation solver

# X: 1st continuous param, 2nd second cat, 3rd third cat
parameters$X <- c(0,0,0)
start1.X1 <- ode(y = P0A, times = t, func = dP_dt_full, parms = parameters)
start2.X1 <- ode(y = P0B, times = t, func = dP_dt_full, parms = parameters)
start3.X1 <- ode(y = P0C, times = t, func = dP_dt_full, parms = parameters)

parameters$X <- c(0,1,0)
start1.X2 <- ode(y = P0A, times = t, func = dP_dt_full, parms = parameters)
start2.X2 <- ode(y = P0B, times = t, func = dP_dt_full, parms = parameters)
start3.X2 <- ode(y = P0C, times = t, func = dP_dt_full, parms = parameters)

parameters$X <- c(0,0,1)
start1.X3 <- ode(y = P0A, times = t, func = dP_dt_full, parms = parameters)
start2.X3 <- ode(y = P0B, times = t, func = dP_dt_full, parms = parameters)
start3.X3 <- ode(y = P0C, times = t, func = dP_dt_full, parms = parameters)


# Check it out: INLA
plot(start1.X1[,1], start1.X1[,2], pch=20, ylim=c(0,1))
points(start1.X2[,1], start1.X2[,2], pch=20, ylim=c(0,1), col="red")
points(start1.X3[,1], start1.X3[,2], pch=20, ylim=c(0,1), col="darkgreen")

plot(p.msm.df[[1]]$time, p.msm.df[[1]]$t11, pch=20, ylim=c(0,1))
points(p.msm.df[[2]]$time, p.msm.df[[2]]$t11, pch=20, ylim=c(0,1), col="red")
points(p.msm.df[[3]]$time, p.msm.df[[3]]$t11, pch=20, ylim=c(0,1), col="darkgreen")

plot(start1.X1[,1], start1.X1[,3], pch=20, ylim=c(0,1))
points(start1.X2[,1], start1.X2[,3], pch=20, ylim=c(0,1), col="red")
points(start1.X3[,1], start1.X3[,3], pch=20, ylim=c(0,1), col="darkgreen")

plot(p.msm.df[[1]]$time, p.msm.df[[1]]$t12, pch=20, ylim=c(0,1))
points(p.msm.df[[2]]$time, p.msm.df[[2]]$t12, pch=20, ylim=c(0,1), col="red")
points(p.msm.df[[3]]$time, p.msm.df[[3]]$t12, pch=20, ylim=c(0,1), col="darkgreen")

plot(start1.X1[,1], start1.X1[,2], pch=20, ylim=c(0,1))
points(start1.X0[,1], start1.X0[,2], pch=20, ylim=c(0,1), col="red")



plot(start1[,1], start1[,2], pch=20, ylim=c(0,1))
points(start1[,1], start1[,3], pch=20, col="red")
points(start1[,1], start1[,4], pch=20, col="darkgreen")

plot(start2[,1], start2[,3], pch=20, ylim=c(0,1))
points(start2[,1], start2[,2], pch=20, col="red")
points(start2[,1], start2[,4], pch=20, col="darkgreen")


# msm
plot(p.msm.df[[1]]$time, p.msm.df[[1]]$t11, pch=20, ylim=c(0,1))
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t12, pch=20, col="red")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t13, pch=20, col="darkgreen")

plot(p.msm.df[[1]]$time, p.msm.df[[1]]$t22, pch=20, ylim=c(0,1))
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t21, pch=20, col="red")
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t31, pch=20, col="darkgreen")

# Compare
plot(start1[,1], start1[,2], pch=20, ylim=c(0,1))
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t11, pch=20, ylim=c(0,1), col="red")

plot(start2[,1], start2[,2], pch=20, ylim=c(0,1))
points(p.msm.df[[1]]$time, p.msm.df[[1]]$t21, pch=20, ylim=c(0,1), col="red")


apply(p.msm.df[[1]] %>% select(t11,t12,t13), 1, sum)



weib.surv$dic$dic
weib.surv$waic$waic

weib.surv.base$dic$dic
weib.surv.base$waic$waic

exp.surv.base$dic$dic
exp.surv.base$waic$waic

###### OLD #######
# Let's plot the covariate coefficients:
tr <- c("11", "12", "13", "21", "22", "23", "31", "32", "33")
coef.real <- c(as.vector(transition_matrix1), as.vector(transition_matrix2), as.vector(transition_matrix3))
input.coef <- data.frame(model = rep("real", length(coef.real)),
                         trans = rep(tr, 3),
                         cov = rep(c("1", "2", "3"), each=9),
                         prob = coef.real,
                         ci_u = rep(NA, length(coef.real)),
                         ci_l = rep(NA, length(coef.real)))

coef.msm <- c(as.vector(hazard.msm(output.msm.cov.CG)$ContinuousCovariate1[,1]),
              as.vector(hazard.msm(output.msm.cov.CG)$CategoricalCovariateLevel2[,1]),
              as.vector(hazard.msm(output.msm.cov.CG)$CategoricalCovariateLevel3[,1]))
ci_l.msm <- c(as.vector(hazard.msm(output.msm.cov.CG)$ContinuousCovariate1[,2]),
              as.vector(hazard.msm(output.msm.cov.CG)$CategoricalCovariateLevel2[,2]),
              as.vector(hazard.msm(output.msm.cov.CG)$CategoricalCovariateLevel3[,2]))
ci_u.msm <- c(as.vector(hazard.msm(output.msm.cov.CG)$ContinuousCovariate1[,3]),
              as.vector(hazard.msm(output.msm.cov.CG)$CategoricalCovariateLevel2[,3]),
              as.vector(hazard.msm(output.msm.cov.CG)$CategoricalCovariateLevel3[,3]))



####





# Recover values from intercept
exp(weib.surv$summary.fixed[rownames(weib.surv$summary.fixed)=="Intercept_S1",1])
exp(weib.surv$summary.fixed[rownames(weib.surv$summary.fixed)=="Intercept_S2",1])
exp(weib.surv$summary.fixed[rownames(weib.surv$summary.fixed)=="Intercept_S6",1])



tr <- c("11", "12", "13", "21", "22", "23", "31", "32", "33")
coef.real <- c(as.vector(transition_matrix1), as.vector(transition_matrix2), as.vector(transition_matrix3))
input.coef <- data.frame(model=rep("real", length(coef.real)),
                         trans=rep(tr, 3),
                         cov=rep(c("1", "2", "3"), each=9),
                         coef=coef.real,
                         ci_u=rep(NA, length(coef.real)),
                         ci_l=rep(NA, length(coef.real)))

coef.msm <- c(as.vector(hazard.msm(output.msm.cov.CG)$ContinuousCovariate1[,1]),
              as.vector(hazard.msm(output.msm.cov.CG)$CategoricalCovariateLevel2[,1]),
              as.vector(hazard.msm(output.msm.cov.CG)$CategoricalCovariateLevel3[,1]))
ci_l.msm <- c(as.vector(hazard.msm(output.msm.cov.CG)$ContinuousCovariate1[,2]),
              as.vector(hazard.msm(output.msm.cov.CG)$CategoricalCovariateLevel2[,2]),
              as.vector(hazard.msm(output.msm.cov.CG)$CategoricalCovariateLevel3[,2]))
ci_u.msm <- c(as.vector(hazard.msm(output.msm.cov.CG)$ContinuousCovariate1[,3]),
              as.vector(hazard.msm(output.msm.cov.CG)$CategoricalCovariateLevel2[,3]),
              as.vector(hazard.msm(output.msm.cov.CG)$CategoricalCovariateLevel3[,3]))

coef.inlaweib <- c()
coef.inlaexp <- c()

a <- qmatrix.msm(output.msm.cov.CG)
as.vector(a$estimates)

output.coef <- data.frame(model=c(rep("msm", length(coef.msm)), 
                                  rep("inla.weib", length(coef.inlaweib)),
                                  rep("inla.exp", length(coef.inlaexp))),
                          trans=rep(tr, 9),
                          cov=rep(rep(c("1", "2", "3"), each=9),3),
                          coef=c(coef.msm, coef.inlaweib, coef.inlaexp),
                          ci_u=rep(NA, length(coef.real)),
                          ci_l=rep(NA, length(coef.real)))

hazard.msm(output.msm.cov.CG)                          


coef <- rbind(input.coef, output.coef)



m <- as.data.frame(cbind(1:60, sapply(1:(ncol(p.msm.df[[1]])-1), FUN=function(x) p.msm.df[[1]][,x+1] <- sapply(p.vector[[1]], FUN=function(z) z[x]))))

pv <- p.vector[[1]]




######## OLD #######

#### Analysis ####
### msm
#Check with dataset existing:
data <- read_csv(here("Data", "BTE", paste0("bte_", "4bM", "_msm_ready.csv")),
                 col_types = cols(.default = col_guess(),
                                  sp_class = col_integer(),
                                  cov_pert_class = col_factor(),
                                  cov_pert_sev = col_factor(),
                                  cov_time_pert = col_double()))[,-1]

# Rework df_select to make it runnable by msm:
data_msm <- df_select %>% 
  pivot_longer(cols = starts_with("Step"),
               names_to = "Time",
               names_prefix = "Step_",
               values_to = "State") %>%
  mutate(Time=(as.numeric(Time)-1)) %>% 
  mutate(ID=rep(1:9999, each=5)) %>% 
  arrange(ID, as.numeric(Time)) %>% 
  select(ID, Time, State, ContinuousCovariate1, CategoricalCovariate)

# Create statetable
msm_state <- statetable.msm(State, ID, data=data_msm)
round(funrar::make_relative(msm_state), 3)

# Define the model (all are possible here)
Q.model <- rep(0.3, 9)
dim(Q.model) <- c(3,3)

# Get initial estimate for Q

Q.init <- crudeinits.msm(State ~ Time, ID, data=data_msm, qmatrix=Q.model)

# Run base msm

output.msm.base.CG <-  msm(State ~ Time, ID, data=data_msm, qmatrix = Q.init, method= "CG",
                           control = list(fnscale = 100000))
output.msm.base.CG
output.msm.base.NM <-  msm(State ~ Time, ID, data=data_msm, qmatrix = Q.init, method= "Nelder-Mead",
                           control = list(fnscale = 100000))
output.msm.base.NM
output.msm.base.BFGS <-  msm(State ~ Time, ID, data=data_msm, qmatrix = Q.init, method= "L-BFGS-B",
                            control = list(fnscale = 100000))
output.msm.base.BFGS


qmatrix.msm(output.msm.base.CG)

pmatrix.msm(output.msm.base.CG, 1)
pmatrix.msm(output.msm.base.CG, 10)

pnext.msm(output.msm.base.CG)

# Run msm with covariates
output.msm.cov.CG <-  msm(State ~ Time, ID, data=data_msm, qmatrix = Q.init, method= "CG",
                           control = list(fnscale = 100000), covariates = ~ CategoricalCovariate + ContinuousCovariate1)
output.msm.cov.CG

output.msm.cov.CG2 <-  msm(State ~ Time, ID, data=data_msm, qmatrix = Q.init, method= "CG",
                          control = list(fnscale = 10000), covariates = ~ CategoricalCovariate + ContinuousCovariate1)
output.msm.cov.CG2


output.msm.cov.NM <-  msm(State ~ Time, ID, data=data_msm, qmatrix = Q.init, method= "Nelder-Mead",
                          control = list(fnscale = 100000), covariates = ~ CategoricalCovariate + ContinuousCovariate1)
output.msm.cov.NM

output.msm.cov.BFGS <-  msm(State ~ Time, ID, data=data_msm, qmatrix = Q.init, method= "L-BFGS-B",
                          control = list(fnscale = 100000), covariates = ~ CategoricalCovariate + ContinuousCovariate1)
output.msm.cov.BFGS

qmatrix.msm(output.msm.cov.CG)
pmatrix.msm(output.msm.cov.CG, 10, ci="normal")
pmatrix.msm(output.msm.cov.CG, 1, covariates = list(CategoricalCovariate="Level1"))
pmatrix.msm(output.msm.cov.CG, 1, covariates = list(CategoricalCovariate="Level2"))
pmatrix.msm(output.msm.cov.CG, 1, covariates = list(CategoricalCovariate="Level3"))

qmatrix.msm(output.msm.cov.CG, covariates = list(CategoricalCovariate="Level3"))
qmatrix.msm(output.msm.cov.CG, covariates = list(CategoricalCovariate="Level2"))
qmatrix.msm(output.msm.cov.CG, covariates = list(CategoricalCovariate="Level1"))
pmatrix.msm(output.msm.cov.CG, 1, covariates = list(CategoricalCovariate="Level1"))
pmatrix.msm(output.msm.cov.CG, 1, covariates = list(CategoricalCovariate="Level2"))
pmatrix.msm(output.msm.cov.CG, 1, covariates = list(CategoricalCovariate="Level3"))

HR.msm.cov.CG <- hazard.msm(output.msm.cov.CG)
HR.msm.cov.CG
##Create some plots p(t)
#Get a list of pmatrices

covariates <- list(c(CategoricalCovariate="Level1"), c(CategoricalCovariate="Level2"), c(CategoricalCovariate="Level3"))
p <- lapply(covariates, FUN=function(c) lapply(1:50, FUN = function(t) pmatrix.msm(output.msm.cov.CG, t, covariates = as.list(c))))

# Convert into dataframe
p.vector <- lapply(p, FUN = function(x) lapply(x, as.vector))
p.msm.df <- lapply(1:3, function(x) data.frame(matrix(ncol = (length(p.vector[[1]][[1]]) + 1), nrow = 50)))
p.msm.df <- lapply(1:3, function(y) setNames(as.data.frame(cbind(1:50, sapply(1:(ncol(p.msm.df[[y]])-1), FUN=function(x) p.msm.df[[y]][,x+1] <- sapply(p.vector[[y]], FUN=function(z) z[x])))),
                                         c("time", "t11", "t12", "t13", "t21", "t22", "t23", "t31", "t32", "t33")))

plot(p.msm.df[[1]]$time, p.msm.df[[1]]$t12)
points(p.msm.df[[2]]$time, p.msm.df[[2]]$t12, col='red')
points(p.msm.df[[3]]$time, p.msm.df[[3]]$t12, col='blue')

# Make the plots

plots.msm.p <- function(p.df, tr.name){
  cbind(p.df[[1]]$time, p.df[[1]][tr.name], p.df[[2]][tr.name], p.df[[3]][tr.name]) %>% 
    as.data.frame() %>% 
    `colnames<-`(c("time", "level1", "level2", "level3")) %>% 
    pivot_longer(cols = c(level1, level2, level3), names_to = "Covariate") %>% 
    ggplot +
    geom_point(aes(x=time, y=value, color=Covariate))
}

plots.msm.p(p.msm.df, "t12")

# For cat 2, p(1,2) and p(3,2) are increased, p(2,3) is decreased
# For cat 3, p(1,2) and p(1,3) are decreased, p(2,3) is decreased

###INLA #####
# Fix data:
data.inla.inter <- data_msm |>
  select(ID, Time, State, CategoricalCovariate, ContinuousCovariate1) |>
  mutate(from=State) |>
  mutate(to=State) |>
  mutate(entry=Time) |>
  relocate(Time, .after = last_col())

ph.from <- c(0, data.inla.inter$to)
data.inla.inter$from <- ph.from[1:(length(ph.from)-1)]

ph.entry <- c(0, data.inla.inter$Time)
data.inla.inter$entry <- ph.entry[1:(length(ph.entry)-1)] 

data.inla <- data.inla.inter |>
  filter(Time>0) |>
  select(-State) |>
  mutate(trans=to-from) |>
  filter(trans != 0) |> 
  ungroup() |>
  select(-ID, -trans)

##State table:
st <- matrix(0, ncol=3, nrow=3)

for(i in 1:nrow(data.inla)) {
  st[data.inla$from[i], data.inla$to[i]] <- st[data.inla$from[i], data.inla$to[i]] + 1
}

##Prepare the survival objects
#First create the empty event dataframes:

event.list <- vector("list", 6)
if(T) {
  for(i in 1:6) {
    event.list[[i]] <- data.frame(matrix(ncol = ncol(data.inla) + 1))
    colnames(event.list[[i]]) <- c(colnames(data.inla), "status")
    event.list[[i]] <- event.list[[i]][-1,]
  }
  
  # Now fill in the full event dataframes for each transition
  for(i in 1:6) { 
    print(i)
    data.tr <- data.inla |>
      filter(from == c(1,1,2,2,3,3)[i])
    
    tr.event <- vector("list", length = nrow(data.tr))
    for(j in 1:nrow(data.tr)) { 
      if(data.tr$to[j] == c(2,3,1,3,1,2)[i]) {
        tr.event[[j]] <- c(data.tr[j,], 1)
        names(tr.event[[j]])[length(tr.event[[j]])] <- "status"
      } else {
        tr.event[[j]] <- c(data.tr[j,], 0)
        tr.event[[j]]$to <- c(2,3,1,3,1,2)[i]
        names(tr.event[[j]])[length(tr.event[[j]])] <- "status"
      }
    }
    event.list[[i]] <- bind_rows(tr.event)
  }
}

Surv.list <- vector("list", 6)
for(i in 1:6) { 
  Surv.list[[i]] <- inla.surv(time = event.list[[i]]$Time, event = event.list[[i]]$status)
}
?inla.surv
s12 <- Surv.list[[1]]
s13 <- Surv.list[[2]]
s21 <- Surv.list[[3]]
s23 <- Surv.list[[4]]
s31 <- Surv.list[[5]]
s32 <- Surv.list[[6]]

inla.weib <- joint(formSurv=list(s12 ~ ContinuousCovariate1+CategoricalCovariate,
                              s13 ~ ContinuousCovariate1+CategoricalCovariate,
                              s21 ~ ContinuousCovariate1+CategoricalCovariate,
                              s23 ~ ContinuousCovariate1+CategoricalCovariate,
                              s31 ~ ContinuousCovariate1+CategoricalCovariate,
                              s32 ~ ContinuousCovariate1+CategoricalCovariate),
                basRisk = rep("weibullsurv", 6),
                dataSurv = list(event.list[[1]],
                                event.list[[2]],
                                event.list[[3]],
                                event.list[[4]],
                                event.list[[5]],
                                event.list[[6]]),
                control = list(config=TRUE))

inla.exp <- joint(formSurv=list(s12 ~ ContinuousCovariate1+CategoricalCovariate,
                                 s13 ~ ContinuousCovariate1+CategoricalCovariate,
                                 s21 ~ ContinuousCovariate1+CategoricalCovariate,
                                 s23 ~ ContinuousCovariate1+CategoricalCovariate,
                                 s31 ~ ContinuousCovariate1+CategoricalCovariate,
                                 s32 ~ ContinuousCovariate1+CategoricalCovariate),
                   basRisk = rep("exponentialsurv", 6),
                   dataSurv = list(event.list[[1]],
                                   event.list[[2]],
                                   event.list[[3]],
                                   event.list[[4]],
                                   event.list[[5]],
                                   event.list[[6]]),
                  control = list(config=TRUE))

summary(inla.weib)
summary(inla.exp)
log(HR.msm.cov.CG$CategoricalCovariateLevel2)
inla.exp$summary.fixed$mean
inla.weib$summary.fixed$mean

plot(inla.test)$Baseline
plot(inla.test)$Outcomes

#REcover values from intercept
exp(inla.weib$summary.fixed[rownames(inla.weib$summary.fixed)=="Intercept_S1",1])
exp(inla.weib$summary.fixed[rownames(inla.weib$summary.fixed)=="Intercept_S2",1])
exp(inla.weib$summary.fixed[rownames(inla.weib$summary.fixed)=="Intercept_S6",1])


#### Analysis -----
## plot covariate coefficients:
#TODO
tr <- c("11", "12", "13", "21", "22", "23", "31", "32", "33")
coef.real <- c(as.vector(transition_matrix1), as.vector(transition_matrix2), as.vector(transition_matrix3))
input.coef <- data.frame(model=rep("real", length(coef.real)),
                         trans=rep(tr, 3),
                         cov=rep(c("1", "2", "3"), each=9),
                         coef=coef.real,
                         ci_u=rep(NA, length(coef.real)),
                         ci_l=rep(NA, length(coef.real)))

coef.msm <- c(as.vector(hazard.msm(output.msm.cov.CG)$ContinuousCovariate1[,1]),
              as.vector(hazard.msm(output.msm.cov.CG)$CategoricalCovariateLevel2[,1]),
              as.vector(hazard.msm(output.msm.cov.CG)$CategoricalCovariateLevel3[,1]))
ci_l.msm <- c(as.vector(hazard.msm(output.msm.cov.CG)$ContinuousCovariate1[,2]),
              as.vector(hazard.msm(output.msm.cov.CG)$CategoricalCovariateLevel2[,2]),
              as.vector(hazard.msm(output.msm.cov.CG)$CategoricalCovariateLevel3[,2]))
ci_u.msm <- c(as.vector(hazard.msm(output.msm.cov.CG)$ContinuousCovariate1[,3]),
              as.vector(hazard.msm(output.msm.cov.CG)$CategoricalCovariateLevel2[,3]),
              as.vector(hazard.msm(output.msm.cov.CG)$CategoricalCovariateLevel3[,3]))

coef.inlaweib <- c()
coef.inlaexp <- c()

a <- qmatrix.msm(output.msm.cov.CG)
as.vector(a$estimates)

output.coef <- data.frame(model=c(rep("msm", length(coef.msm)), 
                                  rep("inla.weib", length(coef.inlaweib)),
                                  rep("inla.exp", length(coef.inlaexp))),
                          trans=rep(tr, 9),
                          cov=rep(rep(c("1", "2", "3"), each=9),3),
                          coef=c(coef.msm, coef.inlaweib, coef.inlaexp),
                          ci_u=rep(NA, length(coef.real)),
                          ci_l=rep(NA, length(coef.real)))

hazard.msm(output.msm.cov.CG)                          


coef <- rbind(input.coef, output.coef)

### Try the first method: Alvares equations:
## Calculate the hazard functions:
hW <- function(t, a, l){
  a * t**(a-1) * l
}
t <- 1:50

h12w <- hW(t, a=inla.weib$summary.hyperpar$mean[1], 
          l=exp(inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "Intercept_S1"]))
h13w <- hW(t, a=inla.weib$summary.hyperpar$mean[2], 
          l=exp(inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "Intercept_S2"]))
h21w <- hW(t, a=inla.weib$summary.hyperpar$mean[3], 
          l=exp(inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "Intercept_S3"]))
h23w <- hW(t, a=inla.weib$summary.hyperpar$mean[4], 
          l=exp(inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "Intercept_S4"]))
h31w <- hW(t, a=inla.weib$summary.hyperpar$mean[5], 
          l=exp(inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "Intercept_S5"]))
h32w <- hW(t, a=inla.weib$summary.hyperpar$mean[6], 
          l=exp(inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "Intercept_S6"]))

h12w.2 <- hW(t=t, a=inla.weib$summary.hyperpar$mean[1], 
          l=exp(inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "Intercept_S1"] +
                  inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "CategoricalCovariateLevel2_S1"]))
h13w.2 <- hW(t=t, a=inla.weib$summary.hyperpar$mean[2], 
          l=exp(inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "Intercept_S2"] +
                  inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "CategoricalCovariateLevel2_S2"]))
h21w.2 <- hW(t=t, a=inla.weib$summary.hyperpar$mean[3], 
          l=exp(inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "Intercept_S3"] +
                  inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "CategoricalCovariateLevel2_S3"]))
h23w.2 <- hW(t=t, a=inla.weib$summary.hyperpar$mean[4], 
          l=exp(inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "Intercept_S4"] +
                  inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "CategoricalCovariateLevel2_S4"]))
h31w.2 <- hW(t=t, a=inla.weib$summary.hyperpar$mean[5], 
          l=exp(inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "Intercept_S5"] +
                  inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "CategoricalCovariateLevel2_S5"]))
h32w.2 <- hW(t=t, a=inla.weib$summary.hyperpar$mean[6], 
          l=exp(inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "Intercept_S6"] +
                  inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "CategoricalCovariateLevel2_S6"]))

h12w.3 <- hW(t=t, a=inla.weib$summary.hyperpar$mean[1], 
            l=exp(inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "Intercept_S1"] +
                    inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "CategoricalCovariateLevel3_S1"]))
h13w.3 <- hW(t=t, a=inla.weib$summary.hyperpar$mean[2], 
            l=exp(inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "Intercept_S2"] +
                    inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "CategoricalCovariateLevel3_S2"]))
h21w.3 <- hW(t=t, a=inla.weib$summary.hyperpar$mean[3], 
            l=exp(inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "Intercept_S3"] +
                    inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "CategoricalCovariateLevel3_S3"]))
h23w.3 <- hW(t=t, a=inla.weib$summary.hyperpar$mean[4], 
            l=exp(inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "Intercept_S4"] +
                    inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "CategoricalCovariateLevel3_S4"]))
h31w.3 <- hW(t=t, a=inla.weib$summary.hyperpar$mean[5], 
            l=exp(inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "Intercept_S5"] +
                    inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "CategoricalCovariateLevel3_S5"]))
h32w.3 <- hW(t=t, a=inla.weib$summary.hyperpar$mean[6], 
            l=exp(inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "Intercept_S6"] +
                    inla.weib$summary.fixed$mean[rownames(inla.weib$summary.fixed) == "CategoricalCovariateLevel3_S6"]))

plot(1:50, h12w, pch=19)
points(1:50, h13w, col="red", pch=19)
points(1:50, h21w, col="red", pch=19)
points(1:50, h23w, col="red", pch=19)
points(1:50, h31w, col="red", pch=19)
points(1:50, h32w, col="red", pch=19)

plot(1:50, h12w, pch=19)
points(1:50, h12w.2, pch=19, col="red")
points(1:50, h12w.3, pch=19, col="green")

# For cat 2, p(1,2) and p(3,2) are increased, p(2,3) is decreased
# For cat 3, p(1,2) and p(1,3) are decreased, p(2,3) is decreased

## Convert to probabilities:
#First using cumsum like in the paper
#p(i,i)=
#p(i,j)
p11.cs <- exp(-cumsum(h12w)-cumsum(h13w))
p22.cs <- exp(-cumsum(h21w)-cumsum(h23w))
p33.cs <- exp(-cumsum(h31w)-cumsum(h32w))
p12.cs <- cumsum(p11.cs*h12w*p22.cs)
p13.cs <- cumsum(p11.cs*h13w*p33.cs)
p21.cs <- cumsum(p22.cs*h21w*p11.cs)
p23.cs <- cumsum(p22.cs*h23w*p33.cs)
p31.cs <- cumsum(p33.cs*h31w*p11.cs)
p32.cs <- cumsum(p33.cs*h32w*p22.cs)

p11.cs+p12.cs+p13.cs
plot(t, p11.cs+p12.cs+p13.cs)
p21.cs+p22.cs+p23.cs

plot(t, p11.cs, pch=19, ylim=c(0,1))
points(t, p12.cs, pch=19, col="red")
points(t, p13.cs, pch=19, col="red")

#Now with own formulas
p11.cs <- exp(-cumsum(h12w)-cumsum(h13w))
p22.cs <- exp(-cumsum(h21w)-cumsum(h23w))
p33.cs <- exp(-cumsum(h31w)-cumsum(h32w))
p12.cs <- 1-exp(-cumsum(h12w))
p13.cs <- 1-exp(-cumsum(h13w))
p21.cs <- 1-exp(-cumsum(h21w))
p23.cs <- 1-exp(-cumsum(h23w))
p31.cs <- 1-exp(-cumsum(h31w))
p32.cs <- 1-exp(-cumsum(h32w))

p11.cs+p12.cs+p13.cs
plot(t, p11.cs+p12.cs+p13.cs)
p21.cs+p22.cs+p23.cs

plot(t, p11.cs, pch=19, ylim=c(0,1))
points(t, p12.cs, pch=19, col="red")
points(t, p13.cs, pch=19, col="red")

### Other method: Differential equations:
## Try for the exponential first:
h12 <- exp(inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "Intercept_S1"])
h13 <- exp(inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "Intercept_S2"])
h21 <- exp(inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "Intercept_S3"])
h23 <- exp(inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "Intercept_S4"])
h31 <- exp(inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "Intercept_S5"])
h32 <- exp(inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "Intercept_S6"])

plot(t, h12w)
points(t, rep(h12, 50), col="red")

h12.2 <- exp(inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "Intercept_S1"] +
               inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "CategoricalCovariateLevel2_S1"])
h13.2 <- exp(inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "Intercept_S2"] +
               inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "CategoricalCovariateLevel2_S2"])
h21.2 <- exp(inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "Intercept_S3"] +
               inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "CategoricalCovariateLevel2_S3"])
h23.2 <- exp(inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "Intercept_S4"] +
               inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "CategoricalCovariateLevel2_S4"])
h31.2 <- exp(inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "Intercept_S5"] +
               inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "CategoricalCovariateLevel2_S5"])
h32.2 <- exp(inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "Intercept_S6"] +
               inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "CategoricalCovariateLevel2_S6"])

h12.3 <- exp(inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "Intercept_S1"] +
               inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "CategoricalCovariateLevel3_S1"])
h13.3 <- exp(inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "Intercept_S2"] +
               inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "CategoricalCovariateLevel3_S2"])
h21.3 <- exp(inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "Intercept_S3"] +
               inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "CategoricalCovariateLevel3_S3"])
h23.3 <- exp(inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "Intercept_S4"] +
               inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "CategoricalCovariateLevel3_S4"])
h31.3 <- exp(inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "Intercept_S5"] +
               inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "CategoricalCovariateLevel3_S5"])
h32.3 <- exp(inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "Intercept_S6"] +
               inla.exp$summary.fixed$mean[rownames(inla.exp$summary.fixed) == "CategoricalCovariateLevel3_S6"])

plot(t, rep(h12, length(t)), ylim=c(0,0.1))
points(t, rep(h12.2, length(t)), col="red")
points(t, rep(h12.3, length(t)), col="green")

# Create Q matrix
Q.inla <- matrix(c(-(h12+h13), h12, h13, h21, -(h21+h23), h23, h31, h32, -c(h31+h32)), nrow=3, ncol=3, byrow = T)
Q.inla2 <- matrix(c(-(h12.2+h13.2), h12.2, h13.2, h21.2, -(h21.2+h23.2), h23.2, h31.2, h32.2, -c(h31.2+h32.2)), nrow=3, ncol=3, byrow = T)
Q.inla3 <- matrix(c(-(h12.3+h13.3), h12.3, h13.3, h21.3, -(h21.3+h23.3), h23.3, h31.3, h32.3, -c(h31.3+h32.3)), nrow=3, ncol=3, byrow = T)

# For cat 2, p(1,2) and p(3,2) are increased, p(2,3) is decreased
# For cat 3, p(1,2) and p(1,3) are decreased, p(2,3) is decreased

# Load the deSolve package
library(deSolve)
P0 <- c(A = 1, B = 0, C = 0)

# Define set of equations (function needed for ode)
dP_dt <- function(t, P, Q) {
  A <- P[1]
  B <- P[2]
  C <- P[3]
  
  dA <- - (Q[1,2] + Q[1,3]) * A + Q[2,1] * B + Q[3,1] * C
  dB <- Q[1,2] * A - (Q[2,1] + Q[2,3]) * B + Q[3,2] * C
  dC <- Q[1,3] * A + Q[2,3] * B - (Q[3,1] + Q[3,2]) * C
  
  list(c(dA, dB, dC))
}

inla.odeA <- ode(y = P0, times = 0:50, func = dP_dt, parms = Q.inla)
inla.odeB <- ode(y = c(A = 0, B = 1, C = 0), times = 0:50, func = dP_dt, parms = Q.inla)
inla.odeC <- ode(y = c(A = 0, B = 0, C = 1), times = 0:50, func = dP_dt, parms = Q.inla)

p.inla.ode1 <- data.frame(time=0:50,
           t11=inla.odeA[,2], t12=inla.odeA[,3], t13=inla.odeA[,4],
           t21=inla.odeB[,2], t22=inla.odeB[,3], t23=inla.odeB[,4],
           t31=inla.odeC[,2], t32=inla.odeC[,3], t33=inla.odeC[,4])

inla.odeA <- ode(y = P0, times = 0:50, func = dP_dt, parms = Q.inla2)
inla.odeB <- ode(y = c(A = 0, B = 1, C = 0), times = 0:50, func = dP_dt, parms = Q.inla2)
inla.odeC <- ode(y = c(A = 0, B = 0, C = 1), times = 0:50, func = dP_dt, parms = Q.inla2)
p.inla.ode2 <- data.frame(time=0:50,
                          t11=inla.odeA[,2], t12=inla.odeA[,3], t13=inla.odeA[,4],
                          t21=inla.odeB[,2], t22=inla.odeB[,3], t23=inla.odeB[,4],
                          t31=inla.odeC[,2], t32=inla.odeC[,3], t33=inla.odeC[,4])

inla.odeA <- ode(y = P0, times = 0:50, func = dP_dt, parms = Q.inla3)
inla.odeB <- ode(y = c(A = 0, B = 1, C = 0), times = 0:50, func = dP_dt, parms = Q.inla3)
inla.odeC <- ode(y = c(A = 0, B = 0, C = 1), times = 0:50, func = dP_dt, parms = Q.inla3)
p.inla.ode3 <- data.frame(time=0:50,
                          t11=inla.odeA[,2], t12=inla.odeA[,3], t13=inla.odeA[,4],
                          t21=inla.odeB[,2], t22=inla.odeB[,3], t23=inla.odeB[,4],
                          t31=inla.odeC[,2], t32=inla.odeC[,3], t33=inla.odeC[,4])

p.inla.ode <- list(p.inla.ode1[-1,], p.inla.ode2[-1,], p.inla.ode3[-1,])

plots.msm.p(p.inla.ode, "t22")

## Matrix:
plot(t, inla.odeA[,2], ylim=c(0,1))
points(t, inla.odeA[,3], col="red")
points(t, inla.odeA[,4], col="red")

testerM <- matrix(c(qmatrix.msm(output.msm.cov.CG)[1,1][1], qmatrix.msm(output.msm.cov.CG)[1,2][1], qmatrix.msm(output.msm.cov.CG)[1,3][1],
                    qmatrix.msm(output.msm.cov.CG)[2,1][1], qmatrix.msm(output.msm.cov.CG)[2,2][1], qmatrix.msm(output.msm.cov.CG)[2,3][1],
                    qmatrix.msm(output.msm.cov.CG)[3,1][1], qmatrix.msm(output.msm.cov.CG)[3,2][1], qmatrix.msm(output.msm.cov.CG)[3,3][1]), ncol=3, nrow=3, byrow=T)

test <- expm(10*testerM)

P <- expm(1*Q.inla)
P2 <- expm(1*Q.inla2)
P3 <- expm(1*Q.inla3)
P <- expm(10*Q.inla)


p.matr <- list(lapply(1:50, FUN=function(x) expm(x*Q.inla)), lapply(1:50, FUN=function(x) expm(x*Q.inla2)), lapply(1:50, FUN=function(x) expm(x*Q.inla3)))

p.matr.vector <- lapply(p.matr, FUN = function(x) lapply(x, as.vector))
p.inla.matr <- lapply(1:3, function(x) data.frame(matrix(ncol = (length(p.matr.vector[[1]][[1]]) + 1), nrow = 50)))
p.inla.matr <- lapply(1:3, function(y) setNames(as.data.frame(cbind(1:50, sapply(1:(ncol(p.matr.df[[y]])-1), FUN=function(x) p.inla.matr[[y]][,x+1] <- sapply(p.matr.vector[[y]], FUN=function(z) z[x])))),
                                             c("time", "t11", "t12", "t13", "t21", "t22", "t23", "t31", "t32", "t33")))

plots.msm.p(p.df, "t33")
plots.msm.p(p.inla.ode, "t33")
plots.msm.p(p.inla.matr, "t33")

plots.msm.p(p.df, "t11")
plots.msm.p(p.inla.ode, "t11")
plots.msm.p(p.inla.matr, "t11")

plots.msm.p(p.df, "t22")
plots.msm.p(p.inla.ode, "t22")
plots.msm.p(p.inla.matr, "t22")

tr <- "t32"
plots.msm.p(p.df, tr)
plots.msm.p(p.inla.ode, tr)
plots.msm.p(p.inla.matr, tr)

### what does msm?

a <- qmatrix.msm(output.msm.cov.CG, covariates = list(CategoricalCovariate="Level2"))$estimates
expm(50*a)

### compare q matrices:
qmatrix.msm(output.msm.cov.CG, covariates = list(CategoricalCovariate="Level1"))$estimates
qmatrix.msm(output.msm.cov.CG, covariates = list(CategoricalCovariate="Level2"))$estimates
qmatrix.msm(output.msm.cov.CG, covariates = list(CategoricalCovariate="Level3"))$estimates

Q.inla
Q.inla2
Q.inla3

# For cat 2, p(1,2) and p(3,2) are increased, p(2,3) is decreased
# For cat 3, p(1,2) and p(1,3) are decreased, p(2,3) is decreased

#### Plot output probabilities ####
## plot covariate coefficients:

tr <- c("11", "12", "13", "21", "22", "23", "31", "32", "33")
prob.real <- c(as.vector(transition_matrix1), as.vector(transition_matrix2), as.vector(transition_matrix3))
input.prob <- data.frame(model=rep("real", length(prob.real)),
                         trans=rep(tr, 3),
                         cov=rep(c("1", "2", "3"), each=9),
                         prob=prob.real,
                         ci_u=rep(NA, length(prob.real)),
                         ci_l=rep(NA, length(prob.real)))

prob.msm <- c(as.vector(pmatrix.msm(output.msm.cov.CG, 1, covariates = list(CategoricalCovariate="Level1"))),
              as.vector(pmatrix.msm(output.msm.cov.CG, 1, covariates = list(CategoricalCovariate="Level2"))),
              as.vector(pmatrix.msm(output.msm.cov.CG, 1, covariates = list(CategoricalCovariate="Level3"))))
prob.alv <- c()
prob.diff <- c()
prob.matr <- c()

output.prob <- data.frame(model=c(rep("msm", length(prob.msm)), 
                                  rep("inla.weib", length(prob.inlaweib)),
                                  rep("inla.exp", length(prob.inlaexp))),
                          trans=rep(tr, 9),
                          cov=rep(rep(c("1", "2", "3"), each=9),3),
                          prob=c(prob.msm, prob.inlaweib, prob.inlaexp),
                          ci_u=rep(NA, length(prob.real)),
                          ci_l=rep(NA, length(prob.real)))

prob <- rbind(input.prob, output.prob)



########
hW(t = 1:50, a = inla.test$summary.hyperpar$mean[1], exp(inla.test$summary.fixed$mean[rownames(inla.test$summary.fixed)=="Intercept_S1"]))
hW(t = 1:50, a = inla.test$summary.hyperpar$mean[1], exp(inla.test$summary.fixed$mean[rownames(inla.test$summary.fixed)=="Intercept_S2"]))

hW(t = 1:50, a = inla.test$summary.hyperpar$mean[1], exp(inla.test$summary.fixed$mean[rownames(inla.test$summary.fixed)=="Intercept_S1"] + 
                                                           inla.test$summary.fixed$mean[rownames(inla.test$summary.fixed)=="CategoricalCovariateLevel2_S1"]))
hW(t = 1:50, a = inla.test$summary.hyperpar$mean[1], exp(inla.test$summary.fixed$mean[rownames(inla.test$summary.fixed)=="Intercept_S1"] + 
                                                           inla.test$summary.fixed$mean[rownames(inla.test$summary.fixed)=="CategoricalCovariateLevel3_S1"]))

rownames(inla.test$summary.fixed) <- paste0(rep(c("I", "L2", "L3"), 6), rep(as.character(1:6), each=3))

inla.test$summary.fixed <- as.data.frame(matrix(0, nrow=6, ncol=7))
covs <- data.frame(I=c("I", "I", "I"), C=c("L1", "L2", "L3"))

#H <- lapply(1:3, lapply(1:6, function(x) ))
H1 <- sapply(1:6, function(x) hW(t = 1:50, a = inla.test$summary.hyperpar$mean[x], 
                                exp(inla.test$summary.fixed$mean[i*3-2])))
H2 <- sapply(1:6, function(x) hW(t = 1:50, a = inla.test$summary.hyperpar$mean[x], 
                                 exp(inla.test$summary.fixed$mean[i*3-2]+
                                   inla.test$summary.fixed$mean[i*3-1])))
H3 <- sapply(1:6, function(x) hW(t = 1:50, a = inla.test$summary.hyperpar$mean[x], 
                                 exp(inla.test$summary.fixed$mean[i*3-2]+
                                       inla.test$summary.fixed$mean[i*3])))

inla.Q1 <- as.list(c(-H1[,1]-H1[,2], H1[,1], H1[,2],
                    H1[,3], -H1[,3]-H1[,4], H1[,4],
                    H1[,5], H1[,6], -H1[,5]-H1[,6]), nrow = 3, ncol = 3)

p11 <- exp(-cumsum(H1[,1])-cumsum(H1[,2]))
p22 <- exp(-cumsum(H1[,3])-cumsum(H1[,4]))
p33 <- exp(-cumsum(H1[,5])-cumsum(H1[,6]))
p11.2 <- exp(-cumsum(H2[,1])-cumsum(H2[,2]))
p22.2 <- exp(-cumsum(H2[,3])-cumsum(H2[,4]))
p33.2 <- exp(-cumsum(H2[,5])-cumsum(H2[,6]))
p11.3 <- exp(-cumsum(H3[,1])-cumsum(H3[,2]))
p22.3 <- exp(-cumsum(H3[,3])-cumsum(H3[,4]))
p33.3 <- exp(-cumsum(H3[,5])-cumsum(H3[,6]))
p12 <- cumsum(p11*H1[,1]*p22)
p12.2 <- cumsum(p11.2*H2[,1]*p22.2)
p12.3 <- cumsum(p11.3*H3[,1]*p22.3)

plot(1:50, p12)
points(1:50, p12.2, col='red')
points(1:50, p12.3, col='blue')

plots.msm.p(p.df, "t12")

plot(1:50, 1-exp(-sum()))
points(1:50, p12.2, col='red')
points(1:50, p12.3, col='blue')


#Try chatgpt

# Load the deSolve package
library(deSolve)

# Initial conditions
P0 <- c(A = 1, B = 0, C = 0)  # Starting in state A

# Time points
times <- seq(1, 50, by = 1)  # Adjust the time range and resolution as needed

# Define the system of differential equations
dP_dt <- function(t, P, H) {
  A <- P[1]
  B <- P[2]
  C <- P[3]
  
  dA <- - (H[t,1] + H[t,2]) * A + H[t,3] * B + H[t,5] * C
  dB <- H[t,1] * A - (H[t,3] + H[t,4]) * B + H[t,6] * C
  dC <- H[t,2] * A + H[t,4] * B - (H[t,5] + H[t,6]) * C
  
  list(c(dA, dB, dC))
}

# Solve the system of differential equations
out <- ode(y = P0, times = times, func = dP_dt, parms = H1)

# Extract probabilities
P_A <- out[, "A"]
P_B <- out[, "B"]
P_C <- out[, "C"]

plot(1:50, p.df[[1]]$t12)
points(1:50, P_B, col='red')

# Check that probabilities sum to 1 at all times
P_sum <- P_A + P_B + P_C
all.equal(P_sum, rep(1, length(P_sum)))  # Should be TRUE

# Plot the results
plot(times, P_A, type = "l", col = "blue", ylim = c(0, 1), ylab = "Probability", xlab = "Time", main = "State Probabilities Over Time")
lines(times, P_B, col = "red")
lines(times, P_C, col = "green")
legend("right", legend = c("State A", "State B", "State C"), col = c("blue", "red", "green"), lty = 1)

