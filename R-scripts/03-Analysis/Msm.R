# Lukas Van Riel
# 2023-09-13
# Runs the markov chain model using the msm package
rm(list = ls())

### Load packages ####
library(tidyverse)
library(msm)
library(minqa)
library(stringr)
library(here)
library(conflicted)

conflicts_prefer(dplyr::filter)

### Data ####
if(F) {
data <- read.csv(here("Data", "BTE", "bte_cov_class.csv"))[,-1]
sum(is.na(data$cov_time_pert))
sum(data$cov_time_pert == -1)

### Prepare data to be compatible with msm required format ###

## Add time column that indicates time since first observation
#First remove rows that lack observation times

data_time <- data |>
  group_by(TESSELLE) |>
  mutate(time=year-min(year)) |>
  arrange(TESSELLE, time)

## Filter out Tesselle with only one measurement
single_meas <- data_time |>
  group_by(TESSELLE) |> 
  summarise(n=n()) |>
  filter(n==1)

data_mult <- subset(data_time, ! TESSELLE %in% single_meas$TESSELLE)

## Check for multiple observations at same time:
diff.previous <- list()
for (i in 2:nrow(data_mult)) {
  print(i)
  if(data_mult$time[i] == data_mult$time[i-1]) {
    diff.previous[[length(diff.previous) + 1]] <- i
  }
}

data_mult$unique <- ifelse(1:nrow(data_mult) %in% diff.previous , FALSE, TRUE)

data_mult_filt <- data_mult |>
  filter(unique) |>
  select(-unique)

write.csv(data_mult_filt, here("Data", "BTE", "bte_msm_ready.csv"))
}

data_msm <- read.csv(here("Data", "BTE", "bte_msm_ready.csv"))[,-1]

### Functions ####

## Function to more easily run msm remotely
##Analyse msm model output 
#TODO: extend this function so it can directly compare different runs
#Function to plot the parameters of a msm function run
#@param: Model has to msm list object as outputted by msm function
plot.msm <- function(model, path=NA) {
  # Extract output parameters from model
  output.params <- data.frame(est=model$estimates.t,
                              ci_l=model$ci[,1],
                              ci_h=model$ci[,2])
  
  # Add correct transitions to dataframe; e.g. 21 means from 2 to 1 
  dim <- nrow(model$qmodel$imatrix) # Extract the extent of the transition matrix
  
  index <- seq(1, dim, by=1)
  transitions <- c()
  for(i in index) {
    for(j in index[-i]) {
      if(model$qmodel$imatrix[i, j]) {transitions <- c(transitions, 10*i + j%%10)}
    }
  }
  output.params$trans <- transitions
  
  if(is.na(path)) {
    ggplot(output.params, aes(x=trans, y=est)) +
      geom_point() +
      geom_errorbar(aes(ymin=ci_l, ymax=ci_h)) +
      labs(title = "Parameter Estimates with Error Bars",
           x = "Transition",
           y = "Estimated q") +
      scale_x_continuous(breaks = seq(10, 100, by = 10)) +
      theme_minimal()
    
  } else {
  
    plot <- ggplot(output.params, aes(x=trans, y=est)) +
      geom_point() +
      geom_errorbar(aes(ymin=ci_l, ymax=ci_h)) +
      labs(title = "Parameter Estimates with Error Bars",
           x = "Transition",
           y = "Estimated q") +
      scale_x_continuous(breaks = seq(10, 100, by = 10)) +
      theme_minimal()
    ggsave(filename = path, plot = plot, device = "pdf")
    }
}

run_remote_msm <- function(data_msm, qmatrix, md = "BFGS", ctrl = 1, cov = NA, name.out.rds) {
  if(! all(c("sp_class", "time", "TESSELLE") %in% colnames(data_msm))) {
    return("Missing information in data.")
  }
  
  if(md %in% c("BFGS", "CG", "Nelder-Mead", "SANN")) {
    msm.model <- tryCatch(msm( sp_class ~ time, subject=TESSELLE, data = data_msm,
                               qmatrix = qmatrix, method= md, control = list(fnscale = ctrl)),
                          error = function(e) NA)
    
  } else if (md %in% c("nlm", "bobyqa", "Fisher")){
    msm.model <- tryCatch(msm( sp_class ~ time, subject=TESSELLE, data = data_msm,
                               qmatrix = qmatrix, opt.method= md, control = list(fnscale = ctrl)),
                          error = function(e) NA)
  }
  
  tryCatch(saveRDS(msm.model, here("Data-Output", "msm", name.out.rds)), error = function(e) NA)
  
  # Plot model parameters
  tryCatch(plot.msm(msm.model, here("Data-Output", "msm",
                                    paste0(str_sub(name.out.rds, end = -5), "_qvalues.pdf"))),
           error = function(e) NA)

}


### Initialise ####


## Create statetables:
msm_state <- statetable.msm(sp_class, TESSELLE, data=data_msm)
round(funrar::make_relative(msm_state), 3)

## Define models through Q matrices

# Multinomial analysis to determine impossible transitions

source(here("R-scripts", "03-Analysis", "multinomial.R"))

if (file.exists(here("Data-Output", "msm", "multinom.rds"))) {
  mnom <- readRDS(here("Data-Output", "msm", "multinom.rds"))
} else {
  mnom <- multinom_model(data_msm)
  saveRDS(mnom, here("Data-Output", "msm", "multinom.rds"))
}

#To determine which transitions are impossible
Q.model <- as.matrix(round(mnom, 3)) 

## Get initial estimate for Q

Q.init <- crudeinits.msm(sp_class ~ time, TESSELLE, data=data_msm, qmatrix=Q.model)


### Run msm ####

# run_remote_msm(data_msm = data_msm, qmatrix = Q.init, ctrl = 1, name.out.rds = "msm.reg.rds")
run_remote_msm(data_msm = data_msm, qmatrix = Q.init, ctrl = 5000000, name.out.rds = "msm.reg.sc5M.rds")

## Lets create nested loops to try out different things:

scaling <- c(5000000, 1)
methods <- c("BFGS", "CG", "Nelder-Mead", "SANN", "nlm", "bobyqa", "Fisher")


for(S in scaling) {
  for(M in methods) {
    run_remote_msm(data_msm = data_msm, qmatrix = Q.init, ctrl = S, md = M,
                   name.out.rds = paste0("msm.", M,".sc", as.character(S), ".rds"))
  }
}


### For future runs ####


##---
#msm.reg <- msm( sp_class ~ time, subject=TESSELLE, data = data_msm, 
#                qmatrix = Q.init) #numerical overflow in calculating likelihood

#Try some other things
# msm.sc <- msm( sp_class ~ time, subject=TESSELLE, data = data_msm, 
#                 qmatrix = Q.init, control = list(fnscale = 5000000))
# 
# 
# data_msm <- data_msm %>% mutate_at(c("cov_CMI", "cov_Tmean"), ~(scale(.) %>% as.vector))
# msm.sc <- msm( sp_class ~ time, subject=TESSELLE, data = data_msm, 
#                qmatrix = Q.init, control = list(fnscale = 5000000), covariates = ~ cov_Tmean)

#msm.cg <- msm( sp_class ~ time, subject=TESSELLE, data = data_msm, 
#                  qmatrix = Q.init, method = "CG")

#msm.sc.cg <- msm( sp_class ~ time, subject=TESSELLE, data = data_msm, 
#               qmatrix = Q.init, control = list(fnscale = 5000000), method = "CG")

#msm.nm <- msm( sp_class ~ time, subject=TESSELLE, data = data_msm, 
#                  qmatrix = Q.init, method = "Nelder-Mead")

#msm.sc.nm <- msm( sp_class ~ time, subject=TESSELLE, data = data_msm, 
#                  qmatrix = Q.init, control = list(fnscale = 5000000), method = "Nelder-Mead")







### Add covariates ###

# msm.sc <- msm( sp_class ~ time, subject=TESSELLE, data = data_msm_filt, 
#                qmatrix = Q.init, control = list(fnscale = 5000000)) #No global minimum

# msm.sc.cov <- msm( sp_class ~ time, subject=TESSELLE, data = data_msm_filt, 
#                qmatrix = Q.init, control = list(fnscale = 5000000),
#                covariates = list("4-6" = ~ cov_soil)) #Crashes Rstudio



if(F) {

### Try with subset of data ###

#Pick a certain percentage of TESSELLE at random
names.subs10 <- names(table(data_msm_filt$TESSELLE))[sample.int(length(table(data_msm_filt$TESSELLE)), 
                                                                0.1 * length(table(data_msm_filt$TESSELLE)), replace = FALSE)]
names.subs20 <- names(table(data_msm_filt$TESSELLE))[sample.int(length(table(data_msm_filt$TESSELLE)), 
                                                                0.2 * length(table(data_msm_filt$TESSELLE)), replace = FALSE)]
names.subs05 <- names(table(data_msm_filt$TESSELLE))[sample.int(length(table(data_msm_filt$TESSELLE)), 
                                                                0.05 * length(table(data_msm_filt$TESSELLE)), replace = FALSE)]
names.subs01 <- names(table(data_msm_filt$TESSELLE))[sample.int(length(table(data_msm_filt$TESSELLE)), 
                                                                0.01 * length(table(data_msm_filt$TESSELLE)), replace = FALSE)]

#Subset data
data_subs10 <- subset(data_msm_filt, TESSELLE %in% names.subs10)
data_subs20 <- subset(data_msm_filt, TESSELLE %in% names.subs20)
data_subs05 <- subset(data_msm_filt, TESSELLE %in% names.subs05)
data_subs01 <- subset(data_msm_filt, TESSELLE %in% names.subs01)

##
msm.sub10 <- msm( sp_class ~ time, subject=TESSELLE, data = data_subs10, 
               qmatrix = Q.init, control = list(fnscale = 5000000)) #Reported estimates are not the maximum likelihood.
msm.sub20 <- msm( sp_class ~ time, subject=TESSELLE, data = data_subs20, 
                  qmatrix = Q.init, control = list(fnscale = 5000000)) #Reported estimates are not the maximum likelihood.
msm.sub05 <- msm( sp_class ~ time, subject=TESSELLE, data = data_subs05, 
                  qmatrix = Q.init, control = list(fnscale = 5000000)) #Reported estimates are not the maximum likelihood.
msm.sub01 <- msm( sp_class ~ time, subject=TESSELLE, data = data_subs01, 
                  qmatrix = Q.init, control = list(fnscale = 5000000)) #Reported estimates are not the maximum likelihood.

plot.msm(msm.sc)
plot.msm(msm.sub10)
plot.msm(msm.sub20)
plot.msm(msm.sub05)
plot.msm(msm.sub01)

#plot
index <- seq(1, 9, by=1)
transitions <- c()
for(i in index) {
  for(j in index[-i]) {
    transitions <- c(transitions, 10*i + j%%10)
  }
}

df <- data.frame(trans=transitions, est=msm.sc$estimates.t, ci_l=msm.sc$ci[,1], ci_h=msm.sc$ci[,2],
                 v10=msm.sub10$estimates.t, v10_l=msm.sub10$ci[,1], v10_h=msm.sub10$ci[,2],
                 v20=msm.sub20$estimates.t, v20_l=msm.sub20$ci[,1], v20_h=msm.sub20$ci[,2],
                 v05=msm.sub05$estimates.t, v05_l=msm.sub05$ci[,1], v05_h=msm.sub05$ci[,2],
                 v01=msm.sub01$estimates.t, v01_l=msm.sub01$ci[,1], v01_h=msm.sub01$ci[,2])

ggplot(df) +
  geom_point(aes(x=trans-0.5, y=est), col="black") +
  geom_errorbar(aes(x=trans-0.5, y=est, ymin=ci_l, ymax=ci_h), col="black") +
  geom_point(aes(x=trans - 0.25, y=v20), col="blue") +
  geom_errorbar(aes(x=trans - 0.25, y=v20, ymin=v20_l, ymax=v20_h), col="blue") +
  geom_point(aes(x=trans + 0.0, y=v10), col="orange") +
  geom_errorbar(aes(x=trans + 0.0, y=v10, ymin=v10_l, ymax=v10_h), col="orange") +
  geom_point(aes(x=trans + 0.25, y=v05), col="red") +
  geom_errorbar(aes(x=trans + 0.25, y=v05, ymin=v05_l, ymax=v05_h), col="red") +
  geom_point(aes(x=trans + 0.5, y=v01), col="purple") +
  geom_errorbar(aes(x=trans + 0.55, y=v01, ymin=v01_l, ymax=v01_h), col="purple") +
  labs(x = "Transition",
       y = "Estimated q") +
  theme_minimal()

ggplot(df) +
  geom_point(aes(x=trans-0.5, y=est), col="black") +
  geom_errorbar(aes(x=trans-0.5, y=est, ymin=ci_l, ymax=ci_h), col="black") +
  geom_point(aes(x=trans - 0.25, y=v20), col="blue") +
  geom_errorbar(aes(x=trans - 0.25, y=v20, ymin=v20_l, ymax=v20_h), col="blue") +
  geom_point(aes(x=trans + 0.0, y=v10), col="orange") +
  geom_errorbar(aes(x=trans + 0.0, y=v10, ymin=v10_l, ymax=v10_h), col="orange") +
  geom_point(aes(x=trans + 0.25, y=v05), col="red") +
  geom_errorbar(aes(x=trans + 0.25, y=v05, ymin=v05_l, ymax=v05_h), col="red") +
  geom_point(aes(x=trans + 0.5, y=v01), col="purple") +
  geom_errorbar(aes(x=trans + 0.55, y=v01, ymin=v01_l, ymax=v01_h), col="purple") +
  labs(x = "Transition",
       y = "Estimated q") +
  xlim(c(11,25)) +
  theme_minimal()

ggplot(df) +
  geom_point(aes(x=trans-0.5, y=est), col="black") +
  geom_errorbar(aes(x=trans-0.5, y=est, ymin=ci_l, ymax=ci_h), col="black") +
  geom_point(aes(x=trans - 0.25, y=v20), col="blue") +
  geom_errorbar(aes(x=trans - 0.25, y=v20, ymin=v20_l, ymax=v20_h), col="blue") +
  geom_point(aes(x=trans + 0.0, y=v10), col="orange") +
  geom_errorbar(aes(x=trans + 0.0, y=v10, ymin=v10_l, ymax=v10_h), col="orange") +
  geom_point(aes(x=trans + 0.25, y=v05), col="red") +
  geom_errorbar(aes(x=trans + 0.25, y=v05, ymin=v05_l, ymax=v05_h), col="red") +
  geom_point(aes(x=trans + 0.5, y=v01), col="purple") +
  geom_errorbar(aes(x=trans + 0.55, y=v01, ymin=v01_l, ymax=v01_h), col="purple") +
  labs(x = "Transition",
       y = "Estimated q") +
  xlim(c(70,98)) +
  theme_minimal()

## Pmatrix for each
round(pmatrix.msm(msm.sc, t=10), 4)
round(pmatrix.msm(msm.sub10, t=10), 4)
round(pmatrix.msm(msm.sub20, t=10), 4)
round(pmatrix.msm(msm.sub05, t=10), 4)
round(pmatrix.msm(msm.sub01, t=10), 4)


##Now try subset with covariates
msm.sub05.cov <- msm( sp_class ~ time, subject=TESSELLE, data = data_subs05, 
                  qmatrix = Q.init, control = list(fnscale = 5000000),
                  covariates = list("4-6" = ~ cov_soil)) #

##

table(data_msm$cov_soil)

head(data_msm[data_msm$cov_soil=="  ",])
}
