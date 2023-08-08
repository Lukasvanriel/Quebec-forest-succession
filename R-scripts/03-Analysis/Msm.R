#### Runs the markov chain model using the msm package ####

### Load packages ###
library(tidyverse)
library(msm)
library(stringr)
library(here)

### Data ###
data <- read.csv(here("Data", "BTE", "bte_sp_class.csv"))[,-1]

### Prepare data to be compatible with msm required format ###

## Add time column that indicates time since first observation
#First remove rows that lack observation times
#TODO: Try to find way to deal with NA's in observation times
data <- data |>
  filter(! is.na(data$AN_PRO_SOU))

data_time <- data |>
  group_by(TESSELLE) |>
  mutate(time=AN_PRO_SOU-min(AN_PRO_SOU)) |>
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


###Add covariate data



#@#@#@#@
data_msm <- data_mult_filt

#write.csv(data_msm, here("Data", "BTE", "bte_msm_ready.csv"))

### Run msm ###

data_msm <- read.csv(here("Data", "BTE", "bte_msm_ready.csv"))[,-1]

## Create statetables:
msm_state <- statetable.msm(sp_class, TESSELLE, data=data_msm)
round(funrar::make_relative(msm_state), 3)

## Define models through Q matrices

# Multinomial analysis to determine impossible transitions

source(here("R-scripts", "03-Analysis", "multinomial.R"))

mnom <- multinom_model(data_msm)

Q.model <- as.matrix(round(mnom, 3)) #To determine which transitions are impossible (seems to be none)

## Get initial estimate for Q

Q.init  <- crudeinits.msm(sp_class ~ time, TESSELLE, data=data_msm, qmatrix=Q.model)

##Run msm

msm.reg <- msm( sp_class ~ time, subject=TESSELLE, data = data_msm, 
                qmatrix = Q.init) #numerical overflow in calculating likelihood

#Try some other things
msm.sc <- msm( sp_class ~ time, subject=TESSELLE, data = data_msm, 
                qmatrix = Q.init, control = list(fnscale = 5000000))

msm.cg <- msm( sp_class ~ time, subject=TESSELLE, data = data_msm, 
                  qmatrix = Q.init, method = "CG")

msm.sc.cg <- msm( sp_class ~ time, subject=TESSELLE, data = data_msm, 
               qmatrix = Q.init, control = list(fnscale = 5000000), method = "CG")

msm.nm <- msm( sp_class ~ time, subject=TESSELLE, data = data_msm, 
                  qmatrix = Q.init, method = "Nelder-Mead")

msm.sc.nm <- msm( sp_class ~ time, subject=TESSELLE, data = data_msm, 
                  qmatrix = Q.init, control = list(fnscale = 5000000), method = "Nelder-Mead")


