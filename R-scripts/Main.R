# Lukas Van Riel
# 2023-09-13
### Load packages ----
library(tidyverse)
library(here)
library(data.table)
library(msm)

rm(list = ls())
### Body ----

# Read in data
data_msm <- read.csv(here("Data", "BTE", "bte_msm_ready.csv"))[,-1]

# Run msm
msm_state <- statetable.msm(sp_class, TESSELLE, data=data_msm)
round(funrar::make_relative(msm_state), 3)

## Define models through Q matrices

# Multinomial analysis to determine impossible transitions
source(here("R-scripts", "03-Analysis", "multinomial.R"))

if (file.exists(here("Data", "msm", "multinom.rds"))) {
  mnom <- readRDS(here("Data", "msm", "multinom.rds"))
} else {
  mnom <- multinom_model(data_msm)
  saveRDS(mnom, here("Data", "msm", "multinom.rds"))
}

Q.model <- as.matrix(round(mnom, 3)) #To determine which transitions are impossible (seems to be none)

## Get initial estimate for Q

Q.init <- crudeinits.msm(sp_class ~ time, TESSELLE, data=data_msm, qmatrix=Q.model)

print("--- *** --- *** ---")

#Run

# msm.sc <- tryCatch(msm( sp_class ~ time, subject=TESSELLE, data = data_msm, 
#                         qmatrix = Q.init, control = list(fnscale = 5000000))
#                    , error = function(e) NA)
# 
# saveRDS(msm.sc, here("R-scripts", "msm_sc.rds"))

print("--- *** --- *** ---")