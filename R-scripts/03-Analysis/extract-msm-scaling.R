# Lukas Van Riel
# 2023-10-31
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

### extract models and check their scalings ####

# Load all folders
# Second line because there was an old folder that did not need to be analysis
zone.list <- list.dirs(path = here("Data-Output", "msm"), full.names = F, recursive = F)
zone.list <- zone.list[zone.list!="OLD"]
  
# Check what zones are present
for (z in zone.list) {
  # Get file list, but drop the multinom.*.rds file
  model.list <- list.files(path = here("Data-Output", "msm", z), full.names = F, recursive = F)
  model.list <- model.list[! grepl("^multinom", model.list)]
  
  print(model.list)
  print(length(model.list))
  if(length(model.list) > 0){
    for(m in model.list) {
      model <- readRDS(here("Data-Output", "msm", z, m))
      if(! is.list(model)) {
        print(here("Data-Output", "msm", z, m))
        print("NO MODEL OUTPUT")
      } else{
        print(here("Data-Output", "msm", z, m))
        print(model$minus2loglik)
        }
    }
  }
}

model.list <- list.files(path = here("Data-Output", "msm", "4bT"), full.names = F, recursive = F)
a <- readRDS(here("Data-Output", "msm", "4bS", "msm_4bS.CG.sc3e+05_Tmean.rds"))
