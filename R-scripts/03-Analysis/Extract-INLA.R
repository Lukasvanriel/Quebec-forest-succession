#### Load packages ####
library(tidyverse)
library(INLA)
library(INLAjoint)
library(here)

### Data ###

inla.all <- readRDS(here("Data-Output", "INLA", "allj.T.C.rds"))

### ### 

names(inla.all)


for (i in seq_along(inla.all)) {
  element_name <- names(inla.all)[i]
  element <- inla.all[[i]]
  
  # Construct the filename
  filename <- here("Data-Output", "INLA", "allj.T.C", paste0("allj.T.C_", element_name, ".RDS"))
  
  # Save the individual element as RDS
  saveRDS(element, filename)
  
  cat("Element", element_name, "saved as", filename, "\n")
}
