# Lukas Van Riel
# 2023-09-13
# Main script to execute when running remote. 

### Load libraries ----
library(here) # Load this here to manage paths

### Initialise ----

if(length(args) == 0){
  print("This script needs the zone and then the covariates values as integers!")
  stop()
}

args <- commandArgs(trailingOnly = TRUE)
cat("Arguments passed:", args, "\n")

source(here("R-scripts", "03-Analysis", "Msm.R"))


### Body ----
zone <- commandArgs(trailingOnly = TRUE)[1]
cat("Selected zone:", zone, "\n")
cat("Selected covariates:", cov.subset, "\n")

S <- ifelse(zone=="4bT" || zone=="4cT", 750000, 
                  ifelse(zone=="4eT", 100000,
                         ifelse(zone=="4fS", 200000,
                                ifelse(zone=="4fT", 500000, 300000)))) # Change scaling here
methods <- c("CG")

covariates <- c('cov_Tmean', 'cov_soil', 'cov_CMI', 'cov_pert_class', 
                'cov_pert_sev', 'cov_time_pert')
cov.subset <- covariates[as.numeric(commandArgs(trailingOnly = TRUE)[-1])]


# Create list of all formulas to use
covariate.formula <- character(0)

for (i in 1:length(cov.subset)) {
  covariate_combinations <- combn(cov.subset, i, simplify = FALSE)
  
  for (combo in covariate_combinations) {
    formula_str <- paste("~", paste(combo, collapse = " + "))
    covariate.formula <- c(covariate.formula, formula_str)
  }
}


# Run model for each formula 
for(C in covariate.formula) {
  print(C)
  file.out <- paste0("msm_", zone, ".", methods,".sc", as.character(S), "_", 
         gsub("cov", "", gsub("[~ _]", "", C)), ".rds")
  print(file.out)
  
  if(! file.exists(here("Data-Output", "msm", zone, file.out))) {
    run_remote_msm(data_msm = data_sc, qmatrix = Q.init, ctrl = S, md = methods,
                   cov = C, zone = zone,
                   name.out.rds = file.out)
  }
}
