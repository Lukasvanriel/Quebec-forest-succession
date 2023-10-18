# Lukas Van Riel
# 2023-09-13
# Main script to execute when running remote. 

### Load libraries ----
library(here) # Load this here to manage paths

### Initialise ----
source(here("R-scripts", "03-Analysis", "Msm.R"))

### Body ----
S = 1 # Change scaling here
methods <- c("CG", "nlm")
covariates <- c('~ cov_Tmean', '~ cov_CMI', '~ cov_soil', 
                '~ cov_CMI + cov_Tmean', '~ cov_Tmean + cov_soil',
                '~ cov_CMI + cov_soil', '~ cov_CMI + cov_Tmean + cov_soil')


covariates <- c('~ cov_pert_class', '~ cov_time_pert', 
                '~ cov_pert_class + cov_Tmean', '~ cov_pert_class + cov_time_pert',
                '~ cov_pert_class + cov_time_pert + cov_Tmean')

covariates <- c('~ cov_Tmean')

for(M in methods) {
  for(C in covariates) {
    run_remote_msm(data_msm = data_sc, qmatrix = Q.init, ctrl = S, md = M,
                   cov = C,
                   name.out.rds = paste0("msm4bM.", M,".sc", as.character(S), "_", 
                                         gsub('[~ cov_]','', C), ".rds"))
  }
}

