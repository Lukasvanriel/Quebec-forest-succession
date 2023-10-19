# Lukas Van Riel
# 2023-09-13
# Main script to execute when running remote. 

### Load libraries ----
library(here) # Load this here to manage paths

### Initialise ----

source(here("R-scripts", "03-Analysis", "Msm.R"))

args <- as.numeric(commandArgs(trailingOnly = TRUE))
cat("Arguments passed:", args, "\n")


### Body ----
S = 300000 # Change scaling here
methods <- c("CG")

covariates <- c('cov_Tmean', 'cov_soil', 'cov_CMI', 'cov_pert_class', 
                'cov_pert_sev', 'cov_time_pert')
cov.subset <- covariates[args]

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
  run_remote_msm(data_msm = data_sc, qmatrix = Q.init, ctrl = S, md = methods,
                 cov = C,
                 name.out.rds = paste0("msm4bM.", M,".sc", as.character(S), "_", 
                                       gsub("cov", "", gsub("[~ _]", "", C)), ".rds"))
}
