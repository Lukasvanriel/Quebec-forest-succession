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

zone <- commandArgs(trailingOnly = TRUE)[1]
### Functions ####

#Function to prepare dataset for msm
prepare_data <- function(dataset){
  ## Add time column that indicates time since first observation
  #First remove rows that lack observation times
  
  data_time <- dataset |>
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

}

time_since <- function(pert, meas) {
  if(is.na(pert)) {
    return(-1)
  } else {return(meas - pert)}
}

determine_perturb_class <- function(pert_string, pert_time, meas_time) {
  if (is.na(pert_string)) {
    pert_class <- 0
    pert_sev <- 0
    time_since_pert <- 200
  } else if(pert_string %in% c("BR", "BRD", "BRU")) {
    pert_class <- 1
    pert_sev <- 2 # A "1" can be added later if we want to look dat partial perturbations
    time_since_pert <- time_since(pert_time, meas_time)
  } else if (pert_string %in% c("CBA","CBT","CEF","CPT","CRB","CRS","CS","CT","ETR","RPS")) {
    pert_class <- 2
    pert_sev <- 2 # A "1" can be added later if we want to look dat partial perturbations
    time_since_pert <- time_since(pert_time, meas_time)
  } else if (pert_string %in% c("ES")) {
    pert_class <- 3
    pert_sev <- 2 # A "1" can be added later if we want to look dat partial perturbations
    time_since_pert <- time_since(pert_time, meas_time)
  } else {
    pert_class <- 0
    pert_sev <- 0
    time_since_pert <- 200
  }
  c(pert_class, pert_sev, time_since_pert)
}

#Function to subset the bte according to several zones the and fix the covariates 
subset_ecotone <- function(data.bte, zone) {
  ## Subset the ecotone to contain only the specified zone
  data.zone <- data.bte %>% 
    filter(SREG_ECO == zone) %>% 
    relocate(TESSELLE, .before=1) %>% 
    relocate(cov_CMI, cov_Tmean, .before=cov_soil) %>% 
    relocate(dom_sp, sp_class, time, .after=GEOCODE)
  
  backup <- data.zone
  ## Now fix covariates

  # Fix Soil
  tab <- unique(data.zone$TESSELLE)
  
  for(t in tab) {
    filtered <- data.zone %>% filter(TESSELLE == t)
    
    soil <- filtered$cov_soil
    #print(soil)
    if(sum(! is.na(soil) > 0)) {
      soil.dom <- as.numeric(names(sort(table(soil), decreasing=TRUE)[1]))
      
      #First fix the NA:
      soil[is.na(soil)] <- soil.dom
      #Now the weird values
      soil[soil != soil.dom] <- soil.dom
      print(soil)
    }
    data.zone$cov_soil[data.zone$TESSELLE == t] <- soil
  }  
  
  # Fix perturbations
  
  for(t in tab) {
    filtered <- data.zone %>% filter(TESSELLE == t)
    
    # First only proceed with tesselle that have at least one perturbation
    if(max(filtered$cov_pert_class) > 0 && sum(!is.na(filtered$AN_ORIGINE)) > 0){
      pert <- cbind(filtered$ORIGINE, filtered$AN_ORIGINE)
      
      pert_filt <- unique(matrix(pert[complete.cases(pert), ], ncol=2))
      
      if(nrow(pert_filt) == 1) {
        pert_filt <- rbind(pert_filt, c("FAKE", 3000)) #Add artificial line that will never be used
      } else(print(t))
      
      pert_filt <- pert_filt[order(pert_filt[,2],decreasing=F),]
      
      p.type <- pert_filt[1, 1]
      p.year <- as.numeric(pert_filt[1, 2])
      
      for(k in 1:nrow(filtered)){
        c.year <- filtered$year[k]
        
        if(c.year > as.numeric(pert_filt[2, 2])) {
          p.type <- pert_filt[2, 1]
          p.year <- as.numeric(pert_filt[2, 2])
        }
        
        if(c.year > p.year){
          filtered$ORIGINE[k] <- p.type
          filtered$AN_ORIGINE[k] <- p.year
        }
      }
      data.zone[data.zone$TESSELLE == t,] <- filtered
    }
  } 
  
  # Now update the cov_columns:
  perturbation.matrix <- t(mapply(determine_perturb_class, data.zone$ORIGINE, 
                                  data.zone$AN_ORIGINE, data.zone$year))
  
  data.zone$cov_pert_class <- factor(perturbation.matrix[,1])
  data.zone$cov_pert_sev <- factor(perturbation.matrix[,2])
  data.zone$cov_time_pert <- perturbation.matrix[,3]
  
  data.zone
  
}

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

#Function to easily run msm from the terminal; Catches errors and makes sure the script doesn't break down
run_remote_msm <- function(data_msm, qmatrix, md = "BFGS", ctrl = 1, cov = "~ 1", zone=NA, name.out.rds) {
  if(! all(c("sp_class", "time", "TESSELLE") %in% colnames(data_msm))) {
    return("Missing information in data.")
  }
  
  if(md %in% c("BFGS", "CG", "Nelder-Mead", "SANN")) {
    msm.model <- tryCatch(msm( sp_class ~ time, subject=TESSELLE, data = data_msm,
                               qmatrix = qmatrix, method= md, control = list(fnscale = ctrl),
                               covariates = as.formula(cov)),
                          error = function(e) NA)
    
  } else if (md %in% c("nlm", "bobyqa", "fisher")){
    msm.model <- tryCatch(msm( sp_class ~ time, subject=TESSELLE, data = data_msm,
                               qmatrix = qmatrix, opt.method= md, control = list(fnscale = ctrl),
                               covariates = as.formula(cov)),
                          error = function(e) NA)
  }
  
  tryCatch(saveRDS(msm.model, here("Data-Output", "msm", zone, name.out.rds)), error = function(e) NA)
  
  # Plot model parameters
  tryCatch(plot.msm(msm.model, here("Data-Output", "msm", zone,
                                    paste0(str_sub(name.out.rds, end = -5), "_qvalues.pdf"))),
           error = function(e) NA)
  
}

### Data ####

#data <- read.csv(here("Data", "BTE", "bte_cov_class.csv"))[,-1]
#write.csv(prepare_data(data), here("Data", "BTE", "bte_msm_ready.csv"))

if(F) {
data.bte <- read.csv(here("Data", "BTE", "bte_msm_ready.csv")) %>% 
  select(-X)

for(z in names(table(data.bte$SREG_ECO)[table(data.bte$SREG_ECO)>10000]) ) {
  data.subset <- subset_ecotone(data.bte, z)
  write.csv(data.subset, here("Data", "BTE", paste0("bte_",z, "_msm_ready.csv")))
}
}

##
args <- commandArgs(trailingOnly = TRUE)

##
data <- read_csv(here("Data", "BTE", paste0("bte_", args[1], "_msm_ready.csv")),
                    col_types = cols(.default = col_guess(),
                                     sp_class = col_integer(),
                                     cov_pert_class = col_factor(),
                                     cov_pert_sev = col_factor(),
                                     cov_time_pert = col_double()))[,-1]# %>% select(-"...X")

# Rescale covariates where necessary. Add column classes
data_sc <- data
data_sc$cov_CMI[is.na(data_sc$cov_CMI)] <- mean(data_sc$cov_CMI, na.rm =T)
data_sc$cov_Tmean[is.na(data_sc$cov_Tmean)] <- mean(data_sc$cov_Tmean, na.rm =T)
data_sc <- data_sc %>% 
  mutate(cov_CMI = (cov_CMI - mean(cov_CMI)) / sd(cov_CMI)) %>% 
  mutate(cov_Tmean = (cov_Tmean - mean(cov_Tmean)) / sd(cov_Tmean)) %>% 
  mutate(cov_time_pert = cov_time_pert / 100) %>% 
  mutate(cov_soil = cov_soil / 10)

### Initialise ####

## Create statetables:
msm_state <- statetable.msm(sp_class, TESSELLE, data=data_sc)
round(funrar::make_relative(msm_state), 3)

## Define models through Q matrices

# Multinomial analysis to determine impossible transitions

source(here("R-scripts", "03-Analysis", "Multinomial.R"))

if (file.exists(here("Data-Output", "msm", zone, paste0("multinom_", args[1], ".rds")))) {
  mnom <- readRDS(here("Data-Output", "msm", zone, paste0("multinom_", args[1], ".rds")))
} else {
  mnom <- multinom_model(data_sc)
  saveRDS(mnom, here("Data-Output", "msm", zone, paste0("multinom_", args[1], ".rds")))
}

#To determine which transitions are impossible
Q.model <- as.matrix(round(mnom, 3)) 

## Get initial estimate for Q

Q.init <- crudeinits.msm(sp_class ~ time, TESSELLE, data=data_sc, qmatrix=Q.model)


### Run msm ####

# run_remote_msm(data_msm = data_msm, qmatrix = Q.init, ctrl = 1, name.out.rds = "msm.reg.rds")
# run_remote_msm(data_msm = data_msm, qmatrix = Q.init, ctrl = 5000000, name.out.rds = "msm.reg.sc5M.rds")


## Lets create nested loops to try out different things:

if(F) {
  #scaling <- c(5000000, 1)
  #methods <- c("BFGS", "CG", "Nelder-Mead", "nlm", "fisher") #Removed SANN and bobyqa
  
  scaling <- c(1)
  methods <- c("Nelder-Mead", "nlm", "fisher")  
  
  for(S in scaling) {
    for(M in methods) {
      run_remote_msm(data_msm = data_msm, qmatrix = Q.init, ctrl = S, md = M,
                     name.out.rds = paste0("msm.", M,".sc", as.character(S), ".rds"))
    }
  }
}

### Try nested loop for covariates

if(F) {
  S <- 500000
  covariates <- c('~ cov_CMI', '~ cov_Tmean', '~ cov_soil', 
                  '~ cov_CMI + cov_Tmean', '~ cov_Tmean + cov_soil',
                  '~ cov_CMI + cov_soil', '~ cov_CMI + cov_Tmean + cov_soil')
  methods <- c("CG", "nlm") # removed bobyqa and SANN
  
  for(M in methods) {
    for(C in covariates) {
      run_remote_msm(data_msm = data_sc, qmatrix = Q.init, ctrl = S, md = M,
                     cov = C,
                     name.out.rds = paste0("msm4bM.", M,".sc", as.character(S), "_", 
                                           gsub('[~ cov_]','', C), ".rds"))
    }
  }
  
}

## Now for perturbation

if(F) {
  S <- 500000
  covariates <- c('~ cov_pert_class', '~ cov_time_pert', 
                  '~ cov_pert_class + cov_Tmean', '~ cov_pert_class + cov_time_pert',
                  '~ cov_pert_class + cov_time_pert + cov_Tmean')
  methods <- c("CG")
  
  for(M in methods) {
    for(C in covariates) {
      run_remote_msm(data_msm = data_sc, qmatrix = Q.init, ctrl = S, md = M,
                     cov = C,
                     name.out.rds = paste0("msm4bM.", M,".sc", as.character(S), "_", 
                                           gsub('[~ cov_]','', C), ".rds"))
    }
  }
  
}





### For future runs ####
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

}
