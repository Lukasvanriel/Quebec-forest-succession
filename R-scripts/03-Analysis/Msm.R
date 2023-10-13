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
  
  bM <- data_mult_filt %>% 
    filter(SREG_ECO == "4bM")
  
  write.csv(bM, here("Data", "BTE", "bte4bM_msm_ready.csv"))
  
#  plot(data_mult_filt$LONGI, data_mult_filt$LATIT)
#  points(d$LONGI, d$LATIT, col="red")
}

if(F) { # Check up on perturbations
  plot(data$LONGI, data$LATIT)
  points(data4bM$LONGI, data4bM$LATIT, col='red')
  
  data4bM <- read.csv(here("Data", "BTE", "bte4bM_msm_ready.csv"))
  length(table(data4bM$TESSELLE))

  data4bM <- data4bM %>% 
    relocate(TESSELLE, .before=1) %>% 
    relocate(cov_CMI, cov_Tmean, .before=cov_soil) %>% 
    relocate(dom_sp, sp_class, time, .after=GEOCODE)
  
  backup <- data4bM
  
  #Fix Soil
  tab <- unique(data4bM$TESSELLE)
  
  for(t in tab) {
    filtered <- data4bM %>% filter(TESSELLE == t)
    
    soil <- filtered$cov_soil
    print(soil)
    if(sum(! is.na(soil) > 0)) {
      soil.dom <- as.numeric(names(sort(table(soil), decreasing=TRUE)[1]))
      
      #First fix the NA:
      soil[is.na(soil)] <- soil.dom
      #Now the weird values
      soil[soil != soil.dom] <- soil.dom
      print(soil)
    }
    data4bM$cov_soil[data4bM$TESSELLE == t] <- soil

  }  
  
  
  #Fix perturbations
  #data4bM <- backup
  
  tab <- unique(data4bM$TESSELLE)
  
  for(t in tab) {
    filtered <- data4bM %>% filter(TESSELLE == t)
    
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
      data4bM[data4bM$TESSELLE == t,] <- filtered
    }
  } 

  b <- data4bM
  
  # Now update the cov_columns:
  perturbation.matrix <- t(mapply(determine_perturb_class, data4bM$ORIGINE, 
                                  data4bM$AN_ORIGINE, data4bM$year))
  
  data4bM$cov_pert_class <- factor(perturbation.matrix[,1])
  data4bM$cov_pert_sev <- factor(perturbation.matrix[,2])
  data4bM$cov_time_pert <- perturbation.matrix[,3]
  
  #Write out result
  write.csv(data4bM, here("Data", "BTE", "bte4bM_msm_ready.csv"))
}  

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

run_remote_msm <- function(data_msm, qmatrix, md = "BFGS", ctrl = 1, cov = "~ 1", name.out.rds) {
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
  
  tryCatch(saveRDS(msm.model, here("Data-Output", "msm", name.out.rds)), error = function(e) NA)
  
  # Plot model parameters
  tryCatch(plot.msm(msm.model, here("Data-Output", "msm",
                                    paste0(str_sub(name.out.rds, end = -5), "_qvalues.pdf"))),
           error = function(e) NA)
  
}


### Initialise ####
data4bM <- read.csv(here("Data", "BTE", "bte4bM_msm_ready.csv"))

## Create statetables:
msm_state <- statetable.msm(sp_class, TESSELLE, data=data4bM)
round(funrar::make_relative(msm_state), 3)

## Define models through Q matrices

# Multinomial analysis to determine impossible transitions

source(here("R-scripts", "03-Analysis", "Multinomial.R"))

if (file.exists(here("Data-Output", "msm", "multinom4bM.rds"))) {
  mnom <- readRDS(here("Data-Output", "msm", "multinom4bM.rds"))
} else {
  mnom <- multinom_model(data4bM)
  saveRDS(mnom, here("Data-Output", "msm", "multinom4bM.rds"))
}

#To determine which transitions are impossible
Q.model <- as.matrix(round(mnom, 3)) 

## Get initial estimate for Q

Q.init <- crudeinits.msm(sp_class ~ time, TESSELLE, data=data4bM, qmatrix=Q.model)


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

data_sc <- data4bM
data_sc$cov_CMI[is.na(data_sc$cov_CMI)] <- mean(data_sc$cov_CMI, na.rm =T)
data_sc$cov_Tmean[is.na(data_sc$cov_Tmean)] <- mean(data_sc$cov_Tmean, na.rm =T)
data_sc <- data_sc %>% 
  mutate(cov_CMI = (cov_CMI - mean(cov_CMI)) / sd(cov_CMI)) %>% 
  mutate(cov_Tmean = (cov_Tmean - mean(cov_Tmean)) / sd(cov_Tmean))


#names.subs10 <- names(table(data_test$TESSELLE))[sample.int(length(table(data_test$TESSELLE)), 
#                                                            0.1 * length(table(data_test$TESSELLE)), replace = FALSE)]

#data_subs10 <- subset(data_test, TESSELLE %in% names.subs10)


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


if(F) {
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
}
