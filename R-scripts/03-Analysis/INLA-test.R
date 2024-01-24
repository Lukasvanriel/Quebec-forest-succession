#### Load packages ####
library(tidyverse)
library(INLA)
library(INLAjoint)
library(here)

### Data ###
#Changes to dataset can be made in Msm.R; Or need to structure it differently
data_msm <- read.csv(here("Data", "BTE", "bte_msm_ready.csv"))[,-1]

#Filter out smaller part of BTE.
#plot(data_msm$LONGI, data_msm$LATIT)
#points(filter(data_msm, SREG_ECO == "4bM")$LONGI, filter(data_msm, SREG_ECO == "4bM")$LATIT, col="red")
#points(filter(data_msm, SREG_ECO == "4bS")$LONGI, filter(data_msm, SREG_ECO == "4bS")$LATIT, col="green")
#points(filter(data_msm, SREG_ECO == "4bT")$LONGI, filter(data_msm, SREG_ECO == "4bT")$LATIT, col="blue")

data_msm_4b <- read.csv(here("Data", "BTE", "bte_msm_ready.csv"))[,-1] %>% 
  filter(SREG_ECO %in% c("4bM", "4bS", "4bT"))

plot(data_msm_4b$LONGI, data_msm_4b$LATIT)
nrow(data_msm_4b) / nrow(data_msm)

### Prepare dataset ###

# Select relevant columns
data.inla.inter <- data_msm_4b |>
  select(TESSELLE, time, sp_class, LONGI, LATIT, 
         cov_CMI, cov_Tmean, cov_soil, cov_pert_class, cov_pert_sev, cov_time_pert) |>
  rename(ID=TESSELLE) |>
  mutate(from=sp_class) |>
  mutate(to=sp_class) |>
  mutate(entry=time) |>
  relocate(from, to, entry, time, .after = ID) 

# Create

ph.from <- c(0, data.inla.inter$to)
data.inla.inter$from <- ph.from[1:(length(ph.from)-1)]

ph.entry <- c(0, data.inla.inter$time)
data.inla.inter$entry <- ph.entry[1:(length(ph.entry)-1)]

if(F) {
  data.inla <- data.inla.inter %>% 
    filter(time>0) %>% 
    mutate(trans=to-from)  %>% 
    group_by(ID) %>% 
    #  mutate(keep = ifelse(trans==0 & row_number() != n(), 0, 1)) %>% 
    #  mutate(entry_upd = ifelse(trans==0 & row_number() == n(), entry[row_number() - 1], entry)) %>% 
    relocate(entry_upd,trans, keep, .before = sp_class)
  
  data.inla[1:15,1:9]
  #if trans=0 and not last line: add entry time to next + delete line or something
} # Do something here if non-trantions need to be cut 

data.inla <- data.inla.inter |>
  filter(time>0) |>
  select(-sp_class) # |> ,-trans
# filter(time != 0) #this used to be trans!=0 because I used to filter out the non-transitions

data.inla.inter[1:14,]
data.inla[1:10,]

### Run INLA ###

##State table:
st <- matrix(0, ncol=9, nrow=9)

for(i in 1:nrow(data.inla)) {
  st[data.inla$from[i], data.inla$to[i]] <- st[data.inla$from[i], data.inla$to[i]] + 1
}
st
# Maybe remove some of the Jack Pine transitions?

##Prepare the survival objects
#First create the empty event dataframes:
helper <- seq(1,9, by=1)
transitions <- c()
for(i in helper) {
  for(j in helper[-i]) {
    transitions <- c(transitions, 10*i + j%%10)
  }
}

event.list <- vector("list", transitions[length(transitions)])
#if(! file.exists(here("Data", "BTE", "INLA", "Event_"), as.character(tr), ".csv", sep="")))
if(F) {
  for(i in transitions) {
    event.list[[i]] <- data.frame(matrix(ncol = ncol(data.inla) + 1))
    colnames(event.list[[i]]) <- c(colnames(data.inla), "status")
    event.list[[i]] <- event.list[[i]][-1,]
  }
  
  # Now fill in the full event dataframes for each transition
  for(tr in transitions) { 
    print(tr)
    data.tr <- data.inla |>
      filter(from == tr %/% 10)
    
    tr.event <- vector("list", length = nrow(data.tr))
    for(i in 1:nrow(data.tr)) { 
      if(data.tr$to[i] == tr %% 10) {
        tr.event[[i]] <- c(data.tr[i,], 1)
        names(tr.event[[i]])[length(tr.event[[i]])] <- "status"
      } else {
        tr.event[[i]] <- c(data.tr[i,], 0)
        tr.event[[i]]$to <- tr %% 10
        names(tr.event[[i]])[length(tr.event[[i]])] <- "status"
      }
    }
    event.list[[tr]] <- bind_rows(tr.event)
    write.csv(event.list[[tr]], 
              paste(here("Data", "BTE", "INLA", "4b", "Event_"), as.character(tr), ".csv", sep=""))
  }
  
}

##Read in event datasets
for(tr in transitions) {
  event.list[[tr]] <- read.csv(paste(here("Data", "BTE", "INLA", "4b", "Event_"), as.character(tr), ".csv", sep="")) %>% 
    subset(select = -X)
}

## Rescale necessary climatic covariates  + Fix missing soil values (for now just take as 0)
fix_covariates <- function(df) {
  if (!is.null(df)) {
    df$cov_Tmean[is.na(df$cov_Tmean)] <- mean(df$cov_Tmean, na.rm = T)
    df$cov_CMI[is.na(df$cov_CMI)] <- mean(df$cov_CMI, na.rm = T)
    
    df$cov_Tmean_scaled <-  (df$cov_Tmean - mean(df$cov_Tmean)) / sd(df$cov_Tmean) 
    df$cov_CMI_scaled <-  (df$cov_CMI - mean(df$cov_CMI)) / sd(df$cov_CMI)
    df$cov_soil_compl <-  ifelse(is.na(df$cov_soil), 0, df$cov_soil)
  }
  return(df)
}

event.list.ext <- lapply(event.list, fix_covariates)


### Run INLA ###

inla.31.T <- joint(formSurv=list(inla.surv(event.list.ext[[31]]$time, event.list.ext[[31]]$status) ~ cov_Tmean),
                   basRisk = rep("weibullsurv", 1),
                   dataSurv = list(event.list.ext[[31]]))

summary(inla.31.T)
plot(inla.31.T)$Baseline
plot(inla.31.T)$Outcomes

saveRDS(inla.31.T, here("Data-Output", "INLA", "31.T.rds"))

a <- readRDS(here("Data-Output", "INLA", "31.T.rds"))

summary(a)
plot(a)$Baseline
plot(a)$Outcomes

