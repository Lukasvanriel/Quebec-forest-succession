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
#nrow(data_msm_4b) / nrow(data_msm)

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

if(F) { # Individual 3-1 transition, cov=T
inla.31.T <- joint(formSurv=list(inla.surv(event.list.ext[[31]]$time, event.list.ext[[31]]$status) ~ cov_Tmean),
                   basRisk = rep("weibullsurv", 1),
                   dataSurv = list(event.list.ext[[31]]))

summary(inla.31.T)
#plot(inla.31.T)$Baseline
#plot(inla.31.T)$Outcomes

saveRDS(inla.31.T, here("Data-Output", "INLA", "31.T.rds"))
}

if(F) {# Each individual transition starting from 3, cov=T
  inla.3.T <- vector("list", 8)
  
  inla.3.T <- lapply(c(31, 32, 34:39), function(x) {
    print(x)
    joint(formSurv= list(inla.surv(time, status) ~ cov_Tmean),
          basRisk = rep("weibullsurv", 1),
          dataSurv = list(event.list.ext[[x]]))})
  
  saveRDS(inla.3.T, here("Data-Output", "INLA", "3.T.rds"))
  
  for(i in 1:8) {print(summary(inla.3.T[[i]]))}
}

if(F) {# Each individual transition starting from 3, cov=T,C,S
  inla.3.T.C.S <- vector("list", 8)
  
  inla.3.T.C.S <- lapply(c(31, 32, 34:39), function(x) {
    print(x)
    joint(formSurv= list(inla.surv(time, status) ~ cov_Tmean + cov_CMI + cov_soil),
          basRisk = rep("weibullsurv", 1),
          dataSurv = list(event.list.ext[[x]]))})
  
  saveRDS(inla.3.T.C.S, here("Data-Output", "INLA", "3.T.C.S.rds"))
  
  for(i in 1:8) {summary(inla.3.T.C.S[[i]])}
}

if(F) {# Joint transitions starting from 3, cov=T
  inla3.j.T <- joint(formSurv=list(inla.surv(event.list.ext[[31]]$time, event.list.ext[[31]]$status) ~ cov_Tmean,
                                   inla.surv(event.list.ext[[32]]$time, event.list.ext[[32]]$status) ~ cov_Tmean,
                                   inla.surv(event.list.ext[[34]]$time, event.list.ext[[34]]$status) ~ cov_Tmean,
                                   inla.surv(event.list.ext[[35]]$time, event.list.ext[[35]]$status) ~ cov_Tmean,
                                   inla.surv(event.list.ext[[36]]$time, event.list.ext[[36]]$status) ~ cov_Tmean,
                                   inla.surv(event.list.ext[[37]]$time, event.list.ext[[37]]$status) ~ cov_Tmean,
                                   inla.surv(event.list.ext[[38]]$time, event.list.ext[[38]]$status) ~ cov_Tmean,
                                   inla.surv(event.list.ext[[39]]$time, event.list.ext[[39]]$status) ~ cov_Tmean),
                     basRisk = rep("weibullsurv", 8),
                     dataSurv = list(event.list.ext[[31]],
                                     event.list.ext[[32]],
                                     event.list.ext[[34]],
                                     event.list.ext[[35]],
                                     event.list.ext[[36]],
                                     event.list.ext[[37]],
                                     event.list.ext[[38]],
                                     event.list.ext[[39]]))
  
  saveRDS(inla3.j.T, here("Data-Output", "INLA", "3j.T.rds"))
  
  
  summary(inla3.j.T)
}

if(F) {# Joint transitions starting from 3, cov=T,C,S
  inla3.j.T.C.S <- joint(formSurv=list(inla.surv(event.list.ext[[31]]$time, event.list.ext[[31]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                   inla.surv(event.list.ext[[32]]$time, event.list.ext[[32]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                   inla.surv(event.list.ext[[34]]$time, event.list.ext[[34]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                   inla.surv(event.list.ext[[35]]$time, event.list.ext[[35]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                   inla.surv(event.list.ext[[36]]$time, event.list.ext[[36]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                   inla.surv(event.list.ext[[37]]$time, event.list.ext[[37]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                   inla.surv(event.list.ext[[38]]$time, event.list.ext[[38]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                   inla.surv(event.list.ext[[39]]$time, event.list.ext[[39]]$status) ~ cov_Tmean + cov_CMI + cov_soil),
                     basRisk = rep("weibullsurv", 8),
                     dataSurv = list(event.list.ext[[31]],
                                     event.list.ext[[32]],
                                     event.list.ext[[34]],
                                     event.list.ext[[35]],
                                     event.list.ext[[36]],
                                     event.list.ext[[37]],
                                     event.list.ext[[38]],
                                     event.list.ext[[39]]))
  
  saveRDS(inla3.j.T.C.S, here("Data-Output", "INLA", "3j.T.C.S.rds"))
  
  summary(inla3.j.T.C.S)
}

if(T) {# Joint transitions for all, cov=T,C,S
  inlaA.j.T.C.S <- joint(formSurv=list(inla.surv(event.list.ext[[12]]$time, event.list.ext[[12]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[13]]$time, event.list.ext[[13]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[14]]$time, event.list.ext[[14]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[15]]$time, event.list.ext[[15]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[16]]$time, event.list.ext[[16]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[17]]$time, event.list.ext[[17]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[18]]$time, event.list.ext[[18]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[19]]$time, event.list.ext[[19]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[21]]$time, event.list.ext[[21]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[23]]$time, event.list.ext[[23]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[24]]$time, event.list.ext[[24]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[25]]$time, event.list.ext[[25]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[26]]$time, event.list.ext[[26]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[27]]$time, event.list.ext[[27]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[28]]$time, event.list.ext[[28]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[29]]$time, event.list.ext[[29]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[31]]$time, event.list.ext[[31]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[32]]$time, event.list.ext[[32]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[34]]$time, event.list.ext[[34]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[35]]$time, event.list.ext[[35]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[36]]$time, event.list.ext[[36]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[37]]$time, event.list.ext[[37]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[38]]$time, event.list.ext[[38]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[39]]$time, event.list.ext[[39]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[41]]$time, event.list.ext[[41]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[42]]$time, event.list.ext[[42]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[43]]$time, event.list.ext[[43]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[45]]$time, event.list.ext[[45]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[46]]$time, event.list.ext[[46]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[47]]$time, event.list.ext[[47]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[48]]$time, event.list.ext[[48]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[49]]$time, event.list.ext[[49]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[51]]$time, event.list.ext[[51]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[52]]$time, event.list.ext[[52]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[53]]$time, event.list.ext[[53]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[54]]$time, event.list.ext[[54]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[56]]$time, event.list.ext[[56]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[57]]$time, event.list.ext[[57]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[58]]$time, event.list.ext[[58]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[59]]$time, event.list.ext[[59]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[61]]$time, event.list.ext[[61]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[62]]$time, event.list.ext[[62]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[63]]$time, event.list.ext[[63]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[64]]$time, event.list.ext[[64]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[65]]$time, event.list.ext[[65]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[67]]$time, event.list.ext[[67]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[68]]$time, event.list.ext[[68]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[69]]$time, event.list.ext[[69]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[71]]$time, event.list.ext[[71]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[72]]$time, event.list.ext[[72]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[73]]$time, event.list.ext[[73]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[74]]$time, event.list.ext[[74]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[75]]$time, event.list.ext[[75]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[76]]$time, event.list.ext[[76]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[78]]$time, event.list.ext[[78]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[79]]$time, event.list.ext[[79]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[81]]$time, event.list.ext[[81]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[82]]$time, event.list.ext[[82]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[83]]$time, event.list.ext[[83]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[84]]$time, event.list.ext[[84]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[85]]$time, event.list.ext[[85]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[86]]$time, event.list.ext[[86]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[87]]$time, event.list.ext[[87]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[89]]$time, event.list.ext[[89]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[91]]$time, event.list.ext[[91]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[92]]$time, event.list.ext[[92]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[93]]$time, event.list.ext[[93]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[94]]$time, event.list.ext[[94]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[95]]$time, event.list.ext[[95]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[96]]$time, event.list.ext[[96]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[97]]$time, event.list.ext[[97]]$status) ~ cov_Tmean + cov_CMI + cov_soil,
                                       inla.surv(event.list.ext[[98]]$time, event.list.ext[[98]]$status) ~ cov_Tmean + cov_CMI + cov_soil),
                         basRisk = rep("weibullsurv", 72),
                         dataSurv = list(event.list.ext[[12]], event.list.ext[[13]], event.list.ext[[14]], event.list.ext[[15]], event.list.ext[[16]], event.list.ext[[17]], event.list.ext[[18]], event.list.ext[[19]],
                                         event.list.ext[[21]], event.list.ext[[23]], event.list.ext[[24]], event.list.ext[[25]], event.list.ext[[26]], event.list.ext[[27]], event.list.ext[[28]], event.list.ext[[29]],
                                         event.list.ext[[31]], event.list.ext[[32]], event.list.ext[[34]], event.list.ext[[35]], event.list.ext[[36]], event.list.ext[[37]], event.list.ext[[38]], event.list.ext[[39]],
                                         event.list.ext[[41]], event.list.ext[[42]], event.list.ext[[43]], event.list.ext[[45]], event.list.ext[[46]], event.list.ext[[47]], event.list.ext[[48]], event.list.ext[[49]],
                                         event.list.ext[[51]], event.list.ext[[52]], event.list.ext[[53]], event.list.ext[[54]], event.list.ext[[56]], event.list.ext[[57]], event.list.ext[[58]], event.list.ext[[59]],
                                         event.list.ext[[61]], event.list.ext[[62]], event.list.ext[[63]], event.list.ext[[64]], event.list.ext[[65]], event.list.ext[[67]], event.list.ext[[68]], event.list.ext[[69]],
                                         event.list.ext[[71]], event.list.ext[[72]], event.list.ext[[73]], event.list.ext[[74]], event.list.ext[[75]], event.list.ext[[76]], event.list.ext[[78]], event.list.ext[[79]],
                                         event.list.ext[[81]], event.list.ext[[82]], event.list.ext[[83]], event.list.ext[[84]], event.list.ext[[85]], event.list.ext[[86]], event.list.ext[[87]], event.list.ext[[89]],
                                         event.list.ext[[91]], event.list.ext[[92]], event.list.ext[[93]], event.list.ext[[94]], event.list.ext[[95]], event.list.ext[[96]], event.list.ext[[97]], event.list.ext[[98]]))
  
  saveRDS(inlaA.j.T.C.S, here("Data-Output", "INLA", "allj.T.C.S.rds"))
  
  summary(inlaA.j.T.C.S)
}
