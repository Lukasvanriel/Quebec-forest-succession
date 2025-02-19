#####Check with Lat/Long

#### Runs the markov chain model using the INLA method ####

#### Load packages ####
library(tidyverse)
library(INLA)
library(INLAjoint)
library(here)

### Data ###
#Changes to dataset can be made in Msm.R; Or need to structure it differently
data_msm2 <- read.csv(here("Data", "BTE", "bte_msm_ready.csv"))[,-1]

### Prepare dataset ###

data.inla.inter <- data_msm2 |>
  select(TESSELLE, time, sp_class, LONGI, LATIT) |>
  mutate(ID=TESSELLE) |>
  mutate(from=sp_class) |>
  mutate(to=sp_class) |>
  mutate(entry=time) |>
  relocate(time, .after = last_col())

ph.from <- c(0, data.inla.inter$to)
data.inla.inter$from <- ph.from[1:(length(ph.from)-1)]

ph.entry <- c(0, data.inla.inter$time)
data.inla.inter$entry <- ph.entry[1:(length(ph.entry)-1)] 

data.inla <- data.inla.inter |>
  filter(time>0) |>
  select(-sp_class) |>
  mutate(trans=to-from) |>
  filter(trans != 0) |> 
  ungroup() |>
  select(-TESSELLE, -trans)

### Run INLA ###

##State table:
st <- matrix(0, ncol=9, nrow=9)

for(i in 1:nrow(data.inla)) {
  st[data.inla$from[i], data.inla$to[i]] <- st[data.inla$from[i], data.inla$to[i]] + 1
}
st

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
              paste(here("Data", "BTE", "INLA", "Event_OLD_"), as.character(tr), ".csv", sep=""))
  }
  
}

##Read in event datasets
for(tr in transitions) {
  event.list[[tr]] <- read.csv(paste(here("Data", "BTE", "INLA", "Event_OLD_"), as.character(tr), ".csv", sep=""),
                               colClasses = c("NULL",NA,NA,NA,NA,NA,NA,NA,NA))
}

## Rescale necessary columns
#TODO: Compare running times of rescaled with unscaled.

#Create survival objects
Surv.list <- vector("list", transitions[length(transitions)])
for(tr in transitions) { 
  Surv.list[[tr]] <- inla.surv(time = event.list[[tr]]$time, event = event.list[[tr]]$status)
}

### Run INLA ###

inla.test1 <- joint(formSurv=list(inla.surv(event.list[[12]]$time, event.list[[12]]$status) ~ LATIT + LONGI),
                   basRisk = rep("weibullsurv", 1),
                   dataSurv = list(event.list[[12]]))

inla1 <- lapply(1:8, function(x) {
  print(x)
  joint(formSurv= list(inla.surv(time, status) ~ LATIT + LONGI),
                  basRisk = rep("weibullsurv", 1),
                  dataSurv = list(event.list[[x+11]]))})
backup <- inla1

for(i in 1:8) {print(summary(inla1[[i]]))}


summary(inla1[[1]])
summary(inla1[[2]])
summary(inla1[[3]])

plot(inla1[[1]])$Baseline
plot(inla1[[1]])$Outcomes
plot(inla1[[2]])$Baseline
plot(inla1[[2]])$Outcomes
plot(inla1[[3]])$Baseline
plot(inla1[[3]])$Outcomes


inla.test2 <- joint(formSurv=list(inla.surv(event.list[[12]]$time, event.list[[12]]$status) ~ LATIT + LONGI,
                                 inla.surv(event.list[[13]]$time, event.list[[13]]$status) ~ LATIT + LONGI),
                 basRisk = rep("weibullsurv", 2),
                 dataSurv = list(event.list[[12]],
                                 event.list[[13]]))

summary(inla.test2)
summary(inla1[[1]])
summary(inla1[[2]])


###OLD


a12 <- Surv.list[[12]]
a13 <- Surv.list[[13]]
a14 <- Surv.list[[14]]
a15 <- Surv.list[[15]]
a16 <- Surv.list[[16]]
a17 <- Surv.list[[17]]
a18 <- Surv.list[[18]]
a19 <- Surv.list[[19]]


inla.12 <- joint(formSurv=list(a12 ~ LATIT + LONGI),
                basRisk = rep("weibullsurv", 1),
                dataSurv = list(event.list[[12]]))

event.list[[12]]

summary(inla.12)

plot(inla.12)$Baseline
plot(inla.12)$Outcomes


inla.12.13 <- joint(formSurv=list(a12 ~ LATIT + LONGI,
                                  a13 ~ LATIT + LONGI),
                    basRisk = rep("weibullsurv", 2),
                    dataSurv = list(event.list[[12]],
                                    event.list[[13]]))

summary(inla.12.13)

plot(inla.12.13)$Baseline
plot(inla.12.13)$Outcomes

inla.1 <- joint(formSurv=list(a12 ~ LATIT + LONGI,
                              a13 ~ LATIT + LONGI,
                              a14 ~ LATIT + LONGI,
                              a15 ~ LATIT + LONGI,
                              a16 ~ LATIT + LONGI,
                              a17 ~ LATIT + LONGI,
                              a18 ~ LATIT + LONGI,
                              a19 ~ LATIT + LONGI),
                basRisk = rep("weibullsurv", 8),
                dataSurv = list(event.list[[12]],
                                event.list[[13]],
                                event.list[[14]],
                                event.list[[15]],
                                event.list[[16]],
                                event.list[[17]],
                                event.list[[18]],
                                event.list[[19]]))

for(i in 2:9) {
  print(i)
  t <- 1*10 + i
  
  inla.12 <- joint(formSurv=list(get(paste0("a", as.character(t))) ~ LATIT + LONGI),
                   basRisk = rep("weibullsurv", 1),
                   dataSurv = list(event.list[[t]]))
}

summary(inla.1)

plot(inla.1)$Baseline
plot(inla.1)$Outcomes

###Try with actual covariates:
#####Check with Lat/Long

#### Runs the markov chain model using the INLA method ####

#### Load packages ####
library(tidyverse)
library(INLA)
library(INLAjoint)
library(here)

### Data ###
#Changes to dataset can be made in Msm.R; Or need to structure it differently
data_msm <- read.csv(here("Data", "BTE", "bte_msm_ready.csv"))[,-1]

### Prepare dataset ###

# Select relevant columns
data.inla.inter <- data_msm |>
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
              paste(here("Data", "BTE", "INLA", "Event_"), as.character(tr), ".csv", sep=""))
  }
  
}

##Read in event datasets
for(tr in transitions) {
  event.list[[tr]] <- read.csv(paste(here("Data", "BTE", "INLA", "Event_"), as.character(tr), ".csv", sep="")) %>% 
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


#TODO: Compare running times of rescaled with unscaled.

#Create survival objects
Surv.list <- vector("list", transitions[length(transitions)])
for(tr in transitions) { 
  Surv.list[[tr]] <- inla.surv(time = event.list.ext[[tr]]$time, event = event.list.ext[[tr]]$status)
}

### Run INLA ###

inla.test1 <- joint(formSurv=list(inla.surv(event.list.ext[[31]]$time, event.list.ext[[31]]$status) ~ LATIT + LONGI),
                    basRisk = rep("weibullsurv", 1),
                    dataSurv = list(event.list.ext[[31]]))

summary(inla.test1)
plot(inla.test1)$Baseline
plot(inla.test1)$Outcomes

inla.test2 <- joint(formSurv=list(inla.surv(event.list.ext[[31]]$time, event.list.ext[[31]]$status) ~ cov_Tmean),
                    basRisk = rep("weibullsurv", 1),
                    dataSurv = list(event.list.ext[[31]]))

summary(inla.test2)
plot(inla.test2)$Baseline
plot(inla.test2)$Outcomes

inla.test3 <- joint(formSurv=list(inla.surv(event.list.ext[[31]]$time, event.list.ext[[31]]$status) ~ cov_Tmean + cov_CMI),
                    basRisk = rep("weibullsurv", 1),
                    dataSurv = list(event.list.ext[[31]]))

summary(inla.test3)
plot(inla.test3)$Baseline
plot(inla.test3)$Outcomes

inla.test4 <- joint(formSurv=list(inla.surv(event.list.ext[[31]]$time, event.list.ext[[31]]$status) ~ cov_Tmean_scaled + cov_CMI_scaled),
                    basRisk = rep("weibullsurv", 1),
                    dataSurv = list(event.list.ext[[31]]))

summary(inla.test4)
plot(inla.test4)$Baseline
plot(inla.test4)$Outcomes


inla3.scaled # is the one with scaled covariates Tmean and CMI
inla3.scaled <- lapply(c(31, 32, 34:39), function(x) {
  print(x)
  joint(formSurv= list(inla.surv(time, status) ~ cov_Tmean_scaled + cov_CMI_scaled),
        basRisk = rep("weibullsurv", 1),
        dataSurv = list(event.list.ext[[x]]))})

for(i in 1:8) {print(summary(inla3.scaled[[i]]))}

inla3 <- lapply(c(31, 32, 34:39), function(x) {
  print(x)
  joint(formSurv= list(inla.surv(time, status) ~ cov_Tmean + cov_CMI),
        basRisk = rep("weibullsurv", 1),
        dataSurv = list(event.list.ext[[x]]))})


for(i in 1:8) {print(summary(inla3[[i]]))}


summary(inla3[[1]])
summary(inla3[[2]])
summary(inla3[[3]])
summary(inla3[[4]])
summary(inla3[[5]])
summary(inla3[[6]])
summary(inla3[[7]])
summary(inla3[[8]])
plot(inla3[[8]])$Baseline
plot(inla3[[8]])$Outcomes



inla.test5 <- joint(formSurv=list(inla.surv(event.list.ext[[31]]$time, event.list.ext[[31]]$status) ~ cov_Tmean + cov_CMI,
                                  inla.surv(event.list.ext[[32]]$time, event.list.ext[[32]]$status) ~ cov_Tmean + cov_CMI),
                    basRisk = rep("weibullsurv", 2),
                    dataSurv = list(event.list.ext[[31]],
                                    event.list.ext[[32]]))

summary(inla.test5)
plot(inla.test5)$Baseline
plot(inla.test5)$Outcomes


inla.test6 <- joint(formSurv=list(inla.surv(event.list.ext[[31]]$time, event.list.ext[[31]]$status) ~ cov_Tmean + cov_CMI,
                                  inla.surv(event.list.ext[[32]]$time, event.list.ext[[32]]$status) ~ cov_Tmean + cov_CMI,
                                  inla.surv(event.list.ext[[34]]$time, event.list.ext[[34]]$status) ~ cov_Tmean + cov_CMI,
                                  inla.surv(event.list.ext[[35]]$time, event.list.ext[[35]]$status) ~ cov_Tmean + cov_CMI,
                                  inla.surv(event.list.ext[[36]]$time, event.list.ext[[36]]$status) ~ cov_Tmean + cov_CMI,
                                  inla.surv(event.list.ext[[37]]$time, event.list.ext[[37]]$status) ~ cov_Tmean + cov_CMI,
                                  inla.surv(event.list.ext[[38]]$time, event.list.ext[[38]]$status) ~ cov_Tmean + cov_CMI,
                                  inla.surv(event.list.ext[[39]]$time, event.list.ext[[39]]$status) ~ cov_Tmean + cov_CMI),
                    basRisk = rep("weibullsurv", 8),
                    dataSurv = list(event.list.ext[[31]],
                                    event.list.ext[[32]],
                                    event.list.ext[[34]],
                                    event.list.ext[[35]],
                                    event.list.ext[[36]],
                                    event.list.ext[[37]],
                                    event.list.ext[[38]],
                                    event.list.ext[[39]]))

summary(inla.test6)
plot(inla.test6)$Baseline
plot(inla.test6)$Outcomes


inla.test7 <- joint(formSurv=list(inla.surv(event.list.ext[[31]]$time, event.list.ext[[31]]$status) ~ cov_soil),
                    basRisk = rep("weibullsurv", 1),
                    dataSurv = list(event.list.ext[[31]]))

summary(inla.test7)
plot(inla.test7)$Baseline
plot(inla.test7)$Outcomes

inla.test8 <- joint(formSurv=list(inla.surv(event.list.ext[[31]]$time, event.list.ext[[31]]$status) ~ cov_soil,
                                  inla.surv(event.list.ext[[32]]$time, event.list.ext[[32]]$status) ~ cov_soil,
                                  inla.surv(event.list.ext[[34]]$time, event.list.ext[[34]]$status) ~ cov_soil),
                    basRisk = rep("weibullsurv", 3),
                    dataSurv = list(event.list.ext[[31]],
                                    event.list.ext[[32]],
                                    event.list.ext[[34]]))

summary(inla.test8)
plot(inla.test8)$Baseline
plot(inla.test8)$Outcomes$S3

inla.test9 <- joint(formSurv=list(inla.surv(event.list.ext[[31]]$time, event.list.ext[[31]]$status) ~ cov_soil,
                                  inla.surv(event.list.ext[[32]]$time, event.list.ext[[32]]$status) ~ cov_soil,
                                  inla.surv(event.list.ext[[34]]$time, event.list.ext[[34]]$status) ~ cov_soil,
                                  inla.surv(event.list.ext[[35]]$time, event.list.ext[[35]]$status) ~ cov_soil,
                                  inla.surv(event.list.ext[[36]]$time, event.list.ext[[36]]$status) ~ cov_soil,
                                  inla.surv(event.list.ext[[37]]$time, event.list.ext[[37]]$status) ~ cov_soil,
                                  inla.surv(event.list.ext[[38]]$time, event.list.ext[[38]]$status) ~ cov_soil,
                                  inla.surv(event.list.ext[[39]]$time, event.list.ext[[39]]$status) ~ cov_soil),
                    basRisk = rep("weibullsurv", 8),
                    dataSurv = list(event.list.ext[[31]],
                                    event.list.ext[[32]],
                                    event.list.ext[[34]],
                                    event.list.ext[[35]],
                                    event.list.ext[[36]],
                                    event.list.ext[[37]],
                                    event.list.ext[[38]],
                                    event.list.ext[[39]]))

summary(inla.test9)
plot(inla.test9)$Baseline
plot(inla.test9)$Outcomes

##Interpretation:

t <- seq(0.1, 1000, by = 1)
riskW <- function(t, lambda, alpha) lambda*alpha*t^(alpha-1)

risk1 <- riskW(t, exp(inla.test8$summary.fixed["Intercept_S1", "mean"]),
               inla.test8$summary.hyperpar$mean[1])

risk1s <- riskW(t, exp(inla.test8$summary.fixed["Intercept_S1", "mean"] +
                         inla.test8$summary.fixed["cov_soil_S1", "mean"]),
                inla.test8$summary.hyperpar$mean[1])


plot(t, risk1, pch=19)


### Smaller database

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

inla.31.T.C <- joint(formSurv=list(inla.surv(event.list.ext[[31]]$time, event.list.ext[[31]]$status) ~ cov_Tmean + cov_CMI),
                   basRisk = rep("weibullsurv", 1),
                   dataSurv = list(event.list.ext[[31]]))

summary(inla.31.T.C)
plot(inla.31.T.C)$Baseline
plot(inla.31.T.C)$Outcomes

inla.31.T.C.S <- joint(formSurv=list(inla.surv(event.list.ext[[31]]$time, event.list.ext[[31]]$status) ~ cov_Tmean + cov_CMI + cov_soil),
                     basRisk = rep("weibullsurv", 1),
                     dataSurv = list(event.list.ext[[31]]))

summary(inla.31.T.C.S)
plot(inla.31.T.C.S)$Baseline
plot(inla.31.T.C.S)$Outcomes

inla.31.T.C.S.sc <- joint(formSurv=list(inla.surv(event.list.ext[[31]]$time, event.list.ext[[31]]$status) ~ cov_Tmean_scaled + cov_CMI_scaled + cov_soil),
                       basRisk = rep("weibullsurv", 1),
                       dataSurv = list(event.list.ext[[31]]))

summary(inla.31.T.C.S.sc)
plot(inla.31.T.C.S.sc)$Baseline
plot(inla.31.T.C.S.sc)$Outcomes

#Create list for each individual transition: 
inla3.T.C.S <- vector("list", 8)
inla3.T.C.S <- lapply(c(31, 32, 34:39), function(x) {
  print(x)
  joint(formSurv= list(inla.surv(time, status) ~ cov_Tmean + cov_CMI + cov_soil),
        basRisk = rep("weibullsurv", 1),
        dataSurv = list(event.list.ext[[x]]))})


#Now try how many we can combine into one joint. 

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



