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

data.inla.inter <- data_msm |>
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

event.list <- list()
for(i in transitions) {
  event.list[[i]] <- data.frame(matrix(ncol = ncol(data.inla)))
  event.list[[i]]$status <- c(NA)
  colnames(event.list[[i]]) <- c(colnames(data.inla), "status")
  event.list[[i]] <- event.list[[i]][-1,-c(1,2)]
}

# Now fill in the full event dataframes for each transition
for(tr in transitions[-c(1,2)]) {
  print(tr)
  data.tr <- data.inla |>
    filter(from == tr %/% 10)
  
  tr.event <- vector("list", length=nrow(data.tr))
  for(i in 1:nrow(data.tr)) { 
    if(data.tr$to[i] == tr %% 10) {
      tr.event[[i]] <- c(data.tr[i,-c(1,2)], 1)
      names(tr.event[[i]])[length(tr.event[[i]])] <- "status"
    } else {
      tr.event[[i]] <- c(data.tr[i,-c(1,2)], 0)
      tr.event[[i]]$to <- tr %% 10
      names(tr.event[[i]])[length(tr.event[[i]])] <- "status"
    }
  }
  event.list[[tr]] <- bind_rows(tr.event)
  write.csv(event.list[[tr]], 
            paste(here("Data", "BTE", "INLA", "Event_"), as.character(tr), ".csv", sep=""))
}

##Read in event datasets
for(tr in transitions) {
  event.list[[tr]] <- read.csv(paste(here("Data", "BTE", "INLA", "Event_"), as.character(tr), ".csv", sep=""))
}

#Create survival objects
Surv10 <- list()
for(tr in transitions) { 
  Surv10[[tr]] <- inla.surv(time = event.list[[tr]]$time, event = event.list[[tr]]$status)
}

### Run INLA ###

a12 <- Surv10[[12]]
a13 <- Surv10[[13]]
inla.12 <- joint(formSurv=list(a12 ~ lat + lon),
                     basRisk = rep("weibullsurv", 1),
                     dataSurv = list(event.list[[12]]))
summary(inla.12)
plot(inla.12)$Baseline
plot(inla.12)$Outcomes

inla.12.13 <- joint(formSurv=list(a12 ~ lat + lon,
                                  a13 ~ lat + lon),
                 basRisk = rep("weibullsurv", 2),
                 dataSurv = list(event10[[12]],
                                 event10[[13]]))



