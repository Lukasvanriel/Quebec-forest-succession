#### Runs the markov chain model using the INLA method ####

#### Load packages ####
library(tidyverse)
library(INLA)
library(INLAjoint)

### Data ###
#Changes to dataset can be made in Msm.R; Or need to structure it differently
data_msm <- read.csv(here("Data", "BTE", "bte_msm_ready.csv"))[,-1]

### Prepare dataset ###

data.inla <- data_msm |>
  select(TESSELLE, time, sp_class, LONGI, LATIT) |>
  mutate(ID=TESSELLE) |>
  mutate(from=sp_class) |>
  mutate(to=sp_class) |>
  mutate(entry=time) |>
  relocate(time, .after = last_col())

ph.from <- c(0, data.inla$to)
data.inla$from <- ph.from[1:(length(ph.from)-1)]

ph.entry <- c(0, data.inla$time)
data.inla$entry <- ph.entry[1:(length(ph.entry)-1)] 

data.inla.ready <- data.inla |>
  filter(time>0) |>
  select(-sp_class) |>
  mutate(trans=to-from) |>
  filter(trans != 0) |> 
  ungroup() |>
  select(-TESSELLE, -trans)

## normalise latitute and longitude

data.inla.ready <- data.inla.ready |>
  mutate(lon=(LONGI-mean(LONGI))/sd(LONGI)) |>
  mutate(lat=(LATIT-mean(LATIT))/sd(LATIT))

### Run INLA ###

##State table:
st <- matrix(0, ncol=9, nrow=9)

for(i in 1:nrow(data.inla.ready)) {
  st[data.inla.ready$from[i], data.inla.ready$to[i]] <- st[data.inla.ready$from[i], data.inla.ready$to[i]] + 1
}
st

#TODO: Add actual INLA models