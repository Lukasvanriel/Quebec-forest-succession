#### Load packages ####
library(tidyverse)
library(INLA)
library(INLAjoint)
library(here)
library(conflicted)
library(msm)

conflicts_prefer(dplyr::filter)

### Data ###
#Changes to dataset can be made in Msm.R; Or need to structure it differently
data_msm2 <- read.csv(here("Data", "BTE", "bte_msm_ready.csv"))[,-1]

length(unique(data_msm2$TESSELLE))

# How many TESSELLE do we have?

data1 <- read.csv(here("Data", "Raw", "SIFORT", "raw_data1.csv"))
data2 <- read.csv(here("Data", "Raw", "SIFORT", "raw_data2.csv"))
data3 <- read.csv(here("Data", "Raw", "SIFORT", "raw_data3.csv"))
data4 <- read.csv(here("Data", "Raw", "SIFORT", "raw_data4.csv"))
data5 <- read.csv(here("Data", "Raw", "SIFORT", "raw_data5.csv"))

length(unique(c(data1$TESSELLE, data2$TESSELLE, data3$TESSELLE, data4$TESSELLE, data5$TESSELLE)))

bteE1 <- data1 %>% filter(SDOM_ECO %in% c("3est", "3ouest", "4est", "4ouest", "5est", "5ouest"))
bteE2 <- data2 %>% filter(SDOM_ECO %in% c("3est", "3ouest", "4est", "4ouest", "5est", "5ouest"))
bteE3 <- data3 %>% filter(SDOM_ECO %in% c("3est", "3ouest", "4est", "4ouest", "5est", "5ouest"))
bteE4 <- data4 %>% filter(SDOM_ECO %in% c("3est", "3ouest", "4est", "4ouest", "5est", "5ouest"))
bteE5 <- data5 %>% filter(SDOM_ECO %in% c("3est", "3ouest", "4est", "4ouest", "5est", "5ouest"))

length(unique(c(bteE1$TESSELLE, bteE2$TESSELLE, bteE3$TESSELLE, bteE4$TESSELLE, bteE5$TESSELLE)))

bte1 <- data1 %>% filter(SDOM_ECO %in% c("4est", "4ouest"))
bte2 <- data2 %>% filter(SDOM_ECO %in% c("4est", "4ouest"))
bte3 <- data3 %>% filter(SDOM_ECO %in% c("4est", "4ouest"))
bte4 <- data4 %>% filter(SDOM_ECO %in% c("4est", "4ouest"))
bte5 <- data5 %>% filter(SDOM_ECO %in% c("4est", "4ouest"))

length(unique(c(bte1$TESSELLE, bte2$TESSELLE, bte3$TESSELLE, bte4$TESSELLE, bte5$TESSELLE)))

## So ideally we could extrapolate up to 0.5M, 1.6M and 4.4M
# 9 classes, so 72 transitions
# 3 continuous variables: 2 climatic + 1 soil
# 1 categorical variable with 4 levels (perturbation)
# 1 add. continuous variables (e.g. time since last pert?)
# 1 add. categorical variable (e.g. perturbation intensity or...)
# 1 frailty cov. (e.g. domain)

### Let's create the datasets ####
# Choose 3 appropriate subregions
table(data_msm2$SREG_ECO)

aT <- data_msm2 %>% filter(SREG_ECO == "4aT")
bT <- data_msm2 %>% filter(SREG_ECO == "4cT")
cT <- data_msm2 %>% filter(SREG_ECO == "4bT")


length(unique(aT$TESSELLE))
length(unique(bT$TESSELLE))
length(unique(bT$TESSELLE))


statetable.msm(sp_class, TESSELLE, aT) # !
statetable.msm(sp_class, TESSELLE, bT) # !
statetable.msm(sp_class, TESSELLE, cT) # !

datasets <- list(aT, bT, cT)
datasets.filt <- lapply(datasets, function(x) {
  x %>% 
    mutate(cov_frail = factor(-round(LONGI,0)),
           TESSELLE = as.character(TESSELLE),
           #sp_class = factor(sp_class),
           cov_pert_class = factor(cov_pert_class),
           cov_pert_sev = factor(cov_pert_sev)) %>% 
    select(TESSELLE, sp_class, time, cov_Tmean, cov_CMI, cov_soil, cov_time_pert, cov_pert_class, cov_pert_sev, cov_frail) %>%
    group_by(TESSELLE) %>%
    filter(n() > 1) %>%
    ungroup()
  })

select_random_TESSELLE <- function(df, N){
  Tsel <- sample(unique(df$TESSELLE), N)  # Randomly sample N unique individuals
  df %>% filter(TESSELLE %in% Tsel)
}

datasets_1k <- lapply(datasets.filt, function(x) select_random_TESSELLE(x, 1000))
datasets_2k <- lapply(datasets.filt, function(x) select_random_TESSELLE(x, 2000))
datasets_4k <- lapply(datasets.filt, function(x) select_random_TESSELLE(x, 4000))
datasets_10k <- lapply(datasets.filt, function(x) select_random_TESSELLE(x, 10000))
datasets_20k <- lapply(datasets.filt, function(x) select_random_TESSELLE(x, 20000))

### Prep datasets for INLA:

expand.trans <- function(line.fr, line.to, transitions = 1:9){
  
  expansions <- lapply(transitions[-line.fr$sp_class], FUN = function(x) {
    lines <- data.frame(ID = line.fr$TESSELLE,
                        Tstart = line.fr$time, 
                        Tstop = line.to$time,
                        from = line.fr$sp_class,
                        to = x,
                        Status = ifelse(x == line.to$sp_class, 1, 0),
                        cov_Tmean = line.fr$cov_Tmean,
                        cov_CMI = line.fr$cov_CMI,
                        cov_soil = line.fr$cov_soil,
                        cov_time_pert = line.fr$cov_time_pert,
                        cov_pert_class = line.fr$cov_pert_class,
                        cov_pert_sev = line.fr$cov_pert_sev,
                        cov_frail = line.fr$cov_frail) 
  } )
  bind_rows(expansions)
}

expand.stay <- function(line.fr, line.to, future = -1, transitions = 1:9){
  inclu <- ifelse(future == -1, list(transitions[-line.fr$sp_class]), list(future))[[1]]
  expansions <- lapply(inclu, FUN = function(x) {
    lines <- data.frame(ID = line.fr$TESSELLE,
                        Tstart = line.fr$time, 
                        Tstop = line.to$time,
                        from = line.fr$sp_class,
                        to = x,
                        Status = 0,
                        cov_Tmean = line.fr$cov_Tmean,
                        cov_CMI = line.fr$cov_CMI,
                        cov_soil = line.fr$cov_soil,
                        cov_time_pert = line.fr$cov_time_pert,
                        cov_pert_class = line.fr$cov_pert_class,
                        cov_pert_sev = line.fr$cov_pert_sev,
                        cov_frail = line.fr$cov_frail) 
  } )
  bind_rows(expansions)
}

expand.all <- function(block){
  expansion <- lapply(2:nrow(block), FUN = function(x) {
    if(block$sp_class[x] == block$sp_class[x - 1]) {
      no.fut.change <- all(block$sp_class[x:nrow(block)] == block$sp_class[x])
      if(no.fut.change) {
        if(x == nrow(block)){
          line.cor <- block[x - 1,]
          r = 2
          while((x - r) > 0 && block$sp_class[x - 1] == block$sp_class[x - r]) {
            line.cor <- block[x - r,]
            r <- r + 1
          }
          expand.stay(line.cor, block[x,], future = -1)}
      } else{
        first.fut <- block$sp_class[x:nrow(block)][which(block$sp_class[x:nrow(block)] != block$sp_class[x])[1]]
        expand.stay(block[x - 1,], block[x,], first.fut) # Only for future state because rest will be taken care of then as well
      }
    } else {
      line.cor <- block[x - 1,]
      r = 2
      while((x - r) > 0 && block$sp_class[x - 1] == block$sp_class[x - r]) {
        line.cor <- block[x - r,]
        r <- r + 1
      }
      expand.trans(line.cor, block[x,])
    }
  })
  bind_rows(expansion)
}

## The data
data_inla <- bind_rows(lapply(unique(datasets_1k[[3]]$TESSELLE), FUN = function(x) {
  expand.all(datasets_1k[[3]] %>% filter(TESSELLE == x))
} ))

datasets_1k_inla <- lapply(datasets_1k, function(df) {
  bind_rows(lapply(unique(df$TESSELLE), FUN = function(x) {
    expand.all(df %>% filter(TESSELLE == x))
  } ))
})

datasets_2k_inla <- lapply(datasets_2k, function(df) {
  bind_rows(lapply(unique(df$TESSELLE), FUN = function(x) {
    expand.all(df %>% filter(TESSELLE == x))
  } ))
})

datasets_4k_inla <- lapply(datasets_4k, function(df) {
  bind_rows(lapply(unique(df$TESSELLE), FUN = function(x) {
    expand.all(df %>% filter(TESSELLE == x))
  } ))
})

datasets_10k_inla <- lapply(datasets_10k, function(df) {
  bind_rows(lapply(unique(df$TESSELLE), FUN = function(x) {
    expand.all(df %>% filter(TESSELLE == x))
  } ))
})

datasets_20k_inla <- lapply(datasets_20k, function(df) {
  bind_rows(lapply(unique(df$TESSELLE), FUN = function(x) {
    expand.all(df %>% filter(TESSELLE == x))
  } ))
})

write_rds(datasets_1k_inla, here("Data", "HardwareRequirements", "HR_inla_1k.RDS"))
write_rds(datasets_2k_inla, here("Data", "HardwareRequirements", "HR_inla_2k.RDS"))
write_rds(datasets_4k_inla, here("Data", "HardwareRequirements", "HR_inla_4k.RDS"))
write_rds(datasets_10k_inla, here("Data", "HardwareRequirements", "HR_inla_10k.RDS"))
write_rds(datasets_20k_inla, here("Data", "HardwareRequirements", "HR_inla_20k.RDS"))

######### Run the models:

data.inla <- read_rds(here("Data", "HardwareRequirements", "HR_inla_2k.RDS"))

df1 <- data.inla[[2]]

## Convert to survival objects:
get_transition_states <- function(k, N) {
  from_state <- (k - 1) %/% (N - 1) + 1
  position <- (k - 1) %% (N - 1)
  
  to_states <- setdiff(1:N, from_state)  # Exclude from-state
  to_state <- to_states[position + 1]   # Pick the corresponding to-state
  
  return(c(from_state, to_state))
}

Nb.states <- 9
Surv.list <- vector("list", Nb.states * (Nb.states - 1))

event.list <- lapply(seq_along(Surv.list), FUN = function(x){
  state.info <- get_transition_states(x, Nb.states)
  df1 %>% 
    filter(from == state.info[1],
           to == state.info[2])
})

for(i in seq_along(Surv.list)) { 
  Surv.list[[i]] <- inla.surv(time = event.list[[i]]$Tstop,
                              truncation = event.list[[i]]$Tstart,
                              event = event.list[[i]]$Status)
}


#### Now set up the runs

# Automate cause otherwise have to repeat 72 times

names(Surv.list) <- paste0("s", seq_along(Surv.list))

# Assign named survival objects to the environment
list2env(Surv.list, envir = .GlobalEnv)


formulas <- lapply(seq_along(Surv.list), function(i) {
  as.formula(paste0("s", i, " ~ cov_Tmean"))
  } )

weib.surv <- joint(formSurv = inla.surv(time = event.list[[1]]$Tstop,
                                        truncation = event.list[[1]]$Tstart,
                                        event = event.list[[1]]$Status) ~ 1,
                   basRisk = rep("exponentialsurv", 1), dataSurv = event.list[[1]],
                   control = list(config=TRUE))
summary(weib.surv)
plot(weib.surv)


weib.surv <- joint(formSurv = formulas,
                   basRisk = rep("weibullsurv", length(formulas)), dataSurv = event.list,
                   control = list(config=TRUE))

s1 <- Surv.list[[1]]
s2 <- Surv.list[[2]]
s3 <- Surv.list[[3]]
s4 <- Surv.list[[4]]
s5 <- Surv.list[[5]]
s6 <- Surv.list[[6]]
s7 <- Surv.list[[7]]
s8 <- Surv.list[[8]]
s9 <- Surv.list[[9]]
s10 <- Surv.list[[10]]
s11 <- Surv.list[[11]]
s12 <- Surv.list[[12]]
s13 <- Surv.list[[13]]
s14 <- Surv.list[[14]]
s15 <- Surv.list[[15]]
s16 <- Surv.list[[16]]
s17 <- Surv.list[[17]]
s18 <- Surv.list[[18]]
s19 <- Surv.list[[19]]
s20 <- Surv.list[[20]]
s21 <- Surv.list[[21]]
s22 <- Surv.list[[22]]
s23 <- Surv.list[[23]]
s24 <- Surv.list[[24]]
s25 <- Surv.list[[25]]
s26 <- Surv.list[[26]]
s27 <- Surv.list[[27]]
s28 <- Surv.list[[28]]
s29 <- Surv.list[[29]]
s30 <- Surv.list[[30]]
s31 <- Surv.list[[31]]
s32 <- Surv.list[[32]]
s33 <- Surv.list[[33]]
s34 <- Surv.list[[34]]
s35 <- Surv.list[[35]]
s36 <- Surv.list[[36]]
s37 <- Surv.list[[37]]
s38 <- Surv.list[[38]]
s39 <- Surv.list[[39]]
s40 <- Surv.list[[40]]
s41 <- Surv.list[[41]]
s42 <- Surv.list[[42]]
s43 <- Surv.list[[43]]
s44 <- Surv.list[[44]]
s45 <- Surv.list[[45]]
s46 <- Surv.list[[46]]
s47 <- Surv.list[[47]]
s48 <- Surv.list[[48]]
s49 <- Surv.list[[49]]
s50 <- Surv.list[[50]]
s51 <- Surv.list[[51]]
s52 <- Surv.list[[52]]
s53 <- Surv.list[[53]]
s54 <- Surv.list[[54]]
s55 <- Surv.list[[55]]
s56 <- Surv.list[[56]]
s57 <- Surv.list[[57]]
s58 <- Surv.list[[58]]
s59 <- Surv.list[[59]]
s60 <- Surv.list[[60]]
s61 <- Surv.list[[61]]
s62 <- Surv.list[[62]]
s63 <- Surv.list[[63]]
s64 <- Surv.list[[64]]
s65 <- Surv.list[[65]]
s66 <- Surv.list[[66]]
s67 <- Surv.list[[67]]
s68 <- Surv.list[[68]]
s69 <- Surv.list[[69]]
s70 <- Surv.list[[70]]
s71 <- Surv.list[[71]]
s72 <- Surv.list[[72]]

weib.surv <- joint(formSurv = list(
  s1 ~ cov_Tmean + cov_pert_class,
  s2 ~ cov_Tmean + cov_pert_class,
  s3 ~ cov_Tmean + cov_pert_class,
  s4 ~ cov_Tmean + cov_pert_class,
  s5 ~ cov_Tmean + cov_pert_class,
  s6 ~ cov_Tmean + cov_pert_class,
  s7 ~ cov_Tmean + cov_pert_class,
  s8 ~ cov_Tmean + cov_pert_class,
  s9 ~ cov_Tmean + cov_pert_class,
  s10 ~ cov_Tmean + cov_pert_class,
  s11 ~ cov_Tmean + cov_pert_class,
  s12 ~ cov_Tmean + cov_pert_class,
  s13 ~ cov_Tmean + cov_pert_class,
  s14 ~ cov_Tmean + cov_pert_class,
  s15 ~ cov_Tmean + cov_pert_class,
  s16 ~ cov_Tmean + cov_pert_class,
  s17 ~ cov_Tmean + cov_pert_class,
  s18 ~ cov_Tmean + cov_pert_class,
  s19 ~ cov_Tmean + cov_pert_class,
  s20 ~ cov_Tmean + cov_pert_class,
  s21 ~ cov_Tmean + cov_pert_class,
  s22 ~ cov_Tmean + cov_pert_class,
  s23 ~ cov_Tmean + cov_pert_class,
  s24 ~ cov_Tmean + cov_pert_class,
  s25 ~ cov_Tmean + cov_pert_class,
  s26 ~ cov_Tmean + cov_pert_class,
  s27 ~ cov_Tmean + cov_pert_class,
  s28 ~ cov_Tmean + cov_pert_class,
  s29 ~ cov_Tmean + cov_pert_class,
  s30 ~ cov_Tmean + cov_pert_class,
  s31 ~ cov_Tmean + cov_pert_class,
  s32 ~ cov_Tmean + cov_pert_class,
  s33 ~ cov_Tmean + cov_pert_class,
  s34 ~ cov_Tmean + cov_pert_class,
  s35 ~ cov_Tmean + cov_pert_class,
  s36 ~ cov_Tmean + cov_pert_class,
  s37 ~ cov_Tmean + cov_pert_class,
  s38 ~ cov_Tmean + cov_pert_class,
  s39 ~ cov_Tmean + cov_pert_class,
  s40 ~ cov_Tmean + cov_pert_class,
  s41 ~ cov_Tmean + cov_pert_class,
  s42 ~ cov_Tmean + cov_pert_class,
  s43 ~ cov_Tmean + cov_pert_class,
  s44 ~ cov_Tmean + cov_pert_class,
  s45 ~ cov_Tmean + cov_pert_class,
  s46 ~ cov_Tmean + cov_pert_class,
  s47 ~ cov_Tmean + cov_pert_class,
  s48 ~ cov_Tmean + cov_pert_class,
  s49 ~ cov_Tmean + cov_pert_class,
  s50 ~ cov_Tmean + cov_pert_class,
  s51 ~ cov_Tmean + cov_pert_class,
  s52 ~ cov_Tmean + cov_pert_class,
  s53 ~ cov_Tmean + cov_pert_class,
  s54 ~ cov_Tmean + cov_pert_class,
  s55 ~ cov_Tmean + cov_pert_class,
  s56 ~ cov_Tmean + cov_pert_class,
  s57 ~ cov_Tmean + cov_pert_class,
  s58 ~ cov_Tmean + cov_pert_class,
  s59 ~ cov_Tmean + cov_pert_class,
  s60 ~ cov_Tmean + cov_pert_class,
  s61 ~ cov_Tmean + cov_pert_class,
  s62 ~ cov_Tmean + cov_pert_class,
  s63 ~ cov_Tmean + cov_pert_class,
  s64 ~ cov_Tmean + cov_pert_class,
  s65 ~ cov_Tmean + cov_pert_class,
  s66 ~ cov_Tmean + cov_pert_class,
  s67 ~ cov_Tmean + cov_pert_class,
  s68 ~ cov_Tmean + cov_pert_class,
  s69 ~ cov_Tmean + cov_pert_class,
  s70 ~ cov_Tmean + cov_pert_class,
  s71 ~ cov_Tmean + cov_pert_class,
  s72 ~ cov_Tmean + cov_pert_class
),
                   basRisk = rep("exponentialsurv", 72), dataSurv = event.list,
                   control = list(config=TRUE))
summary(weib.surv)
plot(weib.surv)


sum(is.na(event.list[[x]]$cov_pert_class))

lapply(event.list, function(x) sum(is.na(x$cov_pert_class)))

