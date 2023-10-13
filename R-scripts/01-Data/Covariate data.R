###TODO: insert description

### Load packages ###
library(tidyverse)
library(here)
library(data.table)

#### Load data ####
bte1_noP <- read.csv(here("Data", "BTE", "bte1_noP.csv"))[,-1]
bte2_noP <- read.csv(here("Data", "BTE", "bte2_noP.csv"))[,-1]
bte3_noP <- read.csv(here("Data", "BTE", "bte3_noP.csv"))[,-1]
bte4_noP <- read.csv(here("Data", "BTE", "bte4_noP.csv"))[,-1]
bte5_noP <- read.csv(here("Data", "BTE", "bte5_noP.csv"))[,-1]

### Create covariate columns ###
##Climate
sg <- readRDS(here("Data", "Raw", "Bioclim", "sg_complete.rds"))
cmi <- readRDS(here("Data", "Raw", "Bioclim", "cmi_complete.rds"))

cmi.tot <- cmi |>
  select(cmi_sum, year, TESSELLE) |>
  mutate(year=as.numeric(year)) |>
  mutate(TESSELLE=as.numeric(TESSELLE)) |>
  rename(cov_CMI=cmi_sum)

sg.tot <- sg |>
  select(an_meanT, year, TESSELLE) |>
  mutate(year=as.numeric(year)) |>
  mutate(TESSELLE=as.numeric(TESSELLE)) |>
  rename(cov_Tmean=an_meanT)

clim.tot <- left_join(cmi.tot, sg.tot)

### Filter out the weird years ###
table(bte1_noP$AN_PRO_SOU)
table(bte2_noP$AN_PRO_SOU)
table(bte3_noP$AN_PRO_SOU)
table(bte4_noP$AN_PRO_SOU)
table(bte5_noP$AN_PRO_SOU)

bte4_noP$AN_PRO_SOU[! bte4_noP$AN_PRO_SOU > 1997]

bte1_noP <- bte1_noP[! is.na(bte1_noP$AN_PRO_SOU),]
bte2_noP <- bte2_noP[! is.na(bte2_noP$AN_PRO_SOU),]
bte3_noP <- bte3_noP[! is.na(bte3_noP$AN_PRO_SOU),]
bte4_noP <- bte4_noP[! is.na(bte4_noP$AN_PRO_SOU),]
bte4_noP <- bte4_noP |>
  filter(AN_PRO_SOU > 1997)

### Extract climatic values for each measurement ###
### TODO: Collect 10 year averages ###

# Subset for faster extractions
clim.tot1 <- clim.tot |>
  filter(year >= (min(bte1_noP$AN_PRO_SOU)), year <= max(bte1_noP$AN_PRO_SOU))
clim.tot2 <- clim.tot |>
  filter(year >= (min(bte2_noP$AN_PRO_SOU)), year <= max(bte2_noP$AN_PRO_SOU))
clim.tot3 <- clim.tot |>
  filter(year >= (min(bte3_noP$AN_PRO_SOU)), year <= max(bte3_noP$AN_PRO_SOU))
clim.tot4 <- clim.tot |>
  filter(year >= (min(bte4_noP$AN_PRO_SOU)), year <= max(bte4_noP$AN_PRO_SOU))
clim.tot5 <- clim.tot |>
  filter(year >= (min(bte5_noP$AN_PRO_SOU)), year <= max(bte5_noP$AN_PRO_SOU))

# Convert to data tables to speed up
clim1.dt <- as.data.table(clim.tot1)
clim2.dt <- as.data.table(clim.tot2)
clim3.dt <- as.data.table(clim.tot3)
clim4.dt <- as.data.table(clim.tot4)
clim5.dt <- as.data.table(clim.tot5)

setkey(clim1.dt, TESSELLE, year)
setkey(clim2.dt, TESSELLE, year)
setkey(clim3.dt, TESSELLE, year)
setkey(clim4.dt, TESSELLE, year)
setkey(clim5.dt, TESSELLE, year)

bte1_noP_dt <- as.data.table(bte1_noP)
bte1_noP_dt <- rename(bte1_noP_dt, year=AN_PRO_SOU)
bte2_noP_dt <- as.data.table(bte2_noP)
bte2_noP_dt <- rename(bte2_noP_dt, year=AN_PRO_SOU)
bte3_noP_dt <- as.data.table(bte3_noP)
bte3_noP_dt <- rename(bte3_noP_dt, year=AN_PRO_SOU)
bte4_noP_dt <- as.data.table(bte4_noP)
bte4_noP_dt <- rename(bte4_noP_dt, year=AN_PRO_SOU)
bte5_noP_dt <- as.data.table(bte5_noP)
bte5_noP_dt <- rename(bte5_noP_dt, year=AN_PRO_SOU)

setkey(bte1_noP_dt, TESSELLE, year)
setkey(bte2_noP_dt, TESSELLE, year)
setkey(bte3_noP_dt, TESSELLE, year)
setkey(bte4_noP_dt, TESSELLE, year)
setkey(bte5_noP_dt, TESSELLE, year)

# Join

bte1_noP_dt <- clim1.dt[bte1_noP_dt, on = c("TESSELLE", "year")]
bte2_noP_dt <- clim2.dt[bte2_noP_dt, on = c("TESSELLE", "year")]
bte3_noP_dt <- clim3.dt[bte3_noP_dt, on = c("TESSELLE", "year")]
bte4_noP_dt <- clim4.dt[bte4_noP_dt, on = c("TESSELLE", "year")]
bte5_noP_dt <- clim5.dt[bte5_noP_dt, on = c("TESSELLE", "year")]

# library(microbenchmark)
# microbenchmark(result.filter <- apply(bte4_noP[1:100,], 1, extract_cmi_filter),
#                result.sub <-  apply(bte4_noP[1:100,], 1, extract_cmi_subset), 
#                result.table <- apply(bte4_noP[1:100,], 1, extract_cmi_dt),
#                times=10)
# 
# start <- Sys.time()
# t <- apply(bte4_noP[1:100,], 1, extract_cmi_dt)
# print( Sys.time() - start )
####
# 
# # Subset for faster extractions
# cmi.tot1 <- cmi.tot |>
#   filter(year >= (min(bte1_noP$AN_PRO_SOU)-10), year < max(bte1_noP$AN_PRO_SOU))
# cmi.tot2 <- cmi.tot |>
#   filter(year >= (min(bte2_noP$AN_PRO_SOU)-10), year < max(bte2_noP$AN_PRO_SOU))
# cmi.tot3 <- cmi.tot |>
#   filter(year >= (min(bte3_noP$AN_PRO_SOU)-10), year < max(bte3_noP$AN_PRO_SOU))
# cmi.tot4 <- cmi.tot |>
#   filter(year >= (min(bte4_noP$AN_PRO_SOU)-10), year < max(bte4_noP$AN_PRO_SOU))
# cmi.tot5 <- cmi.tot |>
#   filter(year >= (min(bte5_noP$AN_PRO_SOU)-10), year < max(bte5_noP$AN_PRO_SOU))
# 
# 
# 
# test <- cmi.tot4 |>
#   filter(TESSELLE == bte4_noP$TESSELLE[1], as.numeric(year) < bte4_noP$AN_PRO_SOU[1], as.numeric(year) >= bte4_noP$AN_PRO_SOU[1] - 10) |>
#   summarise(av.cmi=mean(cmi_sum)) |>
#   pull(av.cmi)
# 
# cmi_vector4 <- rep(0, nrow(bte1_noP))
# for (i in 1:nrow(bte1_noP)) {
#   print(i)
#   cmi_vector4[i] <- mean(cmi.tot4$cmi_sum[as.logical((cmi.tot4$TESSELLE == bte4_noP$TESSELLE[i]) * (cmi.tot4$year >= (bte4_noP$AN_PRO_SOU[1]-10)) * (cmi.tot4$year < bte4_noP$AN_PRO_SOU[i]))])
# }
# for (i in 1:nrow(bte1_noP)) {
#   print(i)
#   cmi_vector4[i] <- cmi.tot4 |>
#     filter(TESSELLE == bte4_noP$TESSELLE[i], as.numeric(year) < bte4_noP$AN_PRO_SOU[i], as.numeric(year) >= bte4_noP$AN_PRO_SOU[1] - 10) |>
#     summarise(av.cmi=mean(cmi_sum)) |>
#     pull(av.cmi)
# }
# 
# cmi_vector4[70]
# mean(cmi.tot4$cmi_sum[as.logical((cmi.tot4$TESSELLE == bte4_noP$TESSELLE[1]) * (cmi.tot4$year >= (bte4_noP$AN_PRO_SOU[1]-10)) * (cmi.tot4$year < bte4_noP$AN_PRO_SOU[1]))])
# 
# (cmi.tot4$TESSELLE == bte4_noP$TESSELLE[1])
# 

##Soil
extract_soil_cov <- function(df) {
  df$number_only <- as.integer(gsub("[^0-9]", "", df$STRATE_SIFORT)) %% 10
}

bte1_noP_dt$cov_soil <- extract_soil_cov(bte1_noP_dt)
bte2_noP_dt$cov_soil <- extract_soil_cov(bte2_noP_dt)
bte3_noP_dt$cov_soil <- extract_soil_cov(bte3_noP_dt)
bte4_noP_dt$cov_soil <- extract_soil_cov(bte4_noP_dt)
bte5_noP_dt$cov_soil <- extract_soil_cov(bte5_noP_dt)

#Ensure soil class stays the same over time



##Perturbations
#TODO: fix missing NA's in perturbation dates

time_since <- function(pert, meas) {
  if(is.na(pert)) {
    return(-1)
    } else {return(meas - pert)}
}

# Function that returns a class and severity depending on the input perturbation string 
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
time_since(NA, 2000)

perturbation.matrix1 <- t(mapply(determine_perturb_class, bte1_noP_dt$ORIGINE, 
                                 bte1_noP_dt$AN_ORIGINE, bte1_noP_dt$year))
perturbation.matrix2 <- t(mapply(determine_perturb_class, bte2_noP_dt$ORIGINE, 
                                 bte2_noP_dt$AN_ORIGINE, bte2_noP_dt$year))
perturbation.matrix3 <- t(mapply(determine_perturb_class, bte3_noP_dt$ORIGINE, 
                                 bte3_noP_dt$AN_ORIGINE, bte3_noP_dt$year))
perturbation.matrix4 <- t(mapply(determine_perturb_class, bte4_noP_dt$ORIGINE, 
                                 bte4_noP_dt$AN_ORIGINE, bte4_noP_dt$year))
perturbation.matrix5 <- t(mapply(determine_perturb_class, bte5_noP_dt$ORIGINE, 
                                 bte5_noP_dt$AN_ORIGINE, bte5_noP_dt$year))

# Add to dataset
bte1_noP_dt$cov_pert_class <- factor(perturbation.matrix1[,1])
bte2_noP_dt$cov_pert_class <- factor(perturbation.matrix2[,1])
bte3_noP_dt$cov_pert_class <- factor(perturbation.matrix3[,1])
bte4_noP_dt$cov_pert_class <- factor(perturbation.matrix4[,1])
bte5_noP_dt$cov_pert_class <- factor(perturbation.matrix5[,1])

bte1_noP_dt$cov_pert_sev <- factor(perturbation.matrix1[,2])
bte2_noP_dt$cov_pert_sev <- factor(perturbation.matrix2[,2])
bte3_noP_dt$cov_pert_sev <- factor(perturbation.matrix3[,2])
bte4_noP_dt$cov_pert_sev <- factor(perturbation.matrix4[,2])
bte5_noP_dt$cov_pert_sev <- factor(perturbation.matrix5[,2])

bte1_noP_dt$cov_time_pert <- perturbation.matrix1[,3]
bte2_noP_dt$cov_time_pert <- perturbation.matrix2[,3]
bte3_noP_dt$cov_time_pert <- perturbation.matrix3[,3]
bte4_noP_dt$cov_time_pert <- perturbation.matrix4[,3]
bte5_noP_dt$cov_time_pert <- perturbation.matrix5[,3]

sum(is.na(bte4_noP_dt$AN_ORIGINE))
sum(is.na(bte4_noP_dt$ORIGINE))

### Also update the combined dataset ###
bte_all <- rbind(bte1_noP_dt, bte2_noP_dt, bte3_noP_dt, bte4_noP_dt, bte5_noP_dt)

###There is a significant fraction of TESSELLE that are lacking Species information
#TODO: Look into which to drop and which to keep
bte <- bte_all |> filter(! is.na(GR_ESS))

#### Write out resulting datasets ####
write.csv(bte1_noP_dt, here("Data", "BTE", "bte1_noP_cov.csv"))
write.csv(bte2_noP_dt, here("Data", "BTE", "bte2_noP_cov.csv"))
write.csv(bte3_noP_dt, here("Data", "BTE", "bte3_noP_cov.csv"))
write.csv(bte4_noP_dt, here("Data", "BTE", "bte4_noP_cov.csv"))
write.csv(bte5_noP_dt, here("Data", "BTE", "bte5_noP_cov.csv"))

write.csv(bte, here("Data", "BTE", "bte_all_cov.csv"))
