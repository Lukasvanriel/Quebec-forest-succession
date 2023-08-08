###TODO: insert description

### Load packages ###
library(here)

#### Load data ####
bte1_noP <- read.csv(here("Data", "BTE", "bte1_noP.csv"))[,-1]
bte2_noP <- read.csv(here("Data", "BTE", "bte2_noP.csv"))[,-1]
bte3_noP <- read.csv(here("Data", "BTE", "bte3_noP.csv"))[,-1]
bte4_noP <- read.csv(here("Data", "BTE", "bte4_noP.csv"))[,-1]
bte5_noP <- read.csv(here("Data", "BTE", "bte5_noP.csv"))[,-1]

#Climate

#TODO: add 

#Use random values for now
library(random) #Remove
bte1_noP$cov_Tmean <- sample.int(20, nrow(bte1_noP), replace = TRUE)
bte2_noP$cov_Tmean <- sample.int(20, nrow(bte2_noP), replace = TRUE)
bte3_noP$cov_Tmean <- sample.int(20, nrow(bte3_noP), replace = TRUE)
bte4_noP$cov_Tmean <- sample.int(20, nrow(bte4_noP), replace = TRUE)
bte5_noP$cov_Tmean <- sample.int(20, nrow(bte5_noP), replace = TRUE)

bte1_noP$cov_CMI <- sample.int(20, nrow(bte1_noP), replace = TRUE)
bte2_noP$cov_CMI <- sample.int(20, nrow(bte2_noP), replace = TRUE)
bte3_noP$cov_CMI <- sample.int(20, nrow(bte3_noP), replace = TRUE)
bte4_noP$cov_CMI <- sample.int(20, nrow(bte4_noP), replace = TRUE)
bte5_noP$cov_CMI <- sample.int(20, nrow(bte5_noP), replace = TRUE)

#Soil

bte1_noP$cov_soil <- str_sub(bte1_noP$STRATE_SIFORT, -3, -2)
bte2_noP$cov_soil <- str_sub(bte2_noP$STRATE_SIFORT, -3, -2)
bte3_noP$cov_soil <- str_sub(bte3_noP$STRATE_SIFORT, -3, -2)
bte4_noP$cov_soil <- str_sub(bte4_noP$STRATE_SIFORT, -3, -2)
bte5_noP$cov_soil <- str_sub(bte5_noP$STRATE_SIFORT, -3, -2)

#Perturbations

# Function that returns a class and severity depending on the input perturbation string 
determine_perturb_class <- function(pert_string) {
  if (is.na(pert_string)) {
    pert_class <- 0
    pert_sev <- 0
  } else if(pert_string %in% c("BR", "BRD", "BRU")) {
      pert_class <- 1
      pert_sev <- 2 # A "1" can be added later if we want to look dat partial perturbations
    } else if (pert_string %in% c("CBA","CBT","CEF","CPT","CRB","CRS","CS","CT","ETR","RPS")) {
      pert_class <- 2
      pert_sev <- 2 # A "1" can be added later if we want to look dat partial perturbations
    } else {
      pert_class <- 0
      pert_sev <- 0
    }
  c(pert_class, pert_sev)
}

perturbation.matrix1 <- t(sapply(bte1_noP$ORIGINE, determine_perturb_class))
perturbation.matrix2 <- t(sapply(bte2_noP$ORIGINE, determine_perturb_class))
perturbation.matrix3 <- t(sapply(bte3_noP$ORIGINE, determine_perturb_class))
perturbation.matrix4 <- t(sapply(bte4_noP$ORIGINE, determine_perturb_class))
perturbation.matrix5 <- t(sapply(bte5_noP$ORIGINE, determine_perturb_class))

# Add to dataset
bte1_noP$cov_pert_class <- factor(perturbation.matrix1[,1])
bte2_noP$cov_pert_class <- factor(perturbation.matrix2[,1])
bte3_noP$cov_pert_class <- factor(perturbation.matrix3[,1])
bte4_noP$cov_pert_class <- factor(perturbation.matrix4[,1])
bte5_noP$cov_pert_class <- factor(perturbation.matrix5[,1])

bte1_noP$cov_pert_sev <- factor(perturbation.matrix1[,2])
bte2_noP$cov_pert_sev <- factor(perturbation.matrix2[,2])
bte3_noP$cov_pert_sev <- factor(perturbation.matrix3[,2])
bte4_noP$cov_pert_sev <- factor(perturbation.matrix4[,2])
bte5_noP$cov_pert_sev <- factor(perturbation.matrix5[,2])

### Also update the combined dataset ###
bte_all <- rbind(bte1_noP, bte2_noP, bte3_noP, bte4_noP, bte5_noP)

###There is a significant fraction of TESSELLE that are lacking Species information
#TODO: Look into which to drop and which to keep
bte <- bte_all |> filter(! is.na(GR_ESS))

#### Write out resulting datasets ####
write.csv(bte1_noP, here("Data", "BTE", "bte1_noP_cov.csv"))
write.csv(bte2_noP, here("Data", "BTE", "bte2_noP_cov.csv"))
write.csv(bte3_noP, here("Data", "BTE", "bte3_noP_cov.csv"))
write.csv(bte4_noP, here("Data", "BTE", "bte4_noP_cov.csv"))
write.csv(bte5_noP, here("Data", "BTE", "bte5_noP_cov.csv"))

write.csv(bte, here("Data", "BTE", "bte_all_cov.csv"))
