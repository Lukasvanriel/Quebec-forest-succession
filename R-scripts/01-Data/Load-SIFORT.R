#### Script to extract and write down SIFORT data ####

#### Load packages ####
library(tidyverse)
library(sf)
library(here)

#### Read in raw data ####
data_raw1 <- sf::st_read(here("Data", "Raw", "TES_PRG_1_DV.gdb.zip")) |> 
  st_drop_geometry()
data_raw2 <- sf::st_read(here("Data", "Raw", "TES_PRG_2_DV.gdb.zip")) |> 
  st_drop_geometry()
data_raw3 <- sf::st_read(here("Data", "Raw", "TES_PRG_3_DV.gdb.zip")) |> 
  st_drop_geometry()
data_raw4 <- sf::st_read(here("Data", "Raw", "TES_PRG_4_DV.gdb.zip")) |> 
  st_drop_geometry()
data_raw5 <- sf::st_read(here("Data", "Raw", "TES_PRG_5_V0.gdb.zip")) |> 
  st_drop_geometry()

#### Clean data ####
#Drop rows without data
data1 <- data_raw1 |> 
  filter(! is.na(data_raw1$LATIT))
data2 <- data_raw2 |> 
  filter(! is.na(data_raw2$LATIT))
data3 <- data_raw3 |> 
  filter(! is.na(data_raw3$LATIT))
data4 <- data_raw4 |> 
  filter(! is.na(data_raw4$LATIT))
data5 <- data_raw5 |> 
  filter(! is.na(data_raw5$LATIT))

#rm(data_raw1, data_raw2, data_raw3, data_raw4, data_raw5)

#Fix problem of swapped longitude and latitude values
data1[data1$LONGI > 0, c("LATIT", "LONGI")] <-  data1[data1$LONGI > 0, c("LONGI", "LATIT")]
data2[data2$LONGI > 0, c("LATIT", "LONGI")] <-  data2[data2$LONGI > 0, c("LONGI", "LATIT")]  
data3[data3$LONGI > 0, c("LATIT", "LONGI")] <-  data3[data3$LONGI > 0, c("LONGI", "LATIT")]  
data4[data4$LONGI > 0, c("LATIT", "LONGI")] <-  data4[data4$LONGI > 0, c("LONGI", "LATIT")]  
data5[data5$LONGI > 0, c("LATIT", "LONGI")] <-  data5[data5$LONGI > 0, c("LONGI", "LATIT")]  

#### Filter data for Boreal-Temperate ecotone ####
filter_bte <- function(data, strict) {
  if (strict) {
    domains <- c("4est", "4ouest")
  } else {domains <- c("3est", "3ouest", "4est", "4ouest", "5est", "5ouest")}
  data |> filter(SDOM_ECO %in% domains)
}

bte1 <- filter_bte(data1, strict=TRUE)
bte2 <- filter_bte(data2, strict=TRUE)
bte3 <- filter_bte(data3, strict=TRUE)
bte4 <- filter_bte(data4, strict=TRUE)
bte5 <- filter_bte(data5, strict=TRUE)

rm(data1, data2, data3, data4, data5)

#### Filter out plantations ####
# Load documentation file on species groups

docu <- read.csv(here("Data", "Documentation", "Groupement-essences.csv"))[,-1] |>
  rename(GR_ESS = Code)

# There are many codes with double or even triple descriptions
docu_filter <- docu[!duplicated(docu[,"GR_ESS"]),]

#TODO: Find better way to filter out multiple occurrences
#Try to look a bit more in detail:
# docu_filt <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), colnames(docu))
# for (ess in names(table(docu$GR_ESS))) {
#   docu_subset <- docu[docu$GR_ESS==ess,]
#   if(nrow(docu_subset) == 1) {
#     docu_filt <- rbind(docu_filt, docu_subset)
#   } else if(nrow(docu_subset) == 2) {
#     if(docu_subset[1,2] == docu_subset[2,2]) {
#       docu_filt <- rbind(docu_filt, docu_subset[1,])
#     } else(docu_filt <- rbind(docu_filt, docu_subset))
#   } else if(nrow(docu_subset) == 3) {
#     
#   }
# }

#Function to filter out all Tesselle that are plantations
filt_plantations <- function(data, docu) {
  ##Based on disturbance information:
  tess_origine <- data |>
    filter(ORIGINE %in% c("P", "PA", "PE", "PL", "PLB", "PLN", "PLR", "PRR")) #TODO
  
  ##Based on Description of species class:
  data_doc <- plyr::join(data, docu, by="GR_ESS", type="left")
  #Test if the word "plantation/Plantation" is written in the description
  tess_doc <- data_doc |> 
    mutate(Plant=grepl("Plantation", Description)) |>
    mutate(plant=grepl("plantation", Description)) |>
    filter(Plant | plant)
  
  ##Create vector by combining both
  tess_union <- c(tess_doc$TESSELLE, 
                  tess_origine$TESSELLE[! tess_origine$TESSELLE %in% tess_doc$TESSELLE])
  
  data_noP <- data |>
    filter(! TESSELLE %in% tess_union)
}

bte1_noP <- filt_plantations(bte1, docu_filter)
bte2_noP <- filt_plantations(bte2, docu_filter)
bte3_noP <- filt_plantations(bte3, docu_filter)
bte4_noP <- filt_plantations(bte4, docu_filter)
bte5_noP <- filt_plantations(bte5, docu_filter)


### Combine all data in one dataframe
# Only keep columns that are present in bte5
bte1_noP <- bte1_noP %>%
  select(intersect(colnames(bte1_noP) , colnames(bte5_noP)))
bte2_noP <- bte2_noP %>%
  select(intersect(colnames(bte2_noP) , colnames(bte5_noP)))
bte3_noP <- bte3_noP %>%
  select(intersect(colnames(bte3_noP) , colnames(bte5_noP)))
bte4_noP <- bte4_noP %>%
  select(intersect(colnames(bte4_noP) , colnames(bte5_noP)))

# Merge all 5
bte_all <- rbind(bte1_noP, bte2_noP, bte3_noP, bte4_noP, bte5_noP)

###There is a significant fraction of TESSELLE that are lacking Species information
#TODO: Look into which to drop and which to keep
bte <- bte_all |> filter(! is.na(GR_ESS))

#### Write out resulting datasets ####
write.csv(bte1_noP, here("Data", "BTE", "bte1_noP.csv"))
write.csv(bte2_noP, here("Data", "BTE", "bte2_noP.csv"))
write.csv(bte3_noP, here("Data", "BTE", "bte3_noP.csv"))
write.csv(bte4_noP, here("Data", "BTE", "bte4_noP.csv"))
write.csv(bte5_noP, here("Data", "BTE", "bte5_noP.csv"))

write.csv(bte, here("Data", "BTE", "bte_all.csv"))
