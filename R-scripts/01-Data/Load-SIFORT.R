#### Script to extract and write down SIFORT data ####

#### Load packages ####
library(tidyverse)
library(sf)
library(here)

#### Read in raw data ####
setwd(dirname(getActiveDocumentContext()$path))

data_raw1 <- sf::st_read("../../../Data/Raw/TES_PRG_1_DV.gdb.zip") |> 
  st_drop_geometry()
data_raw2 <- sf::st_read("../../../Data/Raw/TES_PRG_2_DV.gdb.zip") |> 
  st_drop_geometry()
data_raw3 <- sf::st_read("../../../Data/Raw/TES_PRG_3_DV.gdb.zip") |> 
  st_drop_geometry()
data_raw4 <- sf::st_read("../../../Data/Raw/TES_PRG_4_DV.gdb.zip") |> 
  st_drop_geometry()
data_raw5 <- sf::st_read("../../../Data/Raw/TES_PRG_5_V0.gdb.zip") |> 
  st_drop_geometry()

#### Clean data ####
#Drop rows without forest type attribution
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

rm(data_raw1, data_raw2, data_raw3, data_raw4, data_raw5)

#Fix problem of swapped longitude and latitude values
data1[data1$LONGI > 0, c("LATIT", "LONGI")] <-  data1[data1$LONGI > 0, c("LONGI", "LATIT")]
data2[data2$LONGI > 0, c("LATIT", "LONGI")] <-  data2[data2$LONGI > 0, c("LONGI", "LATIT")]  
data3[data3$LONGI > 0, c("LATIT", "LONGI")] <-  data3[data3$LONGI > 0, c("LONGI", "LATIT")]  
data4[data4$LONGI > 0, c("LATIT", "LONGI")] <-  data4[data4$LONGI > 0, c("LONGI", "LATIT")]  
data5[data5$LONGI > 0, c("LATIT", "LONGI")] <-  data5[data5$LONGI > 0, c("LONGI", "LATIT")]  

#### Filter data for Boreal-Temperate ecotone ####
filter_bte <- function(raw.data, strict) {
  if (strict) {
    domains <- c("4est", "4ouest")
  } else {domains <- c("3est", "3ouest", "4est", "4ouest", "5est", "5ouest")}
  raw.data |> filter(SDOM_ECO %in% domains)
}

bte1 <- filter_bte(data1, strict=TRUE)
bte2 <- filter_bte(data2, strict=TRUE)
bte3 <- filter_bte(data3, strict=TRUE)
bte4 <- filter_bte(data4, strict=TRUE)
bte5 <- filter_bte(data5, strict=TRUE)

sum(is.na(bte4$GR_ESS))
head(bte4[is.na(bte4$GR_ESS),], 20)

plot(bte4[]$LONGI, bte4[]$LATIT)
points(bte4[is.na(bte4$GR_ESS),]$LONGI, bte4[is.na(bte4$GR_ESS),]$LATIT, col="red", pch=0.1)
table(bte4[is.na(bte4$GR_ESS),]$GR_ESS)

#### Filter out plantations ####
#TODO
bte1_noP <- bte1 
bte2_noP <- bte2 
bte3_noP <- bte3
bte4_noP <- bte4 
bte5_noP <- bte5 

#### Write out resulting datasets ####
write.csv(bte1_noP, "../../../Data/BTE/bte1_noP.csv")
write.csv(bte2_noP, "../../../Data/BTE/bte2_noP.csv")
write.csv(bte3_noP, "../../../Data/BTE/bte3_noP.csv")
write.csv(bte4_noP, "../../../Data/BTE/bte4_noP.csv")
write.csv(bte5_noP, "../../../Data/BTE/bte5_noP.csv")

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

# Combine all
bte_all <- rbind(bte1_noP, bte2_noP, bte3_noP, bte4_noP, bte5_noP)

###There is a significant fraction of TESSELLE that are lacking Species information
#TODO: Look into which to drop and which to keep
bte <- bte_all |> filter(! is.na(GR_ESS))

#### Write out combined data ####
write.csv(bte, "../../../Data/BTE/bte_all.csv")
