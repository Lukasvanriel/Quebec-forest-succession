#### Script to extract and write down SIFORT data ####

#### Load packages ####
library(tidyverse)
library(sf)
library(here)

#TODO: Add disclaimer about date accessed + website
#### Read in raw data ####
data_raw1 <- sf::st_read(here("Data", "Raw", "SIFORT", "TES_PRG_1_DV.gdb.zip")) |> 
  st_drop_geometry()
data_raw2 <- sf::st_read(here("Data", "Raw", "SIFORT", "TES_PRG_2_DV.gdb.zip")) |> 
  st_drop_geometry()
data_raw3 <- sf::st_read(here("Data", "Raw", "SIFORT", "TES_PRG_3_DV.gdb.zip")) |> 
  st_drop_geometry()
data_raw4 <- sf::st_read(here("Data", "Raw", "SIFORT", "TES_PRG_4_DV.gdb.zip")) |> 
  st_drop_geometry()
data_raw5 <- sf::st_read(here("Data", "Raw", "SIFORT", "TES_PRG_5_V0.gdb.zip")) |> 
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

empty1 <- data_raw1 |> 
  filter(is.na(data_raw1$LATIT))
empty2 <- data_raw2 |> 
  filter(is.na(data_raw2$LATIT))
empty3 <- data_raw3 |> 
  filter(is.na(data_raw3$LATIT))
empty4 <- data_raw4 |> 
  filter(is.na(data_raw4$LATIT))
empty5 <- data_raw5 |> 
  filter(is.na(data_raw5$LATIT))

plot(empty1$LONGI, empty1$LATIT)
table(empty1$TESSELLE)
empty1$TESSELLE == empty2$TESSELLE

#rm(data_raw1, data_raw2, data_raw3, data_raw4, data_raw5)
write.csv(data1, here("Data", "Raw", "SIFORT", "raw_data1.csv"))
write.csv(data2, here("Data", "Raw", "SIFORT", "raw_data2.csv"))
write.csv(data3, here("Data", "Raw", "SIFORT", "raw_data3.csv"))
write.csv(data4, here("Data", "Raw", "SIFORT", "raw_data4.csv"))
write.csv(data5, here("Data", "Raw", "SIFORT", "raw_data5.csv"))

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

#rm(data1, data2, data3, data4, data5)
#rm(data_raw1, data_raw2, data_raw3, data_raw4, data_raw5)

###There is a significant fraction of TESSELLE that are lacking Species information

bte1 <- bte1 |> filter(! is.na(GR_ESS))
bte2 <- bte2 |> filter(! is.na(GR_ESS))
bte3 <- bte3 |> filter(! is.na(GR_ESS))
bte4 <- bte4 |> filter(! is.na(GR_ESS))
bte5 <- bte5 |> filter(! is.na(GR_ESS))

# 
bte11 <- bte1 |> filter(is.na(GR_ESS))
bte22 <- bte2 |> filter(is.na(GR_ESS))
bte33 <- bte3 |> filter(is.na(GR_ESS))
bte44 <- bte4 |> filter(is.na(GR_ESS))
bte55 <- bte5 |> filter(is.na(GR_ESS))


#TODO: Look into which to drop and which to keep
#Split into tesselle that never have a species description and those that miss some
tess.list <- names(table(c(bte11$TESSELLE, bte22$TESSELLE, bte33$TESSELLE,
                           bte44$TESSELLE, bte55$TESSELLE)))

#Create boolean vector indicating which may be saved based on other information
tess.salv <- sapply(tess.list, FUN = function(t) {
  all.info <- c(bte1[bte1$TESSELLE == t, "GR_ESS"], bte2[bte2$TESSELLE == t, "GR_ESS"],
                     bte3[bte3$TESSELLE == t, "GR_ESS"], bte4[bte4$TESSELLE == t, "GR_ESS"],
                     bte5[bte5$TESSELLE == t, "GR_ESS"])
  if(sum(! is.na(all.info)) > 0){return(TRUE)
  } else {return(FALSE)}
  })
tess.salv

### Check out documentation file on species groups ### 

docu <- read.csv(here("Data", "Documentation", "Groupement-essences.csv")) |>
  rename(GR_ESS = Code)

# There are many codes with double or even triple descriptions
docu_filter <- docu[!duplicated(docu[,"GR_ESS"]),]

#TODO: Find better way to filter out multiple occurrences
codes_present <- names(table(c(bte1$GR_ESS, bte2$GR_ESS, bte3$GR_ESS, bte4$GR_ESS, bte5$GR_ESS)))
docu_mult <- docu |>
  filter(GR_ESS %in% codes_present) |>
  count(GR_ESS) |>
  filter(n > 1)

docu_mult_diff <- as.data.frame(matrix(nrow=0, ncol=2))
for(i in 1:length(docu_mult$GR_ESS)) {
  descr <- docu$Description[docu$GR_ESS == docu_mult$GR_ESS[i]]
  if(length(descr) == 2) {
    if(substring(descr[2], 1, 11) == "(BFEC 3-4) "){descr[2] <- substring(descr[2], 12, nchar(descr[2]))}
    if(! substring(descr[1], 1, 11) == "Plantation " && substring(descr[2], 1, 11) == "Plantation "){
      if(descr[1] != descr[2]) {docu_mult_diff <- rbind(docu_mult_diff, docu[docu$GR_ESS == docu_mult$GR_ESS[i],])}
    }
  } else {docu_mult_diff <- rbind(docu_mult_diff, docu[docu$GR_ESS == docu_mult$GR_ESS[i],])}
}

#Go over all possibly problematic codes
for(i in 1:length(table(docu_mult_diff$GR_ESS))) {
  print(names(table(docu_mult_diff$GR_ESS))[i])
  print(docu$Description[docu$GR_ESS == names(table(docu_mult_diff$GR_ESS))[i]])
  
  # Wait for Enter key to be pressed before continuing
  cat("Press Enter to continue to the next iteration...")
  readline()
}

risk_codes <- c("EPEU", "ERR", "ES", "FIPB", "FIPG", "FIPR", "FTPB", "MEPG", 
                "PEPB", "PEPG", "PEPR", "PGFI", "PGME", "PGPE", "PISE")

for(c in risk_codes){
  print(c)
  print(docu[docu$GR_ESS == c,])
  # Wait for Enter key to be pressed before continuing
  cat("Press Enter to continue to the next iteration...")
  readline()
}
# a <- filter(bte5, GR_ESS == "PISE")

## Correct documentation for multiple occurences
#The first occurence is the correct for almost all group codes
docu_filter <- docu[!duplicated(docu[,"GR_ESS"]),]

#Only for PISE group this needs to be corrected
docu_filter[docu_filter$GR_ESS == "PISE",] <- docu[docu$GR_ESS == "PISE",][2,]

#Several codes need to be dropped from 4 & 5 since then they became plantations
bte4 <- bte4 |>
  filter(! GR_ESS %in% risk_codes[-c(1,2,3,15)])
bte5 <- bte5 |>
  filter(! GR_ESS %in% risk_codes[-c(1,2,3,15)])

# ES contained 2 different descriptions, so best is to rename them
#docu[docu$GR_ESS == "ES",]
#docu[docu$GR_ESS == "EBEB",]
#docu[docu$GR_ESS == "ER",]

bte1$GR_ESS[bte1$GR_ESS == "ES"] <- "EBEB"
bte2$GR_ESS[bte2$GR_ESS == "ES"] <- "EBEB"
bte3$GR_ESS[bte3$GR_ESS == "ES"] <- "EBEB"
bte4$GR_ESS[bte4$GR_ESS == "ES"] <- "ER"
bte5$GR_ESS[bte5$GR_ESS == "ES"] <- "ER"

#### Filter out plantations ####
table(c(bte1$ORIGINE, bte2$ORIGINE, bte3$ORIGINE, bte4$ORIGINE, bte5$ORIGINE))

#Function to filter out all Tesselle that are plantations
filt_plantations <- function(data, docu) {
  ##Based on disturbance information:
  tess_origine <- data |>
    filter(ORIGINE %in% c("P", "PA", "PE", "PL", "PLB", "PLN", "PLR", "PRR",
                          "ENM", "ENS", "ETR", "FR", "REA", "RPS", "VG")) #TODO: Maybe too many e.g FR?
  
  ##Based on Description of species class:
  data_doc <- plyr::join(data, docu, by="GR_ESS", type="left")
  #Test if the word "plantation/Plantation" is written in the description #TODO: Other words to filter on?
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

docu_filter$Description
docu_noP <- docu_filter |> 
  mutate(Plant=grepl("Plantation", Description)) |>
  mutate(plant=grepl("plantation", Description)) |>
  filter(! Plant | plant) |>
  select(Description)

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

# Merge all
bte <- rbind(bte1_noP, bte2_noP, bte3_noP, bte4_noP, bte5_noP)

#### Write out resulting datasets ####
write.csv(bte1_noP, here("Data", "BTE", "bte1_noP.csv"))
write.csv(bte2_noP, here("Data", "BTE", "bte2_noP.csv"))
write.csv(bte3_noP, here("Data", "BTE", "bte3_noP.csv"))
write.csv(bte4_noP, here("Data", "BTE", "bte4_noP.csv"))
write.csv(bte5_noP, here("Data", "BTE", "bte5_noP.csv"))

write.csv(bte, here("Data", "BTE", "bte_all.csv"))
