#### Extract bioclimatic variables 2km2 resolution #####


# bioclimatic raster data were obtained from Natural Resources Canada: ftp://ftp.nrcan.gc.ca/pub/outgoing/NAM_grids
# This script extract bioclimatic for all Qc forest plots
# Return the mean of the bioclimatic for the year of the plot survey as well as over the last 10 years before the year of the plot survey.

### PACKAGES ####

library(raster)
library(sf)
library(zoo)
library(reshape2)
library(tidyr)
library(dplyr)
library(here)
library(data.table)
library(tidytable)


### DATA ####
data_raw4 <- sf::st_read(here("Data", "Raw", "TES_PRG_4_DV.gdb.zip"))

bte1 <- filter_bte(data_raw4, strict=TRUE)
bte2 <- filter_bte(data_raw4, strict=TRUE)
bte3 <- filter_bte(data_raw4, strict=TRUE)
bte4 <- filter_bte(data_raw4, strict=TRUE)
bte5 <- filter_bte(data_raw4, strict=TRUE)

bte1[bte1$LONGI > 0, c("LATIT", "LONGI")] <-  bte1[bte1$LONGI > 0, c("LONGI", "LATIT")]
bte2[bte2$LONGI > 0, c("LATIT", "LONGI")] <-  bte2[bte2$LONGI > 0, c("LONGI", "LATIT")]  
bte3[bte3$LONGI > 0, c("LATIT", "LONGI")] <-  bte3[bte3$LONGI > 0, c("LONGI", "LATIT")]  
bte4[bte4$LONGI > 0, c("LATIT", "LONGI")] <-  bte4[bte4$LONGI > 0, c("LONGI", "LATIT")]  
bte5[bte5$LONGI > 0, c("LATIT", "LONGI")] <-  bte5[bte5$LONGI > 0, c("LONGI", "LATIT")]  

#Did not filter on plantations yet!
rm(data_raw4, data4, bte4_noP, docu, docu_filter)

bte1 <- bte1 %>% 
  st_transform("+proj=longlat +datum=WGS84 +no_defs")
bte2 <- bte2 %>% 
  st_transform("+proj=longlat +datum=WGS84 +no_defs")
bte3 <- bte3 %>% 
  st_transform("+proj=longlat +datum=WGS84 +no_defs")
bte4 <- bte4 %>% 
  st_transform("+proj=longlat +datum=WGS84 +no_defs")
bte5 <- bte5 %>% 
  st_transform("+proj=longlat +datum=WGS84 +no_defs")

new5 <- bte5 |>
  filter(! TESSELLE %in% bte1$TESSELLE)
  

##Functions

retrieveClimateData <- function(years = 1900:2022,
                                info =  c("bio", "cmi", "mint", "maxt", "pcp", "sg"), res = 300,
                                path = here("Data", "Raw", "Bioclim"), geom) {
  
  stopifnot(res %in% c(60, 300))
  ls_clim <- list()
  
  dir.create(path, showWarnings = FALSE)
  
  #natural-resources.canada.ca
  # basurl <- "ftp://ftp.nrcan.gc.ca/pub/outgoing/NAM_grids/zipfiles"
  basurl <- "https://ftp.maps.canada.ca/pub/nrcan_rncan/Climate-archives_Archives-climatologiques/NAM_monthly/monthly_by_var"
  info <- match.arg(info)
  beg <- paste0(basurl, "/")
  end <- paste0("_", res, "arcsec.zip")
  # year available: from 1900 to 2018
  for (year in years) {
    tmp <- tempfile(fileext = ".zip")
    print(tmp)
    print(paste0(beg, info, year, end))
    
    # Download
    curl::curl_download(paste0(beg, info, year, end), tmp)
    unzip(tmp, exdir = path)
    unlink(tmp)
    
    # extract data
    ls_tmp <- extract_climate_data(path = path, info = info,
                                   year = year, geom = geom)
    
    # Save intermediate results just in case it crashes
    saveRDS(ls_tmp, paste0(path, "/", info, year, "_", res, ".rds"))
    ls_clim[[paste0(info, year)]] <- ls_tmp
    unlink(paste0(path, "/", year), recursive = TRUE)
  }
  invisible(NULL)
  
  # Save final results
  saveRDS(ls_clim, paste0(path, "/", info, "_", res, "_extra.rds"))}

### Function to extract data from raster to multipoints ####

# Modified from https://github.com/inSileco/inSilecoDataRetrieval/blob/master/R/get_climate_nam_grids.R

extract_climate_data <- function(path, info, year, geom, pattern = "\\.asc$|\\.tif$") {
  nm_fo <- paste0(path, "/", year)
  fls <- list.files(nm_fo, pattern = pattern, full.names = TRUE)
  
  out <- lapply(lapply(fls, raster), 
                function(x) extract(crop(x, y =  geom), y =  geom))
  
  names(out) <- paste0(gsub(list.files(nm_fo, pattern = pattern), pat = pattern, rep = ""), "_", year)
  
  out
}

###
retrieveClimateData(years = 2021:2022, info = "sg", res = 60, geom = bte4, path = here("Data", "Raw", "Bioclim"))
retrieveClimateData(years = 2019:2022, info = "cmi", res = 60, geom = bte4, path = here("Data", "Raw", "Bioclim"))

helper_sg <- list()
helper_cmi <- list()

for(y in 1960:2018) {
  print(y)
  helper_sg[[as.character(y)]] <- readRDS(here("Data", "Raw", "Bioclim", paste0("sg", as.character(y), "_60.rds")))
  helper_cmi[[as.character(y)]] <- readRDS(here("Data", "Raw", "Bioclim", paste0("cmi", as.character(y), "_60.rds")))
}
for(y in 2019:2020) {
  print(y)
  helper_sg[[as.character(y)]] <- readRDS(here("Data", "Raw", "Bioclim", paste0("sg", as.character(y), "_60.rds")))
}

saveRDS(helper_sg, here("Data", "Raw", "Bioclim", "sg_60.rds"))
saveRDS(helper_cmi, here("Data", "Raw", "Bioclim", "cmi_60.rds"))


### Read in bioclimatic variables ###
sg <- readRDS(here("Data", "Raw", "Bioclim", "sg_60.rds"))
cmi <- readRDS(here("Data", "Raw", "Bioclim", "cmi_60.rds"))

####
# Find NA values (from plots that are located at the margin of the climate raster)
na_sg <- which(is.na(sg[[1]]$sg60_01_1960))
na_cmi <- which(is.na(cmi[[1]]$cmi60_sum_1960))
na_sg_xy <- bte4[na_sg,]
na_cmi_xy <- bte4[na_cmi,]

####
#########
####

# Replace NAs with nearest neighbor

xy_unassign_sg <- st_join(na_sg_xy, bte4[-na_sg,], join = st_nearest_feature)
xy_unassign_cmi <- st_join(na_cmi_xy, bte4[-na_cmi,], join = st_nearest_feature)

unassign_sg <- unlist(lapply(xy_unassign_sg$TESSELLE.y, function(x) which(bte4$TESSELLE %in% x)))
unassign_cmi <- unlist(lapply(xy_unassign_cmi$TESSELLE.y, function(x) which(bte4$TESSELLE %in% x)))

sg <- rapply(sg, function(x) replace(x, na_sg, x[unassign_sg]), how = 'list')
cmi <- rapply(cmi, function(x) replace(x, na_cmi, x[unassign_cmi]), how = 'list')

### From list to df ####
#cmi_backup <- cmi
cmi <- cmi_backup

class(cmi[[1]])
cmi[[1]]

names(cmi) <- 1960:2018
cmi <- lapply(cmi, as.data.frame)
cmi <- lapply(cmi, function(x) {rownames(x) <- bte4$TESSELLE; x})
cmi <- lapply(cmi, function(x) {colnames(x) <- paste0("cmi_", c(1:12, 'sum')); x})
cmi <- do.call(rbind, lapply(cmi, data.frame))
cmi <- tibble::rownames_to_column(cmi)
cmi <- cmi %>%
  tidyr::separate(rowname, c("year", "TESSELLE"), "\\.")

# Convert sg to a tidytable
cmi_ttbl <- as_tidytable(cmi)

# Separate the "rowname" column
cmi_ttbl <- cmi_ttbl %>%
  separate(rowname, into = c("year", "TESSELLE"), sep = "\\.", remove = FALSE)

# Convert back to a data frame
cmi.test <- as.data.frame(cmi_ttbl)

saveRDS(cmi, here("Data", "Raw", "Bioclim", "cmi_complete.rds"))



sg <- readRDS(here("Data", "Raw", "Bioclim", "sg_complete.rds"))

names(sg) <- 1960:2020
sg <- lapply(sg, as.data.frame)
sg <- lapply(sg, function(x) {rownames(x) <- bte4$TESSELLE; x})
sg <- lapply(sg, function(x) {
  colnames(x) <- c("start_gs", "end_gs", "length_gs",
  "pcp_1", "pcp_2","pcp_3", "pcp_4", "gdd_1", "gdd_2", "gdd_3", "gdd_4",
  "an_meanT", "an_minT", "an_maxT", "an_meanT_3", "an_range_T_3")
  x
  } )
sg <- do.call(rbind, lapply(sg, data.frame))
sg <- tibble::rownames_to_column(sg)


# Convert sg to a tidytable
sg_ttbl <- as_tidytable(sg)

# Separate the "rowname" column
sg_ttbl <- sg_ttbl %>%
  separate(rowname, into = c("year", "TESSELLE"), sep = "\\.", remove = FALSE)

# Convert back to a data frame
sg.test <- as.data.frame(sg_ttbl)

saveRDS(sg.test, here("Data", "Raw", "Bioclim", "sg_complete.rds"))






### Are all TESSELLE covered this way? No!

bte1_noP <- read.csv(here("Data", "BTE", "bte1_noP.csv"))[,-1]
bte2_noP <- read.csv(here("Data", "BTE", "bte2_noP.csv"))[,-1]
bte3_noP <- read.csv(here("Data", "BTE", "bte3_noP.csv"))[,-1]
bte4_noP <- read.csv(here("Data", "BTE", "bte4_noP.csv"))[,-1]
bte5_noP <- read.csv(here("Data", "BTE", "bte5_noP.csv"))[,-1]

# Are all TESSELLE from surveys in survey 4 (used to dowload the bioclim data)

extra1 <- bte1_noP$TESSELLE[! bte1_noP$TESSELLE %in% bte4$TESSELLE]
extra2 <- bte2_noP$TESSELLE[! bte2_noP$TESSELLE %in% bte4$TESSELLE]
extra3 <- bte3_noP$TESSELLE[! bte3_noP$TESSELLE %in% bte4$TESSELLE]
extra5 <- bte5_noP$TESSELLE[! bte5_noP$TESSELLE %in% bte4$TESSELLE]

## Combine: 

extra <- extra1
extra <- c(extra, extra2[! extra2 %in% extra])
extra <- c(extra, extra3[! extra3 %in% extra])
extra <- c(extra, extra5[! extra5 %in% extra])

###Add these
geom1 <- data_raw1 |>
  filter(TESSELLE %in% extra1)
rm(data_raw1)

geom2 <- data_raw2 |>
  filter(TESSELLE %in% extra2)
rm(data_raw2)

geom3 <- data_raw3 |>
  filter(TESSELLE %in% extra3)
rm(data_raw3)

geom5 <- data_raw5 |>
  filter(TESSELLE %in% extra5)
rm(data_raw5)

geom1 <- geom1 %>%
  select(intersect(colnames(geom1) , colnames(geom5)))
geom2 <- geom2 %>%
  select(intersect(colnames(geom2) , colnames(geom5)))
geom3 <- geom3 %>%
  select(intersect(colnames(geom3) , colnames(geom5)))
geom5 <- geom5 %>%
  select(intersect(colnames(geom5) , colnames(geom5)))

geom2 <- geom2 |>
  filter(! TESSELLE %in% geom1$TESSELLE)
geom3 <- geom3 |>
  filter(! TESSELLE %in% c(geom1$TESSELLE, geom2$TESSELLE))
geom5 <- geom2 |>
  filter(! TESSELLE %in% c(geom1$TESSELLE, geom2$TESSELLE, geom3$TESSELLE))

geom.comb <- bind_rows(geom1, geom2, geom3, geom5)
geom.comb <- st_as_sf(geom.comb)

geom.comb <- geom.comb %>% 
  st_transform("+proj=longlat +datum=WGS84 +no_defs")

#TODO fix error and run 
retrieveClimateData(years = 1960:2020, info = "sg", res = 60, geom = geom.comb, path = here("Data", "Raw", "Bioclim"))
retrieveClimateData(years = 1960:2018, info = "cmi", res = 60, geom = geom.comb, path = here("Data", "Raw", "Bioclim"))
