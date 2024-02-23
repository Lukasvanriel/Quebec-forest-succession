### Carte de l'aire d'étude ####

library(sf)
library(geodata)
library(RColorBrewer)
library(scales)
library(graphicsutils)

### DATA ####
st_layers("CLASSI_ECO_QC_GDB/CLASSI_ECO_QC.gdb")

ecoregion <- st_read("CLASSI_ECO_QC_GDB/CLASSI_ECO_QC.gdb",
                     layer = "N3_DOM_BIO")

ecoregion <- st_transform(ecoregion, 32188)

ecoregion_simple <- st_simplify(ecoregion, dTolerance = 500, preserveTopology = TRUE)

ecoregion_simple$DOM_BIO <- factor(ecoregion_simple$DOM_BIO, 
                                   levels = ecoregion_simple$DOM_BIO)

# Regions
ecoregion_cut <- ecoregion_simple[1:6,]
ecoregion_cut$DOM_BIO <- factor(ecoregion_cut$DOM_BIO, 
                                   levels = c(ecoregion_cut$DOM_BIO,7))

# study area:
reg <- st_read("CLASSI_ECO_QC_GDB/CLASSI_ECO_QC.gdb",
             layer = "N6_SREG_ECO")

reg <- reg[12,]
colnames(reg)
colnames(ecoregion_cut)

colnames(reg)[1] <- "NOM_DB"

reg <- reg[,colnames(ecoregion_cut)]

reg <- st_transform(reg, 32188)

regions <- rbind(ecoregion_cut, reg)

regions[7,"DOM_BIO"] <- 7

col_reg = rep <- c("darkgreen",rep("#be6f67", 2), "#fac679", rep("#be6f67", 3))


plot(regions["DOM_BIO"],  border = NA, lwd = .5, 
     pal =  alpha(rev(col_reg),.99), main=NA, key.pos=NULL)

legend(-76, 46.5, legend=c("Current", "Next: BTE", "Final: Province"), fill=c("darkgreen", "#fac679", "#be6f67"),
       density=c(1000, 1000, 1000), bty="n", text.font=1, cex=1.5) 

legend(750024, 5478140, legend=c("Current", "Next: Boreal-temperate ecotone", "Final: Province"), fill=c("darkgreen", "#fac679", "#be6f67"),
       density=c(1000, 1000, 1000), bty="n", text.font=1, cex=2.5,
       x.intersp = 0.1,
       y.intersp=c(0.7,0.7,0.7)) 

a$usr
points(750024, 5478140)



legend("top", legend=c("Current", "Next: Boreal-temperate ecotone", "Final: Province"), fill=c("darkgreen", "#fac679", "#be6f67"),
       density=c(1000, 1000, 1000), bty="n", text.font=1, cex=2.4,
       x.intersp = 0.1,
       y.intersp=c(0.7,0.7,0.7)) 


plot(c(1,100),c(1,100))

legend(5028315, 337070, legend=c("Current", "Next: Boreal-temperate ecotone", "Final: Province"), fill=c("darkgreen", "#fac679", "#be6f67"),
       density=c(1000, 1000, 1000), bty="n", text.font=1, cex=2.5,
       x.intersp = 0.1,
       y.intersp=c(0.7,0.7,0.7)) 



# Make legend



### MH --------
### GET NORTH AMERICA 

# get canada boundary map
can <- gadm("CAN", level = 1, resolution = 2, path = "/Users/mariehbrice/Documents/GitHub/classfication_ecologique_qc/carte_ERR/")
us <- gadm("US", level = 0, resolution = 2, path = "/Users/mariehbrice/Documents/GitHub/classfication_ecologique_qc/carte_ERR/")

# convert to sf
can_st <- st_as_sf(can)
us_st <- st_as_sf(us)

can_prj <- st_transform(can_st, 32188)
us_prj <- st_transform(us_st, 32188)

can_simple_prj <- st_simplify(can_prj, dTolerance = 800, preserveTopology = F)
us_simple_prj <- st_simplify(us_prj, dTolerance = 500, preserveTopology = F)

# Subset Quebec

qc <- can_simple_prj %>% subset(NAME_1 %in% c("Québec"))
qc_neigh <- can_simple_prj %>% subset(NAME_1 %in% c("Ontario", "New Brunswick", "Newfoundland and Labrador", "Nova Scotia", "Manitoba", "Nunavut"))

col_reg = c("#555753","#000000", "#8b3f96", "#6f739d", "#1f678b","#67b6bd",  "#a7a980", "#fac679", "#f9a469", "#be6f67")
col_reg = rep <- c("darkgreen",rep("#be6f67", 2), "#fac679", rep("#be6f67", 3))




xy <- data4bM[c("LONGI", "LATIT")]
xy <- st_as_sf(xy, 
               coords = c("LONGI", "LATIT"),
               crs = "+proj=longlat +datum=WGS84")


xy_proj <- st_transform(xy, crs(ecoregion_cut))
test <- st_transform(ecoregion_cut, "+proj=longlat +datum=WGS84")

polys = xy %>% 
  dplyr::summarise() %>%
  st_cast("POLYGON") %>% 
  st_convex_hull()

plot(polys)

test[7,] <- c(7, "study area", 1000, "A", "A", 100, 100, polys$geometry)

test$Study_area <- c("Final", "Final", "Final", "Next", "Final", "Final", "Current")

plot(ecoregion_cut["DOM_BIO"],  border = NA, lwd = .5, 
     pal =  alpha(rev(col_reg),.9), main="Study area")

plot(xy, col="#b2d922", pch = 19, axes = TRUE)

plot(test["DOM_BIO"],  border = NA, lwd = .5, 
     pal =  alpha(rev(col_reg),.99), main=NA, key.pos=NULL)


## Graticule
bb <- st_bbox(qc)
bb <- c(-700000, 4560000, 2000000,  7500000)
grat <- st_graticule(bb, crs = 32188, lon = seq(-100,-50, by = 5))
grat_x <- grat[grat$type == "E",]
grat_y <- grat[grat$type == "N",]

### Cartes de la région 

png("carte_region.png", width = 5.8, height = 5.5, res = 400, units = "in")

par(mar = c(1,1.7,.5,1), oma = c(0,0,0,.1), xpd = F)
plot(st_geometry(qc), border = "grey65", lwd = .9, col = "grey95",
     xlim = c(-380000, 1500000), ylim = c(4900000, 7000000))

box2(1:4, col = "grey55")

par(xpd = F)
plot(st_geometry(us_simple_prj), add = T, border = "grey65", lwd = .65, col = "grey75")
par(xpd = F)
plot(st_geometry(qc_neigh), add = T, border = "grey55", lwd = .65, col = "grey95")
par(xpd = F)
plot(st_geometry(qc), border = "grey35", lwd = 1, add = T)

par(xpd = F)
plot(st_geometry(grat), add = T, col = alpha('grey35', .3), lwd = .6) 

plot(ecoregion_simple["DOM_BIO"], add = T, border = "grey35", lwd = .5, pal =  alpha(rev(col_reg),.5))

axis(1, at = grat_x$x_start, labels = paste(abs(grat_x$degree),"°W"),
     cex.axis = .7, line = -1, tick = FALSE)
axis(2, at = grat_y$y_start, labels = paste(grat_y$degree,"°N"), 
     cex.axis = .7, las = 1, line = -.9, tick = FALSE)

dev.off()