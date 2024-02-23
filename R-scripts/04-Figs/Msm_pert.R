# Lukas Van Riel
# 2023-09-13
# 
rm(list = ls())

### Load packages ####
library(tidyverse)
library(msm)
library(minqa)
library(stringr)
library(here)
library(conflicted)
library(graphicsutils)
library(RColorBrewer)

conflicts_prefer(dplyr::filter)

### Data ####
model <- readRDS(here("Data-Output", "msm", "OLD", "msm4bM.CG.sc5e+05_pertlass+Tmean.rds"))

### Functions ####

plot_trans <- function(pmat, ci = NULL, cols = c("#FDF7F7", "red3", "#060000"), 
                       states_lab = NULL, labels = FALSE, main = NULL) {
  
  col <- colorRampPalette(cols)(200)
  
  if(!is.null(ci)) {
    ci <- matrix(paste0("(", round(ci[,,1],2), ", ", round(ci[,,2],2), ")"), 4)
  }
  
  if(is.null(states_lab)) states_lab = colnames(tr)
  
  n <- nrow(pmat)
  
  # Plot matrix
  image2(pmat, col = col, border = "white", lwd = 2)
  
  # Axis labels
  coordx <- seq(0, 1, len = n)
  coordy <- rev(coordx)
  axis(3, at = coordx, labels = states_lab, font = 2, las=2,
       tick = FALSE, cex.axis = 2.4, line = -1, col.axis = "grey15")
  axis(2, at = coordy, labels = states_lab, font = 2, 
       tick = FALSE, cex.axis = 2.4, las = 1, line = -.8, col.axis = "grey15")
  if(labels) {
    mtext("From", 2, font = 3, at = 1.0, las = 1, line = 1.5, cex = .8)
    mtext("To", 3, font = 3, at = -.05, line = 1.3, cex = .8)
  }
  #c(bottom, left, top, right)
  par(mar = c(5, 11, 11, 2) + 0.1)
  # Main
  mtext(main, 3, line = 1.7, font = 1, cex = .85)
  ?image2
  # Probabilities
  for(i in 1:n) {
    if(is.null(ci)) {
      text(x = coordx, y = coordy[i], labels = round(pmat[i,], 2), 
           col = ifelse(pmat[i,]<.3, "black", "white"), cex = 4, xpd = NA)
    } else {
      text(x = coordx, y = coordy[i]+.1, labels = round(pmat[i,], 2), 
           col = ifelse(pmat[i,]<.5, "black", "white"), cex = 4, xpd = NA)
      text(x = coordx, y = coordy[i]-.2, labels = ci[i,], 
           col = ifelse(pmat[i,]<.5, "black", "white"), cex = .8, xpd = NA)
    }
  }
}

plot_trans(p.base, labels=F, states_lab = 
             c('Pap. Birch', 'Other intol.', 'Y. Birch', 'Maple', 'Other Dec.', 'Bals. Fir', 'B/R Spruce', 'Jack Pine' , 'Other Con.'))

extract.colours <- function(c.df, p.matrix) {
  pmat <- round(p.matrix, 2)
  pmat <- sort(as.vector(pmat))
  print(pmat)
  
  matching_indices <- match(pmat, c.df$val)
  
  # Step 3: Retrieve the corresponding colors from the "colours" column
  matched_colors <- c.df$colours[matching_indices]
}


colours.df <- data.frame(val=(seq(1,33)-17)/100, colours=colorRampPalette(c("blue", "white", "red3"))(33))

plot_trans_diff <- function(pmat, ci = NULL, cols = c("lightblue", "white", "red3"), 
                            states_lab = NULL, labels = FALSE, main = NULL) {
  col <- colorRampPalette(cols)(100)
  
  #col <- cols

  if(!is.null(ci)) {
    ci <- matrix(paste0("(", round(ci[,,1],2), ", ", round(ci[,,2],2), ")"), 4)
  }
  
  if(is.null(states_lab)) states_lab = colnames(tr)
  
  n <- nrow(pmat)
  
  # Plot matrix
  image2(pmat, col  = col, border = "white", lwd = 2)

  # Axis labels
  coordx <- seq(0, 1, len = n)
  coordy <- rev(coordx)
  axis(3, at = coordx, labels = states_lab, font = 2, las=2,
       tick = FALSE, cex.axis = 2.4, line = -1, col.axis = "grey15")
  axis(2, at = coordy, labels = states_lab, font = 2, 
       tick = FALSE, cex.axis = 2.4, las = 1, line = -.8, col.axis = "grey15")
  if(labels) {
    mtext("From", 2, font = 3, at = 1.0, las = 1, line = 1.5, cex = 2.8)
    mtext("To", 3, font = 3, at = -.05, line = 1.3, cex = 2.8)
  }
  
  # Main
  mtext(main, 3, line = 1.7, font = 1, cex = .85)
  
  # Probabilities
  for(i in 1:n) {
    if(is.null(ci)) {
      text(x = coordx, y = coordy[i], labels = round(pmat[i,], 2), 
           col = ifelse(pmat[i,]> -.09, "black", "white"), cex = 4, xpd = NA)
    } else {
      text(x = coordx, y = coordy[i]+.1, labels = round(pmat[i,], 2), 
           col = ifelse(pmat[i,]> -.09, "black", "white"), cex = 4, xpd = NA)
      text(x = coordx, y = coordy[i]-.2, labels = ci[i,], 
           col = ifelse(pmat[i,]> -.09, "black", "white"), cex = .8, xpd = NA)
    }
  }
}
sort(as.vector(round(p.fire-p.base,2)))

plot_trans_diff(p.fire-p.base, labels=F, cols=extract.colours(colours.df, p.fire-p.base), states_lab =
                  c('Pap. Birch', 'Other intol.', 'Y. Birch', 'Maple', 'Other Dec.', 'Bals. Fir', 'B/R Spruce', 'Jack Pine' , 'Other Con.'))

plot_trans_diff(p.fire-p.base, labels=F, states_lab =
                  c('Pap. Birch', 'Other intol.', 'Y. Birch', 'Maple', 'Other Dec.', 'Bals. Fir', 'B/R Spruce', 'Jack Pine' , 'Other Con.'))

plot_trans_diff(p.harv-p.base, labels=F, col=extract.colours(colours.df, p.harv-p.base),
                states_lab = c('Pap. Birch', 'Other intol.', 'Y. Birch', 'Maple', 'Other Dec.', 'Bals. Fir', 'B/R Spruce', 'Jack Pine' , 'Other Con.'))

plot_trans_diff(p.outb-p.base, labels=F, col=extract.colours(colours.df, p.outb-p.base),
                states_lab = c('Pap. Birch', 'Other intol.', 'Y. Birch', 'Maple', 'Other Dec.', 'Bals. Fir', 'B/R Spruce', 'Jack Pine' , 'Other Con.'))

t <- round(p.fire-p.base,2)
plot.matrix(t)


plot_pmatrix <- function(mod, t = 1:40, covar = "mean", ci = "none", 
                         st_col = c("#158282", "#A1BD93","#FEAC19", "#D43650"),
                         states = c("Boreal", "Mixed", "Pioneer", "Temperate"), 
                         main = F, yaxis = T) {
  
  
  
  if(is.list(covar)) {
    p <- list()
    for(i in 1:length(covar)) {
      p_tmp <- pmatrix.msm(mod, t = t, covariates = covar[[i]], ci = ci)
      class(p_tmp) <- "array"
      p[[i]] <- p_tmp
    }
  } else {
    p <- pmatrix.msm(mod, t = t, covariates = covar, ci = ci)
    class(p) <- "array"
  }
  
  
  for(st_from in 1:4) {
    plot0(x = t, xlim = c(0, max(t)), ylim = c(0,1), xaxs = "i", yaxs = "i")
    
    axis(1, labels = F, tcl= -0.5, col = "grey35")
    if(st_from==4) axis(1, tick = F, cex = .9)
    
    axis(2, labels = F, tcl= -0.5, col = "grey35")
    if(yaxis) axis(2, las = 1, tick = F, cex = .9)
    #box2()
    
    if(main) mtext(paste("From", states[st_from]), 3)
    
    if(is.list(p)) {
      for(i in 1:length(p)){
        from <- p[[i]][st_from, , ]
        
        for(st_to in 1:4) {
          lty = 1:3
          lines(t, from[st_to,], col = st_col[st_to], lwd = 1.3, lty = lty[i])
        }  
      }
    } else {
      from <- p[st_from, , ]
      
      for(st_to in 1:4) {
        lines(t, from[st_to,], col = st_col[st_to], lwd = 1.3, lty = 1)
      }  
    }
    
  }
  
}

create_trans_timelapse <- function(from, to, t=50, p=c("p.b", "p.f", "p.h", "p.o")){
  y <- vector(mode = "list", length = length(p))
  for(j in 1:length(y)){
    y_temp <- vector(mode = "double", length = t)
    for(k in 1:length(y_temp)) {
      y_temp[k] <- get(p[j])[[k]][from, to]
    }
    y[[j]] <- y_temp
  }
  
  jpeg(file=here("R-scripts","04-Figs", paste0("base-pert-timelapse-", as.character(from), "-",
                                               as.character(to),".jpg"))) 
  
  plot(1:t,y[[1]], ylim=c(0,max(y[[1]], y[[2]], y[[3]], y[[4]])), pch=19, 
       xlab="Time (yr)", ylab="P")
  
  points(1:t,y[[2]], col="red", pch=19)
  points(1:t,y[[3]], col="darkgreen", pch=19)
  points(1:t,y[[4]], col="purple", pch=19)
  
  dev.off()
}

plot_trans_timelapse <- function(from, to, t=50, p=c("p.b", "p.f", "p.h", "p.o"),
                                 trans.title){
  y <- vector(mode = "list", length = length(p))
  for(j in 1:length(y)){
    y_temp <- vector(mode = "double", length = t)
    for(k in 1:length(y_temp)) {
      y_temp[k] <- get(p[j])[[k]][from, to]
    }
    y[[j]] <- y_temp
  }
  plot(1:t,y[[1]], ylim=c(0,0.27), pch=19, 
       xlab="Time (yr)", ylab="P", cex=3, main=trans.title,
       cex.axis=2.5, cex.main=3)
  
  points(1:t,y[[2]], col="red", pch=19, cex=3)
  points(1:t,y[[3]], col="darkgreen", pch=19, cex=3)
  points(1:t,y[[4]], col="purple", pch=19, cex=3)
}

#### Analysis ------

#### TEST stuff -------
data.bte <- read.csv(here("Data", "BTE", "bte_msm_ready.csv")) %>% 
  select(-X)

data.4bM <- read.csv(here("Data", "BTE", paste0("bte_", "4bM", "_msm_ready.csv")))
data.4bS <- read.csv(here("Data", "BTE", paste0("bte_", "4bS", "_msm_ready.csv")))
data.4bT <- read.csv(here("Data", "BTE", paste0("bte_", "4bT", "_msm_ready.csv")))

m.4bM.pcl <- readRDS(here("Data-Output", "msm", "4bM", "msm_4bM.CG.sc3e+05_pertclass.rds"))
m.4bM.T <- readRDS(here("Data-Output", "msm", "4bM", "msm_4bM.CG.sc3e+05_Tmean.rds"))
m.4bM.s <- readRDS(here("Data-Output", "msm", "4bM", "msm_4bM.CG.sc3e+05_Soil.rds"))
m.4bM.combo <- readRDS(here("Data-Output", "msm", "OLD", "msm4bM.CG.sc5e+05_pertlass+Tmean.rds"))
m.4bM.T.no1 <- readRDS("~/Desktop/msm_4bM_NO1.CG.sc3e+05_Tmean.rds")
m.4bM.T.no1_conv <- readRDS("~/Desktop/msm_4bM_NO1.CG.sc3e+05_Tmean_e-16.rds")

round(m.4bM.T.no1_conv$Qmatrices$baseline,4)
round(m.4bM.T.no1$Qmatrices$baseline,4)
round(m.4bM.T$Qmatrices$baseline,4)
round(pmatrix.msm(m.4bM.T.no1, t=10),4)
round(pmatrix.msm(m.4bM.T, t=10),4)

m.4bS.s <- readRDS(here("Data-Output", "msm", "4bS", "msm_4bS.CG.sc3e+05_Soil.rds"))
m.4bT.s <- readRDS(here("Data-Output", "msm", "4bT", "msm_4bT.CG.sc3e+05_Soil.rds"))

hist(data.4bT$cov_soil)

## Which zones to compare for proposal?
plot(data.bte$LONGI, data.bte$LATIT, pch=19)
points(data.4bM$LONGI, data.4bM$LATIT, pch=19, col="red")
points(data.4bS$LONGI, data.4bS$LATIT, pch=19, col="blue")
points(data.4bT$LONGI, data.4bT$LATIT, pch=19, col="green")

## Check out models
round(m.4bM.pcl$Qmatrices$baseline,3)
round(m.4bM.T$Qmatrices$baseline,3)
round(m.4bM.s$Qmatrices$baseline,3)
round(m.4bM.combo$Qmatrices$baseline,3)

round(pmatrix.msm(m.4bM.pcl, t=10),3)
round(pmatrix.msm(m.4bM.T, t=10),3)
round(pmatrix.msm(m.4bM.s, t=10),3)

qmatrix.msm(m.4bM.pcl)
sojourn.msm(m.4bM.pcl)
pnext.msm(m.4bM.pcl)

# 
round(m.4bM.T$Qmatrices$baseline,3)
round(m.4bM.T$Qmatrices$logbaseline,3) #just the log of the baseline matrix
round(m.4bM.T$Qmatrices$cov_Tmean,3)

round(pmatrix.msm(m.4bM.T, t=10, covariates = list(cov_Tmean=0)),4)
round(pmatrix.msm(m.4bM.T, t=10, covariates = list(cov_Tmean=2)),4)
round(pmatrix.msm(m.4bM.T, t=10, covariates = list(cov_Tmean=-2)),4)

round(pmatrix.msm(m.4bM.T, t=10, covariates = list(cov_Tmean=0)),4)

#Soil for three zones:

round(m.4bM.s$Qmatrices$cov_soil,3)
round(m.4bS.s$Qmatrices$cov_soil,3)
round(m.4bT.s$Qmatrices$cov_soil,3)

round(pmatrix.msm(m.4bM.s, t=10),4)
round(pmatrix.msm(m.4bS.s, t=10),4)
round(pmatrix.msm(m.4bT.s, t=10),4)

m.4bT.s$opt$hessian

list.files("~/Desktop/4bM", pattern = "\\.rds$")
list.files("~/Desktop/4bM")
folders <- Sys.glob("/Users/lukas/Desktop/RProjects/Quebec-forest-succession/Data-Output/msm/4*")
files <- Sys.glob("~/Desktop/4bT/msm*.rds")



for(fo in folders) {
  files <- Sys.glob(Sys.glob(paste0(fo,"/", "msm*.rds")))
  for(f in files) {
    print(f)
    print(readRDS(f))
    print("--------------------------------------------")
  }
}


a <- readRDS("~/Desktop/4bT/msm_4bT.CG.sc750000_soil.rds")
b <- readRDS("~/Desktop/4bT/msm_4bT.CG.sc3e+05_soil.rds")

round(a$Qmatrices$baseline, 4)
round(b$Qmatrices$baseline, 4)

round(pmatrix.msm(a, t=10), 4)
round(pmatrix.msm(b, t=10), 4)

##POSTER##----

#summary(model)
#Used: Tmean + perturb.class: 0=no, 1=burn, 2=cut, 3=outbreak

round(model$Qmatrices$baseline,3)
round(model$Qmatrices$cov_pert_class1,3)
round(model$Qmatrices$baseline,3) + round(model$Qmatrices$cov_pert_class1,3)



q.baseline <- round(qmatrix.msm(model)$estimates, 5)

q.fire <- round(qmatrix.msm(model, 
                                covariates=list(cov_pert_class=1))$estimates, 5)
q.harv <- round(qmatrix.msm(model, 
                                   covariates=list(cov_pert_class=2))$estimates, 5)
q.outb <- round(qmatrix.msm(model, 
                                   covariates=list(cov_pert_class=3))$estimates, 5)


df <- data.frame(cov=factor(levels = c("base", "fire", "harv", "outb")), 
                 trans=numeric(), q=numeric())
for(i in 1:length(q.baseline)) {
  b <- which(q.baseline==q.baseline[i], arr.ind=TRUE)
  df[nrow(df)+1,] <- c("base", b[1]*10+b[2], q=q.baseline[i])
  
  f <- which(q.fire==q.fire[i], arr.ind=TRUE)
  df[nrow(df)+1,] <- c("fire", f[1]*10+f[2], q=q.fire[i])
  
  h <- which(q.harv==q.harv[i], arr.ind=TRUE)
  df[nrow(df)+1,] <- c("harv", h[1]*10+h[2], q=q.harv[i])
  
  o <- which(q.outb==q.outb[i], arr.ind=TRUE)
  df[nrow(df)+1,] <- c("outb", o[1]*10+o[2], q=q.outb[i])
}

ggplot(df) +
  geom_point(aes(x=trans, y=q, col=cov)) +
  labs(x = "Transition",
       y = "Estimated q") +
  theme_minimal()
  
p.av <- round(pmatrix.msm(model, t=10), 5)

p.base <- round(pmatrix.msm(model, t=10, covariates = list(cov_pert_class=0)), 5)
p.fire <- round(pmatrix.msm(model, t=10, covariates = list(cov_pert_class=1)), 5)
p.harv <- round(pmatrix.msm(model, t=10, covariates = list(cov_pert_class=2)), 5)
p.outb <- round(pmatrix.msm(model, t=10, covariates = list(cov_pert_class=3)), 5)

write.csv(p.base, "~/Desktop/pbase_10yr.csv", row.names = F)
write.csv(p.fire, "~/Desktop/pfire_10yr.csv", row.names = F)
write.csv(p.harv, "~/Desktop/pharv_10yr.csv", row.names = F)
write.csv(p.outb, "~/Desktop/poutb_10yr.csv", row.names = F)

write.csv(round(p.fire - p.base, 5), "~/Desktop/delta_pfire_10yr.csv", row.names = F)
write.csv(round(p.harv - p.base, 5), "~/Desktop/delta_pharv_10yr.csv", row.names = F)
write.csv(round(p.outb - p.base, 5), "~/Desktop/delta_poutb_10yr.csv", row.names = F)

p.av <- round(pmatrix.msm(model, t=10), 5)
p.base <- round(pmatrix.msm(model, t=10, covariates = list(cov_pert_class=0)), 5)

p.Tbase <- round(pmatrix.msm(model, t=10, covariates = list(cov_Tmean=0)), 5)

p.T1 <- round(pmatrix.msm(model, t=10, covariates = list(cov_Tmean=1)), 5)
p.Tmin1 <- round(pmatrix.msm(model, t=10, covariates = list(cov_Tmean=-1)), 5)

p.T3 <- round(pmatrix.msm(model, t=10, covariates = list(cov_Tmean=3)), 5)
p.Tmin3 <- round(pmatrix.msm(model, t=10, covariates = list(cov_Tmean=-3)), 5)

write.csv(p.Tbase, "~/Desktop/pTbase_10yr.csv", row.names = F)

write.csv(p.T1, "~/Desktop/pT1_10yr.csv", row.names = F)
write.csv(p.Tmin1, "~/Desktop/pT-1_10yr.csv", row.names = F)

write.csv(p.T3, "~/Desktop/pT3_10yr.csv", row.names = F)
write.csv(p.Tmin3, "~/Desktop/pT-3_10yr.csv", row.names = F)

write.csv(round(p.T1 - p.Tbase, 5), "~/Desktop/delta_T1_10yr.csv", row.names = F)
write.csv(round(p.Tmin1 - p.Tbase, 5), "~/Desktop/delta_T-1_10yr.csv", row.names = F)

write.csv(round(p.T3 - p.Tbase, 5), "~/Desktop/delta_T3_10yr.csv", row.names = F)
write.csv(round(p.Tmin3 - p.Tbase, 5), "~/Desktop/delta_T-3_10yr.csv", row.names = F)


heatmap(p.base)
heatmap(p.harv)

heatmap(p.base)
heatmap(p.base-p.fire)
heatmap(p.base-p.harv)
heatmap(a, Colv=NA, Rowv=NA)
a <- p.base-p.fire

par(mfrow=c(1,2))
heatmap(p.base, Colv=NA, Rowv=NA)
heatmap(p.fire, Colv=NA, Rowv=NA)
heatmap(p.harv, Colv=NA, Rowv=NA)
heatmap(p.outb, Colv=NA, Rowv=NA)

#1: Paper birch
#2: Other shade intolerant species
#3: Yellow birch
#4: Sugar and red maple
#5: Other deciduous
#6: Balsam fir
#7: Red and black spruce
#8: Jack pine
#9: Other coniferous


######PLOT MATRICES#######
if(F) {

# col <- colorRampPalette(c("lightblue", "red3"))(99)[ii_l]
# 
# ii_l <- cut(seq(-0.17, 0, by=0.01), breaks = seq(-0.17, 0, len = 100), 
#           include.lowest = TRUE)
# col_l <- colorRampPalette(c("lightblue", "white"))(99)[ii_l]
# show_col(col_l)
# ii_h <- cut(seq(0, 0.17, by=0.01), breaks = seq(0, 0.17, len = 100), 
#             include.lowest = TRUE)
# col_h <- colorRampPalette(c("white", "red3"))(99)[ii_h]
# show_col(col_h)
# 
# col <- c(col_l, "white", col_h)
# show_col(col)
# 
# colorRampPalette(c("lightblue", "white", "red3"))(35)
# show_col(colorRampPalette(c("lightblue", "white", "red3"))(35))
  

plot_trans_diff(p.base-p.fire, labels=F, col=extract.colours(colours.df, p.base-p.fire),
                states_lab = c('1', '2', '3', '4', '5', '6', '7', '8' , '9'))

#average
plot_trans(p.av, labels=T, states_lab = 
             c('1', '2', '3', '4', '5', '6', '7', '8' , '9'))

#Baseline
plot_trans(p.base, labels=F, states_lab = 
             c('1', '2', '3', '4', '5', '6', '7', '8' , '9'))

par(mar = c(5.1, 15.1, 15.1, 10.1))
plot_trans(p.base, labels=F, states_lab = 
             c('Pap. Birch', 'Other intol.', 'Y. Birch', 'Maple', 'Other Dec.', 'Bals. Fir', 'B/R Spruce', 'Jack Pine' , 'Other Con.'))

#With disturbances
plot_trans(p.fire, labels=T, states_lab = 
             c('1', '2', '3', '4', '5', '6', '7', '8' , '9'))

plot_trans(p.harv, labels=T, states_lab = 
             c('1', '2', '3', '4', '5', '6', '7', '8' , '9'))

plot_trans(p.outb, labels=T, states_lab = 
             c('1', '2', '3', '4', '5', '6', '7', '8' , '9'))

# Differences with disturbances
par(mar = c(5.1, 4.1, 4.1, 10.1))
dev.off()

show_col(colorRampPalette(c("lightblue", "white", "red3"))(33))

max(c(p.base-p.fire, p.base-p.harv, p.base-p.outb))
min(c(p.base-p.fire, p.base-p.harv, p.base-p.outb))

colours.df <- data.frame(val=seq(-0.16, 0.16, by=0.01), colours=colorRampPalette(c("blue", "white", "red3"))(33))



extract.colours <- function(c.df, p.matrix) {
  c.df %>% 
    filter(val<(max(p.matrix)+0.01)) %>% 
    filter(val>(min(p.matrix)-0.01)) %>% 
    pull(colours)
}

par(mar = c(5.1, 15.1, 15.1, 10.1))

plot_trans_diff(p.base-p.fire, labels=F, col=extract.colours(colours.df, p.base-p.fire), states_lab =
                  c('Pap. Birch', 'Other intol.', 'Y. Birch', 'Maple', 'Other Dec.', 'Bals. Fir', 'B/R Spruce', 'Jack Pine' , 'Other Con.'))

plot_trans_diff(p.base-p.harv, labels=F, col=extract.colours(colours.df, p.base-p.harv),
                states_lab = c('Pap. Birch', 'Other intol.', 'Y. Birch', 'Maple', 'Other Dec.', 'Bals. Fir', 'B/R Spruce', 'Jack Pine' , 'Other Con.'))

plot_trans_diff(p.base-p.outb, labels=F, col=extract.colours(colours.df, p.base-p.outb),
                states_lab = c('Pap. Birch', 'Other intol.', 'Y. Birch', 'Maple', 'Other Dec.', 'Bals. Fir', 'B/R Spruce', 'Jack Pine' , 'Other Con.'))

#c('1', '2', '3', '4', '5', '6', '7', '8' , '9')


##
par(mar = c(1, 18, 8, 1))
#a <- matrix(c(0.75, 0.2 , 0.05 ,0.2,0.7,0.1,0.2,0.2,0.6), byrow = T , ncol=3, nrow=3)
plot_trans(a, labels=T, states_lab = 
             c("State 1", "State 2", "State 3"))

dev.off()

}


######CREATE TIME PROBABILITIES#######

pmatrix.msm(model, t = 1, covariates = "mean")

m <- matrix(rep(0,81), nrow=9, ncol=9)

p.b <- vector(mode = "list", length = 50); names(p.b) <- 1:50
p.f <- vector(mode = "list", length = 50); names(p.f) <- 1:50
p.h <- vector(mode = "list", length = 50); names(p.h) <- 1:50
p.o <- vector(mode = "list", length = 50); names(p.o) <- 1:50
for(i in 1:50) {
  p.b[[i]] <- m
  p.f[[i]] <- m
  p.h[[i]] <- m
  p.o[[i]] <- m
}

for(i in 1:length(p.b)) {
  p.b[[i]] <- round(pmatrix.msm(model, t=i, covariates = list(cov_pert_class=0)), 5)
  p.f[[i]] <- round(pmatrix.msm(model, t=i, covariates = list(cov_pert_class=1)), 5)
  p.h[[i]] <- round(pmatrix.msm(model, t=i, covariates = list(cov_pert_class=2)), 5)
  p.o[[i]] <- round(pmatrix.msm(model, t=i, covariates = list(cov_pert_class=3)), 5)
}


##### OWN try
p <- c("p.b", "p.f", "p.h", "p.o")
t <- 1:50



create_trans_timelapse(from=3, to=8)

for(from in 1:9) {
  for(to in 1:9) {
    create_trans_timelapse(from=from, to=to)
  }
}

plot(t,y[[1]], ylim=c(0,max(y[[1]], y[[2]], y[[3]], y[[4]])), pch=19, 
     xlab="Time (yr)", ylab="P", ax.cex=2)

?plot

points(t,y[[2]], col="red", pch=19)
points(t,y[[3]], col="darkgreen", pch=19)
points(t,y[[4]], col="purple", pch=19)


##
t <- 1:40

##
plot0(x = t, xlim = c(0, max(t)), ylim = c(0,1), xaxs = "i", yaxs = "i")

axis(1, labels = F, tcl= -0.5, col = "grey35")
if(st_from==4) axis(1, tick = F, cex = .9)

axis(2, labels = F, tcl= -0.5, col = "grey35")
if(yaxis) axis(2, las = 1, tick = F, cex = .9)
#box2()

if(main) mtext(paste("From", states[st_from]), 3)

if(is.list(p)) {
  for(i in 1:length(p)){
    from <- p[[i]][st_from, , ]
    
    for(st_to in 1:4) {
      lty = 1:3
      lines(t, from[st_to,], col = st_col[st_to], lwd = 1.3, lty = lty[i])
    }  
  }
} else {
  from <- p[st_from, , ]
  
  for(st_to in 1:4) {
    lines(t, from[st_to,], col = st_col[st_to], lwd = 1.3, lty = 1)
  }  
}

##
if(T){
dev.off()
o <- matrix(c(1,2,3,4,5,6,7,8,8),nrow = 3,ncol = 3,byrow = TRUE)

layout(mat = o, heights = c(0.4, 0.4, 0.2))

f <- c(1, 3, 4, 4, 5, 7)
t <- c(7, 4, 1, 6, 4, 9)

for(i in 1:6) {
  plot_trans_timelapse(from=f[i], to=t[i])
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("None", "Fire", "Harvest", "Outbreak"), col = c("black","red", "darkgreen", "purple"), lwd = 9, xpd = TRUE, horiz = TRUE, cex = 2, seg.len=1.5, bty = 'n')
# xpd = TRUE makes the legend plot to the figure
}




p.b <- round(pmatrix.msm(model, t=10, covariates = list(cov_pert_class=0)), 5)
p.f <- round(pmatrix.msm(model, t=10, covariates = list(cov_pert_class=1)), 5)
p.h <- round(pmatrix.msm(model, t=10, covariates = list(cov_pert_class=2)), 5)
p.o <- round(pmatrix.msm(model, t=10, covariates = list(cov_pert_class=3)), 5)

###
tit <- c("Pap. birch to B/R spruce", "Y. birch to Maple", "Maple to Pap. birch",
         "Maple to Bals. fir", "Other Dec. to Maple", "B/R spruce to Other Con.")

op <- par(mfrow = c(2,3),
          oma = c(10,7,0,0) + 0.1,
          mar = c(1,3,6,1) + 0.1)

for(i in 1:6) {
  plot_trans_timelapse(from=f[i], to=t[i], trans.title = tit[i])
}

title(xlab = "Time (yr)",
      ylab = "Transition probability",
      outer = TRUE, cex.lab=5, line = 3)

op <- par(fig = c(0, 1, 0, 1), oma = c(0, 5, 0, 0), mar = c(0, 0, 0, 0), new = TRUE, cex = 1)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("None", "Fire", "Harvest", "Outbreak"), col = c("black","red", "darkgreen", "purple"),
       lwd = 19, xpd = TRUE, horiz = TRUE, cex = 1, seg.len=0.0005, bty = 'n')

dev.off()

create_trans_timelapse(1,2, t=50)

#### Good version I think?

tit <- c("Pap. birch to B/R spruce", "Y. birch to Maple", "Maple to Pap. birch",
         "Maple to Bals. fir", "Other Dec. to Maple", "B/R spruce to Other Con.",
         "Bals. fir to Y. Birch", "B/R spruce to Maple")

op <- par(mfrow = c(2,4),
          oma = c(10,7,0,0) + 0.1,
          mar = c(1,3,6,1) + 0.1)


f <- c(1, 3, 4, 4, 5, 7, 6, 7)
t <- c(7, 4, 1, 6, 4, 9, 3, 4)

for(i in 1:8) {
  plot_trans_timelapse(from=f[i], to=t[i], trans.title = tit[i])
}

title(xlab = "Time (yr)",
      ylab = "Transition probability",
      outer = TRUE, cex.lab=5, line = 3)

op <- par(fig = c(0, 1, 0, 1), oma = c(0, 5, 0, 0), mar = c(0, 0, 0, 0), new = TRUE, cex = 1)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend('bottom',legend = c("None", "Fire", "Harvest", "Outbreak"), col = c("black","red", "darkgreen", "purple"),
       lwd = 19, xpd = TRUE, horiz = TRUE, cex = 1, seg.len=0.0005, bty = 'n')

dev.off()
