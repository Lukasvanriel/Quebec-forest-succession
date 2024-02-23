#### Load packages ####
library(tidyverse)
library(INLA)
library(INLAjoint)
library(here)

### Data ###

inla.all <- readRDS(here("Data-Output", "INLA", "3j.T.C.S.rds"))

inla.summ <- readRDS(here("Data-Output", "INLA", "3j.T.C.S", "3j.T.C.S_summary.fixed.RDS"))

inla.hype <- readRDS(here("Data-Output", "INLA", "3j.T.C.S", "3j.T.C.S_summary.hyperpar.RDS"))

inla.mlik <- readRDS(here("Data-Output", "INLA", "3j.T.C.S", "3j.T.C.S_mlik.RDS"))

inla.dic <- readRDS(here("Data-Output", "INLA", "3j.T.C.S", "3j.T.C.S_dic.RDS"))

inla.waic <- readRDS(here("Data-Output", "INLA", "3j.T.C.S", "3j.T.C.S_waic.RDS"))

inla.cpuu <- readRDS(here("Data-Output", "INLA", "3j.T.C.S", "3j.T.C.S_cpu.used.RDS"))


### Functions ###
#Create own functions to emulate what existing INLA functions can do

inl.summary <- function(sum.fix, sum.hyp, mlik, dic, waic, cpu.time) {
  N.cov <- (nrow(sum.fix) / nrow(sum.hyp)) - 1 #; print(N.cov)
  for(i in 1:nrow(sum.hyp)) {
    out <- rbind(round(sum.hyp[i, 1:5],4),
                 rep(0.1234, 5),
                 round(sum.fix[((i-1)*(N.cov + 1) + 2):(i*(N.cov + 1)), 1:5],4))
    rownames(out)[1:2] <- c(paste0("Weibull (shape)_S", as.character(i)), paste0("Weibull (scale)_S", as.character(i)))
    
    cat(noquote(paste0("Survival outcome (S", as.character(i), ")")), fill = T)
    print(out)
    cat(noquote(" "), fill = T)
  }
  
  colnames(mlik) <- ""
  print(mlik, fill=T)
  cat("", fill = T)
  
  cat(noquote(paste0("Deviance Information Criterion: ", as.character(round(dic$dic,1)))), fill = T)
  cat(noquote(paste0("Widely applicable Bayesian information criterion: ", as.character(round(waic$waic,1)))), fill = T)
  cat(noquote(paste0("Computation time: ", as.character(round(cpu.time[4], 2)))), fill = T)
}
inl.summary(inla.summ, inla.hype, inla.mlik, inla.dic, inla.waic, inla.cpuu)

set.back.inla <- function(string){
  object <- list()
  
  object <- c(object, list(summary.fixed = readRDS(here("Data-Output", "INLA", string, paste0(string, "_summary.fixed.RDS")))))
  object <- c(object, list(marginals.fixed = readRDS(here("Data-Output", "INLA", string, paste0(string, "_marginals.fixed.RDS")))))
  object <- c(object, list(summary.random = readRDS(here("Data-Output", "INLA", string, paste0(string, "_summary.random.RDS")))))
  object <- c(object, list(marginals.random = readRDS(here("Data-Output", "INLA", string, paste0(string, "_marginals.random.RDS")))))
  object <- c(object, list(summary.hyperpar = readRDS(here("Data-Output", "INLA", string, paste0(string, "_summary.hyperpar.RDS")))))
  object <- c(object, list(marginals.hyperpar = readRDS(here("Data-Output", "INLA", string, paste0(string, "_marginals.hyperpar.RDS")))))
  object <- c(object, list(summary.linear.predictor = readRDS(here("Data-Output", "INLA", string, paste0(string, "_summary.linear.predictor.RDS")))))
  object <- c(object, list(marginals.linear.predictor = readRDS(here("Data-Output", "INLA", string, paste0(string, "_marginals.linear.predictor.RDS")))))
  object <- c(object, list(summary.fitted.values = readRDS(here("Data-Output", "INLA", string, paste0(string, "_summary.fitted.values.RDS")))))
  object <- c(object, list(marginals.fitted.values = readRDS(here("Data-Output", "INLA", string, paste0(string, "_marginals.fitted.values.RDS")))))
  
  if(file.exists(here("Data-Output", "INLA", string, paste0(string, "_summary.lincomb.RDS")))) {
    object <- c(object, list(summary.lincomb = readRDS(here("Data-Output", "INLA", string, paste0(string, "_summary.lincomb.RDS")))))
  } else {
    object <- c(object, list(summary.lincomb = NULL))
  }
  #object <- c(object, list(summary.lincomb = readRDS(here("Data-Output", "INLA", string, paste0(string, "_summary.lincomb.RDS")))))
  
  object <- c(object, list(marginals.lincomb = readRDS(here("Data-Output", "INLA", string, paste0(string, "_marginals.lincomb.RDS")))))
  object <- c(object, list(selection = readRDS(here("Data-Output", "INLA", string, paste0(string, "_selection.RDS")))))
  object <- c(object, list(dic = readRDS(here("Data-Output", "INLA", string, paste0(string, "_dic.RDS")))))
  object <- c(object, list(cpo = readRDS(here("Data-Output", "INLA", string, paste0(string, "_cpo.RDS")))))
  object <- c(object, list(po = readRDS(here("Data-Output", "INLA", string, paste0(string, "_po.RDS")))))
  object <- c(object, list(residuals = readRDS(here("Data-Output", "INLA", string, paste0(string, "_residuals.RDS")))))
  object <- c(object, list(waic = readRDS(here("Data-Output", "INLA", string, paste0(string, "_waic.RDS")))))
  object <- c(object, list(mlik = readRDS(here("Data-Output", "INLA", string, paste0(string, "_mlik.RDS")))))
  
  if(file.exists(here("Data-Output", "INLA", string, paste0(string, "_neffp.RDS")))) {
    object <- c(object, list(neffp = readRDS(here("Data-Output", "INLA", string, paste0(string, "_neffp.RDS")))))
  } else {
    object <- c(object, list(neffp = NULL))
  }
  #object <- c(object, list(neffp = readRDS(here("Data-Output", "INLA", string, paste0(string, "_neffp.RDS")))))
  
  object <- c(object, list(mode = readRDS(here("Data-Output", "INLA", string, paste0(string, "_mode.RDS")))))
  object <- c(object, list(call = readRDS(here("Data-Output", "INLA", string, paste0(string, "_call.RDS")))))
  
  if(file.exists(here("Data-Output", "INLA", string, paste0(string, "_formula.RDS")))) {
    object <- c(object, list(formula = readRDS(here("Data-Output", "INLA", string, paste0(string, "_formula.RDS")))))
  } else {
    object <- c(object, list(neffp = NULL))
  }
  #object <- c(object, list(formula = readRDS(here("Data-Output", "INLA", string, paste0(string, "_formula.RDS")))))
  
  object <- c(object, list(nhyper = readRDS(here("Data-Output", "INLA", string, paste0(string, "_nhyper.RDS")))))
  object <- c(object, list(cpu.used = readRDS(here("Data-Output", "INLA", string, paste0(string, "_cpu.used.RDS")))))
  
  class(object) <- c("INLAjoint", "inla")
  return(object)
}

set.back.inla2 <- function(string){
  object <- list()
  
  file.list <- list.files(path=here("Data-Output", "INLA", string), pattern=".RDS", all.files=TRUE, 
                          full.names=F)
  
  for(file in file.list) {
    file <- substr(file, 10,nchar(file))
    file <- substr(file, 1, nchar(file) - 4)
    
    object[[file]] <- readRDS(here("Data-Output", "INLA", string, paste0(string, "_", file, ".RDS")))
    print(names(object)[length(object)])
  }
 
  class(object) <- c("INLAjoint", "inla")
  return(object)
}

test.f <- set.back.inla2("allj.T.C")


test <- set.back.inla("3j.T.C.S")






### Documentation ###

?joint()
inla.doc("weibullsurv")

vignette("INLAjoint")

### ### 

inlA <- set.back.inla2("allj.T.C")

a <- summary(inlA)

inlA$summary.fixed
inlA$summary.hyperpar

plot.cov <- plot(inlA)$Outcomes
plot.hyp <- plot(inlA)$Baseline

plot.cov$S22




plot(inla.all)$Baseline
plot(inla.all)$Outcomes



###### ##@##
print(summary.inlaJOINT)
getAnywhere("summary.inla")
getAnywhere("plot.inla")





##### OLD #######


### Basics ###

summary(inla3j)

plot(inla3j)$Baseline
plot(inla3j)$Outcomes

### Visualise ###
Wb <- function(t, alpha, lambda) {
  alpha * t**(alpha-1) * lambda
}

t <- seq(0,200, by = 1)

plot(t, Wb(t, 1.5, 0.01))

inla3j$summary.fixed["Intercept_S1", "mean"]
inla3j$summary.hyperpar$mean

#hazards

r31 <- Wb(t, inla3j$summary.hyperpar$mean[1], exp(inla3j$summary.fixed["Intercept_S1", "mean"]))
r32 <- Wb(t, inla3j$summary.hyperpar$mean[2], exp(inla3j$summary.fixed["Intercept_S2", "mean"]))
r34 <- Wb(t, inla3j$summary.hyperpar$mean[3], exp(inla3j$summary.fixed["Intercept_S3", "mean"]))
r35 <- Wb(t, inla3j$summary.hyperpar$mean[4], exp(inla3j$summary.fixed["Intercept_S4", "mean"]))
r36 <- Wb(t, inla3j$summary.hyperpar$mean[5], exp(inla3j$summary.fixed["Intercept_S5", "mean"]))
r37 <- Wb(t, inla3j$summary.hyperpar$mean[6], exp(inla3j$summary.fixed["Intercept_S6", "mean"]))
r38 <- Wb(t, inla3j$summary.hyperpar$mean[7], exp(inla3j$summary.fixed["Intercept_S7", "mean"]))
r39 <- Wb(t, inla3j$summary.hyperpar$mean[8], exp(inla3j$summary.fixed["Intercept_S8", "mean"]))


plot(t, r31, pch=19, cex = 0.5)
r31.u <- Wb(t, inla3j$summary.hyperpar$mean[1], exp(inla3j$summary.fixed["Intercept_S1", "0.975quant"]))
r31.l <- Wb(t, inla3j$summary.hyperpar$mean[1], exp(inla3j$summary.fixed["Intercept_S1", "0.025quant"]))
r31.uu <- Wb(t, inla3j$summary.hyperpar$"0.975quant"[1], exp(inla3j$summary.fixed["Intercept_S1", "0.975quant"]))
r31.ll <- Wb(t, inla3j$summary.hyperpar$"0.025quant"[1], exp(inla3j$summary.fixed["Intercept_S1", "0.025quant"]))
points(t, r31.u, pch=19, cex = 0.5, col="red")
points(t, r31.l, pch=19, cex = 0.5, col="red")
points(t, r31.uu, pch=19, cex = 0.5, col="blue")
points(t, r31.ll, pch=19, cex = 0.5, col="blue")

plot_hazard <- function(tr, years) {
  t <- seq(0, years, by = 1)
  r <- Wb(t, inla3j$summary.hyperpar$mean[tr], exp(inla3j$summary.fixed[paste0("Intercept_S", as.character(tr)), "mean"]))
  
  r.u <- Wb(t, inla3j$summary.hyperpar$mean[tr], exp(inla3j$summary.fixed[paste0("Intercept_S", as.character(tr)), "0.975quant"]))
  r.l <- Wb(t, inla3j$summary.hyperpar$mean[tr], exp(inla3j$summary.fixed[paste0("Intercept_S", as.character(tr)), "0.025quant"]))
  r.uu <- Wb(t, inla3j$summary.hyperpar$"0.975quant"[tr], exp(inla3j$summary.fixed[paste0("Intercept_S", as.character(tr)), "0.975quant"]))
  r.ll <- Wb(t, inla3j$summary.hyperpar$"0.025quant"[tr], exp(inla3j$summary.fixed[paste0("Intercept_S", as.character(tr)), "0.025quant"]))
  
  plot(t, r, pch=19, cex = 0.5)
  points(t, r.u, pch=19, cex = 0.5, col="red")
  points(t, r.l, pch=19, cex = 0.5, col="red")
  points(t, r.uu, pch=19, cex = 0.5, col="blue")
  points(t, r.ll, pch=19, cex = 0.5, col="blue")
}

plot_hazard(1, 100)


## Full equation

data_msm <- read.csv(here("Data", "BTE", "bte_msm_ready.csv"))[,-1]

data_msm_4b <- read.csv(here("Data", "BTE", "bte_msm_ready.csv"))[,-1] %>% 
  filter(SREG_ECO %in% c("4bM", "4bS", "4bT"))


plot_hazard_full <- function(tr, years) {
  t <- seq(0, years, by = 1)
  r <- Wb(t, inla3j$summary.hyperpar$mean[tr], exp(inla3j$summary.fixed[paste0("Intercept_S", as.character(tr)), "mean"] +
                                                     inla3j$summary.fixed[paste0("cov_Tmean_S", as.character(tr)), "mean"] * mean(data_msm_4b$cov_Tmean, na.rm = T) +
                                                     inla3j$summary.fixed[paste0("cov_CMI_S", as.character(tr)), "mean"] * mean(data_msm_4b$cov_CMI, na.rm = T) +
                                                     inla3j$summary.fixed[paste0("cov_soil_S", as.character(tr)), "mean"] * mean(data_msm_4b$cov_soil, na.rm = T)))
  plot(t, r, pch=19, cex = 0.5)
}

plot_hazard_full(8, 100)











h31 <- Wb(t, inla3j$summary.hyperpar$mean[1], 
   exp(inla3j$summary.fixed["Intercept_S1", "mean"] + 
         inla3j$summary.fixed["cov_Tmean_S1", "mean"] * mean(data_msm_4b$cov_Tmean, na.rm = T) +
         inla3j$summary.fixed["cov_CMI_S1", "mean"] * mean(data_msm_4b$cov_CMI, na.rm = T)+ 
         inla3j$summary.fixed["cov_soil_S1", "mean"] * mean(data_msm_4b$cov_soil, na.rm = T)))

plot(t, h31)





###
cumsum(r31)

plot(t, cumsum(r31), pch=19, cex = 0.5)
plot(t, cumsum(r31), pch=19, cex = 0.5)

plot(t, p33, pch=19, cex = 0.5)
plot(t, p31, pch=19, cex = 0.5)

###probablities
p33 <- exp(-cumsum(r31)-cumsum(r32)-cumsum(r34)-cumsum(r35)-cumsum(r36)-cumsum(r37)-cumsum(r38)-cumsum(r39))
p31 <- exp(r31*p33)







plot(t, r32, pch=19, cex = 0.5)
plot(t, r34, pch=19, cex = 0.5)
plot(t, r35, pch=19, cex = 0.5)
plot(t, r36, pch=19, cex = 0.5)
plot(t, r37, pch=19, cex = 0.5)
plot(t, r38, pch=19, cex = 0.5)
plot(t, r39, pch=19, cex = 0.5)

plot(t, Wb(t, inla3j$summary.hyperpar$mean[1], 1))
plot(t, Wb(t, inla3j$summary.hyperpar$mean[2], 1))
plot(t, Wb(t, inla3j$summary.hyperpar$mean[3], 1))
plot(t, Wb(t, inla3j$summary.hyperpar$mean[4], 1))
plot(t, Wb(t, inla3j$summary.hyperpar$mean[5], 1))
plot(t, Wb(t, inla3j$summary.hyperpar$mean[6], 1))
plot(t, Wb(t, inla3j$summary.hyperpar$mean[7], 1))
plot(t, Wb(t, inla3j$summary.hyperpar$mean[8], 1))

summary(inla3j)
inla3j$waic$waic
inla3j$dic$dic

plot_inla_coef <- function(inla.out) {
  plots <- list()
  par <- inla.out$summary.fixed
  hyper <- inla.out$summary.hyperpar
  
  coef <- setNames(data.frame(matrix(ncol = 5, nrow = nrow(par) + nrow(hyper))), c("name", "mean", "lb", "ub", "hyperpar"))
  
  for(i in 1:nrow(par)) {
    coef[i,] <- c(rownames(par)[i], par[i, "mean"], par[i, "0.025quant"], par[i, "0.975quant"], FALSE)
  }
  
  for(j in 1:nrow(hyper)) {
    coef[nrow(par)+j,] <- c(rownames(hyper)[j], hyper$mean[j], hyper$"0.025quant"[j], hyper$"0.975quant"[j], TRUE)
  }
  coef[, "hyperpar"] <- as.logical(coef[, "hyperpar"])
  print(coef)
  
  ##Add per from transition group
  
  plots[["all"]] <- ggplot(data=coef, aes(x=1:nrow(coef))) +
    geom_point(aes(y=mean)) +
    geom_errorbar(aes(ymin=lb, ymax=ub))
  
  plots[["coef"]] <- ggplot(data=coef %>% filter(hyperpar==FALSE), aes(x=name)) +
    geom_point(aes(y=mean)) +
    geom_errorbar(aes(ymin=lb, ymax=ub))
  
  plots[["hyper"]] <- ggplot(data=coef %>% filter(hyperpar==TRUE), aes(x=name)) +
    geom_point(aes(y=mean)) +
    geom_errorbar(aes(ymin=lb, ymax=ub))
  
  plots
}

a <- plot_inla_coef(inla3j)
a[["coef"]]






a[["all"]]
a[["hyper"]]

inla3j$names.fixed

priors.used(inla3j)
summary(inla3j)

plot(exp(cumsum(r31)))
plot(exp(cumsum(r32)))

