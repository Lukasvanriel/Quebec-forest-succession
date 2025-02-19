#### Load packages ####
library(tidyverse)
library(INLA)
library(INLAjoint)
library(here)

### Functions ####
set.back.inla <- function(string){
  object <- list()
  
  file.list <- list.files(path=here("Data-Output", "INLA", string), pattern=".RDS", all.files=TRUE, 
                          full.names=F)
  # skip .args.RDS (too heavy for my laptop; fix later if necessary)
  if(paste0(string, "_.args.RDS") %in% file.list) file.list <- file.list[file.list != paste0(string, "_.args.RDS")]
  
  for(file in file.list) {
    file <- substr(file, nchar(string) + 2, nchar(file))
    file <- substr(file, 1, nchar(file) - 4)
    object[[file]] <- readRDS(here("Data-Output", "INLA", string, paste0(string, "_", file, ".RDS")))
    print(names(object)[length(object)])
  }
  
  class(object) <- c("INLAjoint", "inla")
  return(object)
}

Wb <- function(t, X, l, a, b) {
  l * a * t^(a - 1) * exp(sum(b * X))
}

h <- function(string, time, params, Nstates){
  f <- match(substr(string, 1,1), LETTERS)
  t <- match(substr(string, 2,2), LETTERS)
  
  if(t == f) {return(-1)}
  
  helper <- seq(1, Nstates, by=1)
  transitions <- c()
  for(i in helper) {
    for(j in helper[-i]) {
      transitions <- c(transitions, 10*i + j%%10)
    }
  }
  info.transitions <- data.frame(inla_trans=seq(1:(Nstates*(Nstates-1))), from = transitions %/% 10, to = transitions %% 10)
  
  transition_number <- info.transitions %>% 
    filter(from==f, to==t) %>% 
    pull(inla_trans)
  
  # Based on the transition number, lets now get the necessary parameters to feed to the Weibull hazard
  
  X <- params$X
  l <- exp(params$beta[names(params$beta) %in% paste0(c("Interc"), transition_number)])
  a <- params$alpha[transition_number]
  b <- params$beta[names(params$beta) %in% paste0(c("Tmean", "CMI", "Fire", "Harv", "Pest"), transition_number)]
  
  Wb(time, X, l, a, b)
}

generate_Q <- function(t, params, Nstates){
  all_transitions <- paste0(rep(LETTERS[1:Nstates], each=Nstates), rep(LETTERS[1:Nstates], Nstates))
  h_all <- lapply(all_transitions, FUN = function(x) h(x, t, params, Nstates))
  h_all <- matrix(h_all, ncol = Nstates, nrow = Nstates, byrow=TRUE)
  
  sapply(1:Nstates, function(i) {
    h_all[i,i] <<- list(-Reduce(`+`, lapply(seq(1, Nstates)[-i], function(j) h_all[i, j][[1]])))
    })
  
  h_all
}

inla_trans_prob <- function(parameters, tmax, N){
  Qall <- generate_Q(t = 0:tmax, parameters, Nstates = 9)
  Q <- lapply(1:(length(Qall[[1]])), function(i) {
    matrix(sapply(Qall, function(m) m[[i]]), nrow = nrow(Qall), ncol = ncol(Qall)) })
  
  P <- vector("list", length = length(tmax))
  p <- diag(N)
  sapply(1:(tmax + 1), function(x) {
    P[[x]] <<- p %*% expm(Q[[x]])
    p <<- P[[x]]
  } )
  
  p.vector <- lapply(P, function(p) {
    p <- as.vector(t(p))
    names(p) <- paste0(rep(LETTERS[1:N], each=N), rep(LETTERS[1:N], N))
    p
  } )
  
  bind_rows(p.vector)
}

### Other ####
info.classes <- data.frame(Class=LETTERS[1:9], Name=c("Paper birch", "Shade int.", "Yellow birch",
                                               "Maples", "Other dec.", "Bals. fir",
                                               "B/R spruce", "Jack pine", "Other con."))

### Body ####
inla.allj.T.C.Pc <- set.back.inla("allj.T.C.Pc")
summary(inla.allj.T.C.Pc)

pars <- inla.allj.T.C.Pc$summary.fixed

## Create parameter object
info.classes
helper <- seq(1, 9, by=1)
transitions <- c()
for(i in helper) {
  for(j in helper[-i]) {
    transitions <- c(transitions, 10*i + j%%10)
  }
}
info.transitions <- data.frame(inla_trans=seq(1:(9*(9-1))), from = transitions %/% 10, to = transitions %% 10)

tmax <- 60

parameters <- list()

parameters$alpha <- inla.allj.T.C.Pc$summary.hyperpar$mean
names(parameters$alpha) <- paste0("alpha", as.character(1:length(parameters$alpha)))

parameters$beta <- inla.allj.T.C.Pc$summary.fixed$mean
names(parameters$beta) <- paste0(rep(c("Interc", "Tmean", "CMI", "Fire", "Harv", "Pest"), each=1, length(parameters$alpha)),
                                 rep(1:length(parameters$alpha), each=6))

parameters$X <- rep(0,5)

## Generate the necessary matrix objects:

parameters$X <- rep(0,5)
p.df0 <- inla_trans_prob(parameters, tmax = 60, N = 9)

parameters$X <- c(0,0,1,0,0)
p.df1 <- inla_trans_prob(parameters, tmax = 60, N = 9)

parameters$X <- c(0,0,0,1,0)
p.df2 <- inla_trans_prob(parameters, tmax = 60, N = 9)

parameters$X <- c(0,0,0,0,1)
p.df3 <- inla_trans_prob(parameters, tmax = 60, N = 9)


### Plot some stuff:

plot(0:tmax, p.df0$AA, pch=20, ylim=c(0,1))
points(0:tmax, p.df0$AB, pch=20, col="red")
points(0:tmax, p.df0$AC, pch=20, col="red")
points(0:tmax, p.df0$AD, pch=20, col="red")
points(0:tmax, p.df0$AE, pch=20, col="red")
points(0:tmax, p.df0$AF, pch=20, col="red")
points(0:tmax, p.df0$AG, pch=20, col="red")
points(0:tmax, p.df0$AH, pch=20, col="red")
points(0:tmax, p.df0$AI, pch=20, col="red")


plot(0:tmax, p.df0$AA, pch=20, ylim=c(0,1))
points(0:tmax, p.df1$AA, pch=20, col="red")
points(0:tmax, p.df2$AA, pch=20, col="blue")
points(0:tmax, p.df3$AA, pch=20, col="green")


plot(0:tmax, p.df0$AA, pch=20, ylim=c(0,1))
points(0:tmax, p.df1$AA, pch=20, col="red")
points(0:tmax, p.df2$AA, pch=20, col="blue")
points(0:tmax, p.df3$AA, pch=20, col="green")

plot(0:tmax, p.df0$FH, pch=20, ylim=c(0,1))
points(0:tmax, p.df1$FH, pch=20, col="red")
points(0:tmax, p.df2$FH, pch=20, col="blue")
points(0:tmax, p.df3$FH, pch=20, col="green")

plot(0:tmax, p.df0$FG, pch=20, ylim=c(0,1))
points(0:tmax, p.df1$FG, pch=20, col="red")
points(0:tmax, p.df2$FG, pch=20, col="blue")
points(0:tmax, p.df3$FG, pch=20, col="green")



pars[(46*6 + 1):(46*6 + 6),]

info.classes
info.transitions








#pdf.backup

a <- generate_Q(t = 0:tmax, parameters, Nstates = 9)

Q <- lapply(1:(length(a[[1]])), function(i) {matrix(sapply(a, function(m) m[[i]]), nrow = nrow(a), ncol = ncol(a))})

P <- vector("list", length = length(tmax))
p <- diag(9)
sapply(1:(tmax + 1), function(x) {
  P[[x]] <<- p %*% expm(Q[[x]])
  p <<- P[[x]]
} )

p.vector <- lapply(P, function(p) {
  p <- as.vector(t(p))
  Nstates <- 9
  names(p) <- paste0(rep(LETTERS[1:Nstates], each=Nstates), rep(LETTERS[1:Nstates], Nstates))
  p
} )

p.df <- bind_rows(p.vector)
#### OLD #########
#outdated way to create summary myself
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
#inl.summary(inla.summ, inla.hype, inla.mlik, inla.dic, inla.waic, inla.cpuu)

##Manual way to restore split inla output

set.back.inla <- function(string){
  object <- list()
  
  file.list <- list.files(path=here("Data-Output", "INLA", string), pattern=".RDS", all.files=TRUE, 
                          full.names=F)
  # skip .args.RDS (too heavy for my laptop; fix later if necessary)
  if(paste0(string, "_.args.RDS") %in% file.list) file.list <- file.list[file.list != paste0(string, "_.args.RDS")]
  
  for(file in file.list) {
    file <- substr(file, nchar(string) + 2, nchar(file))
    file <- substr(file, 1, nchar(file) - 4)
    object[[file]] <- readRDS(here("Data-Output", "INLA", string, paste0(string, "_", file, ".RDS")))
    print(names(object)[length(object)])
  }
 
  class(object) <- c("INLAjoint", "inla")
  return(object)
}

### Documentation ####
if(F) {
?joint()
inla.doc("weibullsurv")
vignette("INLAjoint")

print(summary.inlaJOINT)
getAnywhere("summary.inla")
getAnywhere("plot.inla")
}
### Data ####

inla.allj.T.C.Pc <- set.back.inla("allj.T.C.Pc")
BTE_inla <- set.back.inla("BTE_allj.T.C.Pc")

inla.allj.T.C.Pc$summary.fixed

### Visualise INLAjoint model output ####
### Basics: ####
# Get summary of first 9 transitions
summary(inla.allj.T.C.Pc)
summary(BTE_inla)

?msm::pmatrix.msm()

# Create plots of marginal posteriors
pl.hyp <- plot(inla.allj.T.C.Pc)$Baseline
pl.cov <- plot(inla.allj.T.C.Pc)$Outcomes

pl.cov$S22
pl.cov$S2
pl.hyp

# Gather all parameters in dataframes
inla.allj.T.C.Pc$summary.fixed
inla.allj.T.C.Pc$summary.hyperpar
hist(inla.allj.T.C.Pc$summary.hyperpar$mean)

# Can I recover the scale parameters from the intercepts? YES!
exp(inla.allj.T.C.Pc$summary.fixed$mean[rownames(inla.allj.T.C.Pc$summary.fixed)=="Intercept_S1"])
exp(inla.allj.T.C.Pc$summary.fixed$mean[rownames(inla.allj.T.C.Pc$summary.fixed)=="Intercept_S2"])
exp(inla.allj.T.C.Pc$summary.fixed$mean[rownames(inla.allj.T.C.Pc$summary.fixed)=="Intercept_S3"])
exp(inla.allj.T.C.Pc$summary.fixed$mean[rownames(inla.allj.T.C.Pc$summary.fixed)=="Intercept_S9"])

summary(inla.allj.T.C.Pc$summary.hyperpar$mean)

### Advanced: ####
# Create helper objects
info.classes <- data.frame(Nb=seq(1:9), Name=c("Paper birch", "Shade int.", "Yellow birch",
                                               "Maples", "Other dec.", "Bals. fir",
                                               "B/R spruce", "Jack pine", "Other con."))

helper <- seq(1,9, by=1)
transitions <- c()
for(i in helper) {
  for(j in helper[-i]) {
    transitions <- c(transitions, 10*i + j%%10)
  }
}
info.transitions <- data.frame(inla_trans=seq(1:72), from=transitions %/% 10, to=transitions %% 10)

## Function to create dataframe containing all the transition parameters
# Transitions need to be Weibull for this to work I think
extract_inla_param <- function(inla.obj){
  output.df <- data.frame(inla_trans=seq(1:inla.obj$nhyper))
  
  rownames(inla.obj$summary.hyperpar) <- NULL
  hyperp <- apply(inla.obj$summary.hyperpar, 1, function(x) as.list(x))
  output.df$alpha <- hyperp
  
  N.par <- nrow(inla.obj$summary.fixed) / inla.obj$nhyper
  rownames(inla.obj$summary.fixed) <- sapply(rownames(inla.obj$summary.fixed), FUN=function(x) {ifelse(substring(x, nchar(x) - 1, nchar(x) - 1) == "S",
                                                                                                       paste0(substring(x, 1, nchar(x)-1), "0", substring(x, nchar(x), nchar(x))),
                                                                                                       x)} )
  param <- names(table(substring(rownames(inla.obj$summary.fixed), 1, nchar(rownames(inla.obj$summary.fixed)) - 4)))
  param <- c("Intercept", param[-match("Intercept", param)])
  
  for(i in 1:N.par) {
    sub.param.df <- inla.obj$summary.fixed[grepl(param[i], rownames(inla.obj$summary.fixed)),]
    p <- apply(sub.param.df, 1, function(x) as.list(x))
    output.df[[ncol(output.df) + 1]] <- p
    colnames(output.df)[ncol(output.df)] <- param[i]
  }
  
  output.df
}

inla.params <- extract_inla_param(inla.allj.T.C.Pc)


# Join the dataframes containing from/to info and parameters
full.df <- right_join(info.transitions, inla.params)

#write_rds(full.df, here("Data-Output", "INLA", "allj.T.C.Pc", "full_df.RDS"))
tester <- read_rds(here("Data-Output", "INLA", "allj.T.C.Pc", "full_df.RDS"))

 ## Plot parameters with CI:
## Function that plots the individual coefficients for the specified parameter,
## starting from the indicated class
# inla.output requires to contain 
plot_inla_param <- function(inla.output, info.classes, fr, parameter) {
  inla.output %>% 
    filter(from==fr) %>% 
    select(from, to, all_of(parameter)) %>% 
    unnest_wider(parameter) %>% 
    rename(c("median"="0.5quant", "lquant"="0.025quant", "uquant"="0.975quant")) %>% 
    merge(rename(info.classes, c("to"="Nb", "Name_to"="Name"))) %>% 
    mutate(Name_to=factor(Name_to, levels = Name_to)) %>% 
    ggplot(aes(x=Name_to)) +
    geom_point(aes(y=median)) +
    geom_errorbar(aes(ymin=lquant, ymax=uquant), width=0.2)+
    ggtitle(paste0("coefficients for ", substring(parameter, 5, nchar(parameter)),", transitions starting from ", info.classes$Name[fr])) +
    labs(y= "", x = "to")
}

plot_inla_param(full.df, info.classes, fr= 1, "cov_Tmean")
plot_inla_param(full.df, info.classes, fr= 2, "cov_Tmean")
plot_inla_param(full.df, info.classes, fr= 3, "cov_Tmean")

plot_inla_param(full.df, info.classes, fr= 5, "cov_CMI")

plot_inla_param(full.df, info.classes, fr= 5, "alpha")

plot_inla_param(full.df, info.classes, fr= 9, "Intercept")


### ####

plot_inla_covariates <- function(inla.output, info.classes, fr, d = 0.3) {
  dodge <- position_dodge(width = d)
  
  inla.output %>% 
    filter(from==fr) %>% 
    select(from, to, colnames(full.df)[grepl("cov_", colnames(full.df))]) %>% 
    unnest_wider(colnames(full.df)[grepl("cov_", colnames(full.df))], names_sep = ".") %>% 
    rename_with(~str_replace(., '0.025quant', 'lquant')) %>% 
    rename_with(~str_replace(., '0.975quant', 'uquant')) %>% 
    rename_with(~str_replace(., '0.5quant', 'median')) %>% 
    merge(rename(info.classes, c("to"="Nb", "Name_to"="Name"))) %>% 
    mutate(Name_to=factor(Name_to, levels = Name_to)) %>% 
    pivot_longer(cols = colnames(.)[grepl("cov_", colnames(.))]) %>% 
    separate(name, into=c("covariate", "param"), sep="\\.") %>%
    filter(param %in% c("median", "uquant", "lquant")) %>% 
    pivot_wider(names_from = param) %>% 
    ggplot(aes(x=Name_to)) +
    geom_point(aes(x=Name_to, y=median, colour=covariate), position = dodge) +
    geom_errorbar(aes(ymin=lquant, ymax=uquant, colour=covariate), position=dodge, width=0.2) +
    ggtitle(paste("Covariate parameters for transitions starting from", info.classes$Name[fr])) +
    labs(y= "", x = "to")
}

plot_inla_covariates(full.df, info.classes, 1)

plot_inla_covariates(full.df, info.classes, 9)

png(filename="~/Desktop/covariates9.png", width = 1500, height = 1100)
plot_inla_covariates(full.df, info.classes, 9)
dev.off()

## Plot hazards
years <- seq(0,60, by = 1)

Wb <- function(t, alpha, lambda) {
  alpha * t**(alpha-1) * lambda
}
Wb_H <- function(t, alpha, lambda) {
  t**(alpha) * lambda
}

aa <- cumsum(Wb(years, 1.2, 0.01))
bb <- Wb_H(years, 1.2, 0.01)

plot(years, Wb(years, 1.5, 0.01))

# Baseline
plot(years, Wb(years, inla.allj.T.C.Pc$summary.hyperpar$mean[1], 
           exp(inla.allj.T.C.Pc$summary.fixed["Intercept_S1", "mean"])))


plot(years, Wb(years, inla.allj.T.C.Pc$summary.hyperpar$mean[1], 
           exp(inla.allj.T.C.Pc$summary.fixed["Intercept_S1", "mean"] +
                 inla.allj.T.C.Pc$summary.fixed["cov_pert_class1_S1", "mean"])))

plot(years, Wb(years, inla.allj.T.C.Pc$summary.hyperpar$mean[1], 
           exp(inla.allj.T.C.Pc$summary.fixed["Intercept_S1", "mean"] +
                 inla.allj.T.C.Pc$summary.fixed["cov_pert_class2_S2", "mean"])))

plot(years, Wb(years, inla.allj.T.C.Pc$summary.hyperpar$mean[1], 
           exp(inla.allj.T.C.Pc$summary.fixed["Intercept_S1", "mean"] +
                 inla.allj.T.C.Pc$summary.fixed["cov_pert_class3_S3", "mean"])))


# plot(t, Wb(t, inla.allj.T.C.Pc$summary.hyperpar$mean[1], 
#            exp(inla.allj.T.C.Pc$summary.fixed["Intercept_S1", "mean"] + 
#                  inla.allj.T.C.Pc$summary.fixed["cov_Tmean_S1", "mean"] * mean(Tmean) )))

## Try some transitions probabilities over time

get_risks <- function(inla.output, t=50){
  years <- seq(0, t, by = 1)
  risks <- list()
  for(i in 1:nrow(inla.output)) {
    risks[[i]] <- list()
    risks[[i]]$baseline <- Wb(years, inla.output$alpha[[i]]$mean, 
                              exp(inla.output$Intercept[[i]]$mean))
    risks[[i]]$fire <- Wb(years, inla.output$alpha[[i]]$mean, 
                          exp(inla.output$Intercept[[i]]$mean + full.df$cov_pert_class1[[i]]$mean))
    risks[[i]]$harvest <- Wb(years, inla.output$alpha[[i]]$mean, 
                             exp(inla.output$Intercept[[i]]$mean + full.df$cov_pert_class2[[i]]$mean))
    risks[[i]]$pest <- Wb(years, inla.output$alpha[[i]]$mean, 
                          exp(inla.output$Intercept[[i]]$mean + full.df$cov_pert_class3[[i]]$mean))
  }
  names(risks) <- as.character(1:nrow(inla.output))
  risks
}

risks <- get_risks(full.df, 60)


# Transition probabilities to stay in state
trans_prob_timelapse_stay <- function(inla.output, state, t=50, risks=NULL) {
  if(is.null(risks) || length(risks[[1]]$baseline) != (t+1)) risks <- get_risks(inla.output, t)
  
  trans.out <- inla.output %>% 
    filter(from == state) %>% 
    pull(inla_trans)
  
  stay.risks <- list()
  
  cumulatives <- sapply(trans.out, FUN=function(x) cumsum(risks[[x]]$baseline))
  stay.risks$baseline <- apply(as.matrix(cumulatives), 1, sum)
  cumulatives <- sapply(trans.out, FUN=function(x) cumsum(risks[[x]]$fire))
  stay.risks$fire <- apply(as.matrix(cumulatives), 1, sum)
  cumulatives <- sapply(trans.out, FUN=function(x) cumsum(risks[[x]]$harvest))
  stay.risks$harvest <- apply(as.matrix(cumulatives), 1, sum)
  cumulatives <- sapply(trans.out, FUN=function(x) cumsum(risks[[x]]$pest))
  stay.risks$pest <- apply(as.matrix(cumulatives), 1, sum)
  
  stay.risks
}

test <- trans_prob_timelapse_stay(full.df, 2, t=60)
plot(exp(-test$baseline), pch=20)
points(exp(-test$fire), pch=20, col="red")
points(exp(-test$harvest), pch=20, col="blue")
points(exp(-test$pest), pch=20, col="darkgreen")



p11 <- exp(- cumsum(risks[[1]]$baseline) - cumsum(risks[[2]]$baseline) - cumsum(risks[[3]]$baseline) - 
             cumsum(risks[[4]]$baseline) - cumsum(risks[[5]]$baseline) - cumsum(risks[[6]]$baseline) -
             cumsum(risks[[7]]$baseline) - cumsum(risks[[8]]$baseline))
p11_fire <- exp(- cumsum(risks[[1]]$fire) - cumsum(risks[[2]]$fire) - cumsum(risks[[3]]$fire) - 
                  cumsum(risks[[4]]$fire) - cumsum(risks[[5]]$fire) - cumsum(risks[[6]]$fire) -
                  cumsum(risks[[7]]$fire) - cumsum(risks[[8]]$fire))
p11_harv <- exp(- cumsum(risks[[1]]$harvest) - cumsum(risks[[2]]$harvest) - cumsum(risks[[3]]$harvest) - 
                  cumsum(risks[[4]]$harvest) - cumsum(risks[[5]]$harvest) - cumsum(risks[[6]]$harvest) -
                  cumsum(risks[[7]]$harvest) - cumsum(risks[[8]]$harvest))
p11_pest <- exp(- cumsum(risks[[1]]$pest) - cumsum(risks[[2]]$pest) - cumsum(risks[[3]]$pest) - 
                  cumsum(risks[[4]]$pest) - cumsum(risks[[5]]$pest) - cumsum(risks[[6]]$pest) -
                  cumsum(risks[[7]]$pest) - cumsum(risks[[8]]$pest))

plot(years, p11, pch=20)
points(years, p11_fire, col="red", pch=20)
points(years, p11_harv, col="blue", pch=20)
points(years, p11_pest, col="darkgreen", pch=20)
points(years, exp(-test), col="orange")

#Probabilities to change state
??pmatrix.piecewise.msm()
msm::pmatrix.piecewise.msm

cumsum(risks[[20]]$baseline)

trans_prob_timelapse_leave <- function(inla.output, state.fr, state.to, t=50, risks=NULL) {
  if(is.null(risks) || length(risks[[1]]$baseline) != (t+1)) risks <- get_risks(inla.output, t)
  
  trans.number <- full.df %>% 
    filter(from == state.fr, to == state.to) %>% 
    pull(inla_trans)
  
  leave.risks <- list()
  
  leave.risks$baseline <- cumsum(exp(-trans_prob_timelapse_stay(inla.output, state.fr, t=t, risks=risks)$baseline) * risks[[trans.number]]$baseline)
  leave.risks$fire <- cumsum(exp(-trans_prob_timelapse_stay(inla.output, state.fr, t=t, risks=risks)$fire) * risks[[trans.number]]$fire)
  leave.risks$harvest <- cumsum(exp(-trans_prob_timelapse_stay(inla.output, state.fr, t=t, risks=risks)$harvest) * risks[[trans.number]]$harvest)
  leave.risks$pest <- cumsum(exp(-trans_prob_timelapse_stay(inla.output, state.fr, t=t, risks=risks)$pest) * risks[[trans.number]]$pest)
  
  leave.risks
}

trans_prob_timelapse_leave_alv <- function(inla.output, state.fr, state.to, t=50, risks=NULL) {
  if(is.null(risks) || length(risks[[1]]$baseline) != (t+1)) risks <- get_risks(inla.output, t)
  
  trans.number <- full.df %>% 
    filter(from == state.fr, to == state.to) %>% 
    pull(inla_trans)
  
  leave.risks <- list()
  
  leave.risks$baseline <- cumsum(exp(-trans_prob_timelapse_stay(inla.output, state.fr, t=t, risks=risks)$baseline) * 
                                   risks[[trans.number]]$baseline * 
                                   exp(-trans_prob_timelapse_stay(inla.output, state.to, t=t, risks=risks)$baseline))
  leave.risks$fire <- cumsum(exp(-trans_prob_timelapse_stay(inla.output, state.fr, t=t, risks=risks)$fire) * 
                               risks[[trans.number]]$fire * 
                               exp(-trans_prob_timelapse_stay(inla.output, state.to, t=t, risks=risks)$fire))
  leave.risks$harvest <- cumsum(exp(-trans_prob_timelapse_stay(inla.output, state.fr, t=t, risks=risks)$harvest) * 
                                  risks[[trans.number]]$harvest * 
                                  exp(-trans_prob_timelapse_stay(inla.output, state.to, t=t, risks=risks)$harvest))
  leave.risks$pest <- cumsum(exp(-trans_prob_timelapse_stay(inla.output, state.fr, t=t, risks=risks)$pest) * 
                               risks[[trans.number]]$pest * 
                               exp(-trans_prob_timelapse_stay(inla.output, state.to, t=t, risks=risks)$pest))
  
  leave.risks
}

a <- trans_prob_timelapse_leave(full.df, 1, 4, t=60)
data.frame(b=a$baseline, f=a$fire, h=a$harvest, p=a$pest, t=0:60) %>%
  pivot_longer(cols = c(b, f, h, p)) %>% 
  ggplot() + geom_point(aes(x=t, y=value, color=name))


#Should this not have the exponential?
c <- trans_prob_timelapse_leave_alv(full.df, 3, 1, t=60)
data.frame(b=c$baseline, f=c$fire, h=c$harvest, p=c$pest, t=0:60) %>%
  pivot_longer(cols = c(b, f, h, p)) %>% 
  mutate(value.alt=(1-exp(-value))) %>% 
  ggplot() + geom_point(aes(x=t, y=value.alt, color=name))



#integrate()
inla.posterior.sample(inla.allj.T.C.Pc, 10)

#predict(inla_model, newdata = data.frame(covariate = covariate, time = time)

tester <- trans_prob_timelapse_leave(full.df, 1, 9, t=60)
plot(years, 1-exp(-tester$baseline))


plot(years, exp(-trans_prob_timelapse_stay(full.df, 6, t=60)$baseline), pch=20, cex =3, col="black", 
     ylab="Transition probability", xlab="time [years]", main="From balsam fir")
points(years, 1-exp(-trans_prob_timelapse_leave(full.df, 6, 2, t=60)$baseline), pch=20,cex =3, col="red")
points(years, 1-exp(-trans_prob_timelapse_leave(full.df, 6, 3, t=60)$baseline), pch=20,cex =3, col="red")
points(years, 1-exp(-trans_prob_timelapse_leave(full.df, 6, 4, t=60)$baseline), pch=20,cex =3, col="red")
points(years, 1-exp(-trans_prob_timelapse_leave(full.df, 6, 5, t=60)$baseline), pch=20,cex =3, col="red")
points(years, 1-exp(-trans_prob_timelapse_leave(full.df, 6, 1, t=60)$baseline), pch=20,cex =3, col="red")
points(years, 1-exp(-trans_prob_timelapse_leave(full.df, 6, 7, t=60)$baseline), pch=20,cex =3, col="red")
points(years, 1-exp(-trans_prob_timelapse_leave(full.df, 6, 8, t=60)$baseline), pch=20,cex =3, col="red")
points(years, 1-exp(-trans_prob_timelapse_leave(full.df, 6, 9, t=60)$baseline), pch=20,cex =3, col="red")

plot(years, 1-exp(-trans_prob_timelapse_leave(full.df, 1, 2, t=60)$baseline), pch=20, cex =3, col="black", 
     ylab="Transition probability", xlab="time [years]", main="Paper birch to  Other conif")
points(years, 1-exp(-trans_prob_timelapse_leave(full.df, 1, 2, t=60)$fire), pch=20, cex =3, col="red")
points(years, 1-exp(-trans_prob_timelapse_leave(full.df, 1, 2, t=60)$harvest), pch=20, cex =3, col="darkgreen")
points(years, 1-exp(-trans_prob_timelapse_leave(full.df, 1, 2, t=60)$pest), pch=20, cex =3, col="purple")


png(filename="~/Desktop/timelaps71.png", width = 800, height = 600)
plot(years, 1-exp(-trans_prob_timelapse_leave(full.df, 7, 1, t=60)$baseline), pch=20, col="black", 
     main=paste0(info.classes$Name[7], " to ", info.classes$Name[1]))
points(years, 1-exp(-trans_prob_timelapse_leave(full.df, 7, 1, t=60)$fire), pch=20, col="red")
points(years, 1-exp(-trans_prob_timelapse_leave(full.df, 7, 1, t=60)$harvest), pch=20, col="blue")
points(years, 1-exp(-trans_prob_timelapse_leave(full.df, 7, 1, t=60)$pest), pch=20, col="darkgreen")
dev.off()

png(filename="~/Desktop/timelapse1.png", width = 800, height = 600)
plot(years, exp(-trans_prob_timelapse_stay(full.df, 1, t=60)$baseline), pch=20, col="black", 
     main=paste0("Probability to stay in ", info.classes$Name[1]))
points(years, exp(-trans_prob_timelapse_stay(full.df, 1, t=60)$fire), pch=20, col="red")
points(years, exp(-trans_prob_timelapse_stay(full.df, 1, t=60)$harvest), pch=20, col="blue")
points(years, exp(-trans_prob_timelapse_stay(full.df, 1, t=60)$pest), pch=20, col="darkgreen")
dev.off()


exp(-trans_prob_timelapse_stay(full.df, 1, t=60)$baseline)[11]
exp(-trans_prob_timelapse_stay(full.df, 2, t=60)$baseline)[11]
exp(-trans_prob_timelapse_stay(full.df, 3, t=60)$baseline)[11]
exp(-trans_prob_timelapse_stay(full.df, 4, t=60)$baseline)[11]
exp(-trans_prob_timelapse_stay(full.df, 5, t=60)$baseline)[11]
exp(-trans_prob_timelapse_stay(full.df, 6, t=60)$baseline)[11]
exp(-trans_prob_timelapse_stay(full.df, 7, t=60)$baseline)[11]
exp(-trans_prob_timelapse_stay(full.df, 8, t=60)$baseline)[11]
exp(-trans_prob_timelapse_stay(full.df, 9, t=60)$baseline)[11]

# test
plot(years, exp(-trans_prob_timelapse_stay(full.df, 9, t=60, risks)$baseline))
plot(years, cumsum(trans_prob_timelapse_stay(full.df, 1, t=60, risks)$baseline *
                     risks[[1]]$baseline *
                     trans_prob_timelapse_stay(full.df, 2, 60, risks)$baseline))

plot(years, cumsum(trans_prob_timelapse_stay(full.df, 1, t=60, risks)$baseline *
                     risks[[1]]$baseline))

aa <- trans_prob_timelapse_stay(full.df, 1, t=60, risks)$baseline
bb <- Wb_H(years, inla.output$alpha[[1]]$mean, exp(inla.output$Intercept[[1]]$mean)) + Wb_H(years, inla.output$alpha[[2]]$mean, exp(inla.output$Intercept[[2]]$mean))  + Wb_H(years, inla.output$alpha[[3]]$mean, exp(inla.output$Intercept[[3]]$mean)) + Wb_H(years, inla.output$alpha[[6]]$mean, exp(inla.output$Intercept[[6]]$mean))+Wb_H(years, inla.output$alpha[[4]]$mean, exp(inla.output$Intercept[[4]]$mean)) + Wb_H(years, inla.output$alpha[[7]]$mean, exp(inla.output$Intercept[[7]]$mean))+ Wb_H(years, inla.output$alpha[[5]]$mean, exp(inla.output$Intercept[[5]]$mean)) + Wb_H(years, inla.output$alpha[[8]]$mean, exp(inla.output$Intercept[[8]]$mean))
bb
aa
exp(-aa)
exp(-bb)
plot(0:60, exp(-aa))
points(0:60, exp(-bb), col="red")

###Test from paper:
p_stay <- function(inla.output, t0, t1, state, risks) {
  t <- seq(t0, t1, by=1)
  trans.out <- inla.output %>% 
    filter(from == state) %>% 
    pull(inla_trans)
  
  ind.risks <- sapply(trans.out, FUN=function(x) risks[[x]]$baseline[t0+1])
  p.tot <- 1-sum(ind.risks) 
  
  if(t0==t1) {return(p.tot)} else{p.tot <- p.tot * p_stay(inla.output, t0+1, t1, state, risks)}
}

plot(0:60, sapply(0:60, FUN=function(x) p_stay(full.df, 0, x, state = 1, risks)))
points(0:60, exp(-trans_prob_timelapse_stay(inla.output = full.df, state = 1, t = 60, risks = risks)$baseline), col="red")





tt <- p_stay(full.df, 0, 10, state = 1, risks)

print(p_stay(full.df, 0, 10, state = 9, risks))




xaa <- 1 - sum(risks[[1]]$baseline[1], risks[[2]]$baseline[1], risks[[3]]$baseline[1], risks[[4]]$baseline[1], 
        risks[[5]]$baseline[1], risks[[6]]$baseline[1], risks[[7]]$baseline[1], risks[[8]]$baseline[1])
bb <- 1 - sum(risks[[1]]$baseline[2], risks[[2]]$baseline[2], risks[[3]]$baseline[2], risks[[4]]$baseline[2], 
        risks[[5]]$baseline[2], risks[[6]]$baseline[2], risks[[7]]$baseline[2], risks[[8]]$baseline[2])
cc <- 1 - sum(risks[[1]]$baseline[3], risks[[2]]$baseline[3], risks[[3]]$baseline[3], risks[[4]]$baseline[3], 
        risks[[5]]$baseline[3], risks[[6]]$baseline[3], risks[[7]]$baseline[3], risks[[8]]$baseline[3])
  
aa*bb*cc

exp(-aa)
exp(-bb)

plot(years, exp(-aa))
points(years, exp(-bb), col="red")

## Create matrix of P over T

p.matr <- matrix(0, 9, 9)
risks <- get_risks(full.df, t=50)
for(i in 1:nrow(p.matr)){
  for(j in 1:ncol(p.matr)){
    if(i == j) {p.matr[i,j] <- exp(-trans_prob_timelapse_stay(full.df, j, t=60, risks)$baseline)[11]
    } else {p.matr[i,j] <- trans_prob_timelapse_leave(full.df, i, j, t=60, risks)$baseline[11]}
  }
}

m <- apply(p.matr, 1, sum)

heatmap(p.matr, symm = T)

p.matr <- matrix(0, 9, 9)
risks <- get_risks(full.df, t=50)
for(i in 1:nrow(p.matr)){
  for(j in 1:ncol(p.matr)){
    if(i == j) {p.matr[i,j] <- exp(-trans_prob_timelapse_stay(full.df, j, t=60, risks)$baseline)[11]
    } else {p.matr[i,j] <- trans_prob_timelapse_leave_alv(full.df, i, j, t=60, risks)$baseline[11]}
  }
}

m <- apply(p.matr, 1, sum)

p.matr.exp <- matrix(0, 9, 9)
risks <- get_risks(full.df, t=50)
for(i in 1:nrow(p.matr.exp)){
  for(j in 1:ncol(p.matr.exp)){
    if(i == j) {p.matr.exp[i,j] <- exp(-trans_prob_timelapse_stay(full.df, j, t=60, risks)$baseline)[11]
    } else {p.matr.exp[i,j] <- 1-exp(-trans_prob_timelapse_leave_alv(full.df, i, j, t=60, risks)$baseline[11])}
  }
}
m <- apply(p.matr.exp, 1, sum)





##### OLD #######

### Visualise ###
r31 <- Wb(t, inla3j$summary.hyperpar$mean[1], exp(inla3j$summary.fixed["Intercept_S1", "mean"]))

plot(t, Wb(t, 1.5, 0.01))

inla3j$summary.fixed["Intercept_S1", "mean"]
inla3j$summary.hyperpar$means

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


