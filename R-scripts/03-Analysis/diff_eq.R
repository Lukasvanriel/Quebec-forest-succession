#### Load packages ####
library(tidyverse)
library(INLA)
library(INLAjoint)
library(here)
library(deSolve)

### Functions ####

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

## Define the hazard functions


dP_dt <- function(from) {}

### Data ####
#Read in inla object
inla.allj.T.C.Pc <- set.back.inla("allj.T.C.Pc")

inla.allj.T.C.Pc$aic$aic
### Analysis ####
# Set initial conditions
P0 <- c(P_A = 1, P_B = 0, P_C = 0, 
        P_D = 0, P_E = 0, P_F = 0,
        P_G = 0, P_H = 0, P_I = 0)  # Starting in state A

t <- seq(0,60, by=1)

parameters <- list()

parameters$alpha <- inla.allj.T.C.Pc$summary.hyperpar$mean
names(parameters$alpha) <- as.character(1:72)
#names(parameters$alpha) <- names(parameters$alpha

parameters$beta <- inla.allj.T.C.Pc$summary.fixed$mean
names(parameters$beta) <- paste0(rep(c("I", "T", "C", "f", "h", "p"), each=1, 72), rep(1:72, each=6))

parameters$X <- 0
#names(parameters$X)

# Define the system of differential equations

#Needs at least alpha, beta and X in params
h <- function(string, time, params){
  #convert the letters to a transition number between 1 and 72
  f <- match(substr(string, 1,1), LETTERS)
  t <- match(substr(string, 2,2), LETTERS)
  
  helper <- seq(1,9, by=1)
  transitions <- c()
  for(i in helper) {
    for(j in helper[-i]) {
      transitions <- c(transitions, 10*i + j%%10)
    }
  }
  info.transitions <- data.frame(inla_trans=seq(1:72), from=transitions %/% 10, to=transitions %% 10)
  
  transition_number <- info.transitions %>% 
    filter(from==f, to==t) %>% 
    pull(inla_trans)
  
  print(transition_number)
  
  #BAsed on the transition number, lets now get the necessary parameters to feed to the Weibull hazard
  
  X <- 0
  l <- exp(params$beta[names(params$beta) %in% paste0(c("I"), transition_number)])
  a <- params$alpha[transition_number]
  b <- params$beta[names(params$beta) %in% paste0(c("T", "C", "f", "h", "p"), transition_number)]
  
  Wb(time, X, l, a, b)
}

Wb <- function(t, X, l, a, b) {
  l * a * t^(a - 1) * exp(sum(X * b))
}

# Params needed: lambda, alpha and beta and X
dP_dt <- function(from, t, P0, params){
  P <- P0
  
  transitions <- setdiff(LETTERS[1:9], from)
  
  leave <- sum(sapply(paste0(from, transitions), h), t=t, params=params) * P[from]
  return <- sum(mapply(function(x, y) h(x,t, params) * P[y], paste0(transitions, f), transitions, SIMPLIFY = T,))
  
  -leave + return
}





#OLDER VERSION
# dP_dt <- function(from, t, P0, params){
#   P <- P0
#   
#   transitions <- setdiff(LETTERS[1:9], from)
#   
#   leave <- sum(sapply(paste0(from, transitions), h), t=t, params=params) * P[from]
#   return <- sum(mapply(function(x, y) h(x) * P(y), paste0(t_l, f), t_l,
#                       MoreArgs = as.list(c(t, params)), SIMPLIFY = T,))
#   
#   -leave + return
# }

dP <- function(t, P0, params) {
  dPA <- dP_dt("A", t, P0, params)
  dPB <- dP_dt("B", t, P0, params)
  dPC <- dP_dt("C", t, P0, params)
  dPD <- dP_dt("D", t, P0, params)
  dPE <- dP_dt("E", t, P0, params)
  dPF <- dP_dt("F", t, P0, params)
  dPG <- dP_dt("G", t, P0, params)
  dPH <- dP_dt("H", t, P0, params)
  dPI <- dP_dt("I", t, P0, params)
  
  list(c(dPA, dPB, dPC, dPD, dPE, dPF, dPG, dPH, dPI))
}

out <- ode(y = P0, times = t, func = dP, parms = c(t=t, P0=P0, params=parameters))
out <- ode(y = P0, times = t, func = dP, parms = parameters)

?ode()

# ODE seems to require a specific format!:

dP_dt_full <- function(t, P0, params) {
  #PREP: define initial conditions and get all get all hij(t)
  P_A <- P0[1]
  P_B <- P0[2]
  P_C <- P0[3]
  P_D <- P0[4]
  P_E <- P0[5]
  P_F <- P0[6]
  P_G <- P0[7]
  P_H <- P0[8]
  P_I <- P0[9]
  
  all_transitions <- paste0(rep(LETTERS[1:9], each=9), rep(LETTERS[1:9], 9))
  # h_all <- sapply(all_transitions, FUN = function(x) {ifelse(substring(x, 1, 1) == substring(x, 2, 2), 
  #                                                            return(0), return(h(x, t, parameters)[2]))})
  h_all <- sapply(all_transitions, FUN = function(x) h(x, t, parameters)[2])
  dim(h_all) <- c(9, 9) 
  h_all[1,1] <- h_all[2,2] <- h_all[3,3] <- h_all[4,4] <- h_all[5,5] <- h_all[6,6] <- h_all[7,7] <- h_all[8,8] <- h_all[9,9] <- 0

  #EQUATIONS
  
  dPA <- - sum(h_all[1,], na.rm = T) * P_A + 
    h_all[2,1] * P_B + h_all[3,1] * P_C + h_all[4,1] * P_D + h_all[5,1] * P_E + 
    h_all[6,1] * P_F + h_all[7,1] * P_G + h_all[8,1] * P_H + h_all[9,1] * P_I
  
  dPB <- - sum(h_all[2,], na.rm = T) * P_B + 
    h_all[1,2] * P_A + h_all[3,2] * P_C + h_all[4,2] * P_D + h_all[5,2] * P_E + 
    h_all[6,2] * P_F + h_all[7,2] * P_G + h_all[8,2] * P_H + h_all[9,2] * P_I
  dPC <- - sum(h_all[3,], na.rm = T) * P_C + 
    h_all[1,3] * P_A + h_all[2,3] * P_B + h_all[4,3] * P_D + h_all[5,3] * P_E + 
    h_all[6,3] * P_F + h_all[7,3] * P_G + h_all[8,3] * P_H + h_all[9,3] * P_I
  dPD <- - sum(h_all[4,], na.rm = T) * P_D + 
    h_all[1,4] * P_A + h_all[2,4] * P_B + h_all[3,4] * P_C + h_all[5,4] * P_E + 
    h_all[6,4] * P_F + h_all[7,4] * P_G + h_all[8,4] * P_H + h_all[9,4] * P_I
  dPE <- - sum(h_all[5,], na.rm = T) * P_E + 
    h_all[1,5] * P_A + h_all[2,5] * P_B + h_all[3,5] * P_C + h_all[4,5] * P_D + 
    h_all[6,5] * P_F + h_all[7,5] * P_G + h_all[8,5] * P_H + h_all[9,5] * P_I
  dPF <- - sum(h_all[6,], na.rm = T) * P_F + 
    h_all[1,6] * P_A + h_all[2,6] * P_B + h_all[3,6] * P_C + h_all[4,6] * P_D + 
    h_all[5,6] * P_E + h_all[7,6] * P_G + h_all[8,6] * P_H + h_all[9,6] * P_I
  dPG <- - sum(h_all[7,], na.rm = T) * P_G + 
    h_all[1,7] * P_A + h_all[2,7] * P_B + h_all[3,7] * P_C + h_all[4,7] * P_D + 
    h_all[5,7] * P_E + h_all[6,7] * P_F + h_all[8,7] * P_H + h_all[9,7] * P_I
  dPH <- - sum(h_all[8,], na.rm = T) * P_H + 
    h_all[1,8] * P_A + h_all[2,8] * P_B + h_all[3,8] * P_C + h_all[4,8] * P_D + 
    h_all[5,8] * P_E + h_all[6,8] * P_F + h_all[7,8] * P_G + h_all[9,8] * P_I
  dPI <- - sum(h_all[9,], na.rm = T) * P_I + 
    h_all[1,9] * P_A + h_all[2,9] * P_B + h_all[3,9] * P_C + h_all[4,9] * P_D + 
    h_all[5,9] * P_E + h_all[6,9] * P_F + h_all[7,9] * P_G + h_all[8,9] * P_H
  
  
  #OUTPUT AS LIST
  list(c(dPA, dPB, dPC, dPD, dPE, dPF, dPG, dPH, dPI))
}



dP_dt_full(t = t, P0=P0, params = parameters)

test <- ode(y = P0, times = t, func = dP_dt_full, parms=parameters)
test
?ode()
sum(c(1,2,3,NA), na.rm = T)
h("AB", t, parameters)[[2]][[1]]

a <- h("AB", t, parameters)[2]


###BACKUP dP_dt_full####
dP_dt_full <- function(t, P0, params) {
  #PREP: define initial conditions and get all get all hij(t)
  P_A <- P0[1]
  P_B <- P0[2]
  P_C <- P0[3]
  P_D <- P0[4]
  P_E <- P0[5]
  P_F <- P0[6]
  P_G <- P0[7]
  P_H <- P0[8]
  P_I <- P0[9]
  
  all_transitions <- paste0(rep(LETTERS[1:9], each=9), rep(LETTERS[1:9], 9))
  h_all <- sapply(all_transitions, FUN = function(x) {ifelse(substring(x, 1, 1) == substring(x, 2, 2), 
                                               return(0), return(as.list(h(x, t, parameters))))})
  print(h_all)
  dim(h_all) <- c(9, 9) 
  
  #EQUATIONS
  
  dPA <- - sum(apply(h_all[1,], FUN= function(x) return(unlist(x)[2])), na.rm = T) * P_A + 
    unlist(h_all[2,1])[2] * P_B + unlist(h[3,1])[2] * P_C + unlist(h[4,1])[2] * P_D + unlist(h[5,1])[2] * P_E + 
    unlist(h_all[6,1])[2] * P_F + unlist(h[7,1])[2] * P_G + unlist(h[8,1])[2] * P_H + unlist(h[9,1])[2] * P_I
  
  dPB <- - sum(apply(h_all[2,], FUN= function(x) return(unlist(x)[2])), na.rm = T) * P_B + 
    unlist(h_all[1,2])[2] * P_B + unlist(h_all[3,2])[2] * P_C + unlist(h_all[4,2])[2] * P_D + unlist(h_all[5,2])[2] * P_E + 
    unlist(h_all[6,2])[2] * P_F + unlist(h_all[7,2])[2] * P_G + unlist(h_all[8,2])[2] * P_H + unlist(h_all[9,2])[2] * P_I
  dPC<- - sum(apply(h_all[3,], FUN= function(x) return(unlist(x)[2])), na.rm = T) * P_C + 
    unlist(h_all[1,3])[2] * P_B + unlist(h_all[2,3])[2] * P_C + unlist(h_all[4,3])[2] * P_D + unlist(h_all[5,3])[2] * P_E + 
    unlist(h_all[6,3])[2] * P_F + unlist(h_all[7,3])[2] * P_G + unlist(h_all[8,3])[2] * P_H + unlist(h_all[9,3])[2] * P_I
  dPD <- - sum(apply(h_all[4,], FUN= function(x) return(unlist(x)[2])), na.rm = T) * P_D + 
    unlist(h_all[1,4])[2] * P_B + unlist(h_all[2,4])[2] * P_C + unlist(h_all[3,4])[2] * P_D + unlist(h_all[5,4])[2] * P_E + 
    unlist(h_all[6,4])[2] * P_F + unlist(h_all[7,4])[2] * P_G + unlist(h_all[8,4])[2] * P_H + unlist(h_all[9,4])[2] * P_I
  dPE <- - sum(apply(h_all[5,], FUN= function(x) return(unlist(x)[2])), na.rm = T) * P_E + 
    unlist(h_all[1,5])[2] * P_B + unlist(h_all[2,5])[2] * P_C + unlist(h_all[3,5])[2] * P_D + unlist(h_all[4,5])[2] * P_E + 
    unlist(h_all[6,5])[2] * P_F + unlist(h_all[7,5])[2] * P_G + unlist(h_all[8,5])[2] * P_H + unlist(h_all[9,5])[2] * P_I
  dPF <- - sum(apply(h_all[6,], FUN= function(x) return(unlist(x)[2])), na.rm = T) * P_F + 
    unlist(h_all[1,6])[2] * P_B + unlist(h_all[2,6])[2] * P_C + unlist(h_all[3,6])[2] * P_D + unlist(h_all[4,6])[2] * P_E + 
    unlist(h_all[5,6])[2] * P_F + unlist(h_all[7,6])[2] * P_G + unlist(h_all[8,6])[2] * P_H + unlist(h_all[9,6])[2] * P_I
  dPG <- - sum(apply(h_all[7,], FUN= function(x) return(unlist(x)[2])), na.rm = T) * P_G + 
    unlist(h_all[1,7])[2] * P_B + unlist(h_all[2,7])[2] * P_C + unlist(h_all[3,7])[2] * P_D + unlist(h_all[4,7])[2] * P_E + 
    unlist(h_all[5,7])[2] * P_F + unlist(h_all[6,7])[2] * P_G + unlist(h_all[8,7])[2] * P_H + unlist(h_all[9,7])[2] * P_I
  dPH <- - sum(apply(h_all[8,], FUN= function(x) return(unlist(x)[2])), na.rm = T) * P_H + 
    unlist(h_all[1,8])[2] * P_B + unlist(h_all[2,8])[2] * P_C + unlist(h_all[3,8])[2] * P_D + unlist(h_all[4,8])[2] * P_E + 
    unlist(h_all[5,8])[2] * P_F + unlist(h_all[6,8])[2] * P_G + unlist(h_all[7,8])[2] * P_H + unlist(h_all[9,8])[2] * P_I
  dPI <- - sum(apply(h_all[9,], FUN= function(x) return(unlist(x)[2])), na.rm = T) * P_I + 
    unlist(h_all[1,9])[2] * P_B + unlist(h_all[2,9])[2] * P_C + unlist(h_all[3,9])[2] * P_D + unlist(h_all[4,9])[2] * P_E + 
    unlist(h_all[5,9])[2] * P_F + unlist(h_all[6,9])[2] * P_G + unlist(h_all[7,9])[2] * P_H + unlist(h_all[8,9])[2] * P_I
  
  
  #OUTPUT AS LIST
  list(c(dPA, dPB, dPC, dPD, dPE, dPF, dPG, dPH, dPI))
}


#####TESTING STUFF#########

##For the combined function:
#prep:
a <- paste0(rep(LETTERS[1:9], each=9), rep(LETTERS[1:9], 9))

b <- sapply(a, FUN = function(x) {ifelse(substring(x, 1, 1) == substring(x, 2, 2), 
                                         return(0), return(as.list(h(x, t, parameters))))})

dim(b) <- c(9,9)
unlist(b[5,5])[2]
b[1,1]

#For the h function
substr("AB", 2,2)
plot(unlist(b[1,9])[2])
unlist(b[1,1])[2]

tn <- 7
ifelse(tn%/%10==0, paste0("0", as.character(tn)), as.character(tn))

rownames(inla.allj.T.C.Pc$summary.fixed)

rep(c("I", "T", "C", "f", "h", "p"), each=1, 72)
rep(1:72, each=6)

paste0(rep(c("I", "T", "C", "f", "h", "p"), each=1, 72), rep(1:72, each=6))

parameters$alpha[names(parameters$alpha) == as.character(h("AB", t, parameters))]

parameters$beta[names(parameters$beta) %in% paste0(c("I", "T", "C", "f", "h", "p"), 12)]

h("AB", t, parameters)



t_n <- 12
test <- exp(parameters$beta[names(parameters$beta) %in% paste0(c("I"), t_n)]) * parameters$alpha[t_n] * t^(parameters$alpha[t_n] - 1)
plot(t, test)

points(t, h("BE", t, parameters), col="red")

#For the dP_dt function
dPA <- -(h("AB") + h("AC") + h("AD") + h("AE") + h("AF") + h("AG") + h("AH") + h("AI")) * P["A"] +
  h("BA") * P["B"] + h("CA") * P["C"] + h("DA") * P["D"] + h("EA") * P["E"] + h("FA") * P["F"] + h("GA") * P["G"] +
  h("HA") * P["H"] + h("IA") * P["I"]

f <- "A"
t_l <- setdiff(LETTERS[1:9], f)
paste0(f, t_l)

leave <- sum(sapply(paste0(f, t_l), h))

sum(sapply(1:8, FUN=print(paste0(t_l, f)[x])))

sum(sapply(1:8, function(x) {h(paste0(t_l, f)[x]) * P(t_l[x])}))




paste0(t_l, f)[3]
return <- 0
  
