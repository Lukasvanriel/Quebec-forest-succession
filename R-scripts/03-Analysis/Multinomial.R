
library(nnet)

#TODO: review and update

##Function that performs a multinomial calculation on a dataset
#data: Must be in a shape that is supported by the msm package. 
#most importantly it must contain a column "sp_class"
multinom_model <- function(data) {
  data_multinom <- list()
  model <- list()
  probabilities <- list()
  
  for (i in 1:max(data$sp_class)) { 
    data_place <- data |> 
      filter(sp_class == i, time == 0)  |>
      select(TESSELLE)
    
    data_multinom[[i]] <- subset(data, TESSELLE %in% data_place$TESSELLE)
    if(! i %in% data_multinom[[i]]$sp_class){
      data_multinom[[i]]$sp_factor <- as.factor(data_multinom[[i]]$sp_class)
      data_multinom[[i]]$sp_factor <- factor(data_multinom[[i]]$sp_factor,
                                             levels = c(levels(data_multinom[[i]]$sp_factor), as.character(i)))
      data_multinom[[i]]$sp_factor <- relevel(data_multinom[[i]]$sp_factor,
                                              ref = as.character(i))
    } else {
    data_multinom[[i]]$sp_factor <- relevel(as.factor(data_multinom[[i]]$sp_class),
                                            ref = as.character(i))
    }
    
    #There seems to be a lot of missing Jackpine, so I artificially had to fix some of these transitions
    if(i!=8 || length(table(data_multinom[[i]]$sp_factor)) > 1){
    model[[i]] <- multinom(sp_factor ~ time, data = data_multinom[[i]])
    c <- coef(model[[i]])
    missing.prob <- c()
    names <- as.character(seq(1,max(data$sp_class)))
    if(nrow(c) < (max(data$sp_class)-1)) {
      ref <- names[-i]
      missing.coef <- ref[! ref %in% rownames(c)]
      missing.prob <- rep(0, length(missing.coef))
      attributes(missing.prob)$names <- missing.coef
    }
    
    exp.coef <- exp(rowSums(c))
    prob <- exp.coef / (1 + sum(exp.coef))
    probabilities[[i]] <- c(1-sum(prob), prob, missing.prob)
    attributes(probabilities[[i]])$names[1] <- as.character(i)
    } else{ 
      probabilities[[i]] <- c(0.914, 0.01, 0.007, 0.000, 0.000, 0.01, 0.01, 0.036, 0.013)
      attributes(probabilities[[i]])$names <- as.list(c("8", "1", "2", "3", "4", "5", "6", "7", "9"))
      }
    probabilities[[i]] <- probabilities[[i]][names]
  }
  as.data.frame(do.call(rbind, probabilities))
}

######

data <- read_csv(here("Data", "BTE", paste0("bte_", "4cM", "_msm_ready.csv")),
                 col_types = cols(.default = col_guess(),
                                  sp_class = col_integer(),
                                  cov_pert_class = col_factor(),
                                  cov_pert_sev = col_factor(),
                                  cov_time_pert = col_double()))[,-1]# %>% select(-"...X")

# Rescale covariates where necessary. Add column classes
data_sc <- data
data_sc$cov_CMI[is.na(data_sc$cov_CMI)] <- mean(data_sc$cov_CMI, na.rm =T)
data_sc$cov_Tmean[is.na(data_sc$cov_Tmean)] <- mean(data_sc$cov_Tmean, na.rm =T)
data_sc <- data_sc %>% 
  mutate(cov_CMI = (cov_CMI - mean(cov_CMI)) / sd(cov_CMI)) %>% 
  mutate(cov_Tmean = (cov_Tmean - mean(cov_Tmean)) / sd(cov_Tmean)) %>% 
  mutate(cov_time_pert = cov_time_pert / 100) %>% 
  mutate(cov_soil = cov_soil / 10)

#

mnom <- multinom_model(data_sc)
round(mnom,3)
sum(mnom[8,])

#

data_multinom <- list()
model <- list()
probabilities <- list()

i=8

data_place <- data_sc |> 
  filter(sp_class == i, time == 0)  |>
  select(TESSELLE)

data_multinom[[i]] <- subset(data_sc, TESSELLE %in% data_place$TESSELLE)
data_multinom[[i]]$sp_factor <- relevel(as.factor(data_multinom[[i]]$sp_class),
                                        ref = as.character(i))

table(data_multinom[[i]]$sp_class)




a <- mnom

data_sc$sp_class



