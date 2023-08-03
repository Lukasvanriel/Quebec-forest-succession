
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
    data_multinom[[i]]$sp_factor <- relevel(as.factor(data_multinom[[i]]$sp_class),
                                            ref = as.character(i))
    
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
    probabilities[[i]] <- probabilities[[i]][names]
  }
  as.data.frame(do.call(rbind, probabilities))
}
