#### Define forest classes according to dominant species ####

### Function that returns an integer vector dependend on provided classification (dataframe)

#sp_list: vector of the species at different polygons
#classes: dataframe with 2 columns: 
#   cl: vector of the unique species found in sp_list
#   msm: vector of integers to be associated with species in cl
determine_sp_class <- function(sp_list, classes) {
  int_class <- integer(length(sp_list))
  for (i in 1:length(sp_list)) {
    if (sp_list[i] %in% classes$cl) {
      int_class[i] <- classes$msm[classes$cl == sp_list[i]]
    }
  }
  int_class[int_class == 0] <- max(int_class) + 1
  int_class
}

### Create necessary dataframes

#1: Paper birch
#2: Other shade intolerant species
#3: Yellow birch
#4: Sugar and red maple
#5: Other deciduous
#6: Balsam fir
#7: Red and black spruce
#8: Jack pine
#9: Other coniferous

class.list <- list()
class.list[[1]] <- c("BB", 1)
class.list[[2]] <- c("FI", 2)
class.list[[3]] <- c("PE", 2)
class.list[[4]] <- c("BJ", 3)
class.list[[5]] <- c("SB", 6)
class.list[[6]] <- c("BP", 1)
class.list[[7]] <- c("ER", 4)
class.list[[8]] <- c("EE", 7)
class.list[[9]] <- c("EN", 7)
class.list[[10]] <- c("EP", 7)
class.list[[11]] <- c("ES", 4)
class.list[[12]] <- c("PG", 8)
class.list[[13]] <- c("S", 6)
class.list[[14]] <- c("E", 7)
class.list[[15]] <- c("RX", 9)
class.list[[16]] <- c("FX", 5)
class.list[[17]] <- c("SS", 6)
class.list[[18]] <- c("EO", 4)
class.list[[19]] <- c("SE", 6)
class.list[[20]] <- c("TO", 5)
class.list[[21]] <- c("FT", 5)
class.list[[22]] <- c("RB", 3)
class.list[[23]] <- c("EB", 9)

class.list[[24]] <- c("SC", 6)
class.list[[25]] <- c("PB", 9)
class.list[[26]] <- c("FN", 5)
class.list[[27]] <- c("CS", 5)
class.list[[28]] <- c("RF", 2)
class.list[[29]] <- c("ML", 9)
class.list[[30]] <- c("EC", 7) #?
class.list[[31]] <- c("SP", 2)
class.list[[32]] <- c("RE", 7)
class.list[[33]] <- c("CE", 5)
class.list[[34]] <- c("RP", 9)
class.list[[35]] <- c("CC", 5)
class.list[[36]] <- c("SF", 2) #x#
class.list[[37]] <- c("ME", 9)
class.list[[38]] <- c("EM", 7)
class.list[[39]] <- c("MR", 9)
class.list[[40]] <- c("MF", 2)
class.list[[41]] <- c("EF", 5)
class.list[[42]] <- c("PR", 9)
class.list[[43]] <- c("PA", 2)

class.list[[44]] <- c("RS", 9)
class.list[[45]] <- c("RC", 9)
class.list[[46]] <- c("CB", 5)
class.list[[47]] <- c("PI", 9)
class.list[[48]] <- c("FH", 5)
class.list[[49]] <- c("PT", 2)
class.list[[50]] <- c("CM", 5)
class.list[[51]] <- c("EU", 7)
class.list[[52]] <- c("SM", 6)
class.list[[53]] <- c("RM", 9)
class.list[[54]] <- c("CR", 5)
class.list[[55]] <- c("HG", 5)
class.list[[56]] <- c("CP", 5)
class.list[[57]] <- c("FO", 9)
class.list[[58]] <- c("RZ", 9)
class.list[[59]] <- c("GG", 9)
class.list[[60]] <- c("RG", 9)

msm.class.normal <- as.data.frame(do.call(rbind, class.list))
colnames(msm.class.normal) <- c("cl", "msm")
msm.class.normal$msm <- as.integer(msm.class.normal$msm)
