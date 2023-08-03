#### Runs the markov chain model using the msm package ####

### Data ###
#Load data
data <- read.csv("../../../Data/BTE/bte_sp_class.csv")[,-1]

# Create dataframe with msm compatible format