#### Script to group species classes ####

bte1_noP <- read.csv("../../../Data/BTE/bte1_noP.csv")[,-1]
bte2_noP <- read.csv("../../../Data/BTE/bte2_noP.csv")[,-1]
bte3_noP <- read.csv("../../../Data/BTE/bte3_noP.csv")[,-1]
bte4_noP <- read.csv("../../../Data/BTE/bte4_noP.csv")[,-1]
bte5_noP <- read.csv("../../../Data/BTE/bte5_noP.csv")[,-1]

### Only keep columns that are present in bte5 ###
bte1_noP <- bte1_noP %>%
  select(intersect(colnames(bte1_noP) , colnames(bte5_noP)))
bte2_noP <- bte2_noP %>%
  select(intersect(colnames(bte2_noP) , colnames(bte5_noP)))
bte3_noP <- bte3_noP %>%
  select(intersect(colnames(bte3_noP) , colnames(bte5_noP)))
bte4_noP <- bte4_noP %>%
  select(intersect(colnames(bte4_noP) , colnames(bte5_noP)))

### Combine data ####
bte_all <- rbind(bte1_noP, bte2_noP, bte3_noP, bte4_noP, bte5_noP)

###There is a significant fraction of TESSELLE that are lacking Species information
#TODO: Look into which to drop and which to keep
bte <- bte_all |> filter(! is.na(GR_ESS))

###Extract the most prevalent classes
most_prev <- bte |> 
  count(GR_ESS) |> arrange(desc(n)) |> head(100)

#Also considering only the dominant species
bte <- bte |>
  mutate(dom_sp = substr(GR_ESS,1, 2))

table(bte$dom_sp)
