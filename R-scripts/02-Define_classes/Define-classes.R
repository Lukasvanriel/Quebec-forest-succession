#### Script to group species classes ####

#### Load data ####
bte1_noP <- read.csv("../../../Data/BTE/bte1_noP.csv")[,-1]
bte2_noP <- read.csv("../../../Data/BTE/bte2_noP.csv")[,-1]
bte3_noP <- read.csv("../../../Data/BTE/bte3_noP.csv")[,-1]
bte4_noP <- read.csv("../../../Data/BTE/bte4_noP.csv")[,-1]
bte5_noP <- read.csv("../../../Data/BTE/bte5_noP.csv")[,-1]

bte <- read.csv("../../../Data/BTE/bte_all.csv")[,-1]

#### Extract the most prevalent classes ####
most_prev <- bte |> 
  count(GR_ESS) |> arrange(desc(n)) |> head(100)

#Also considering only the dominant species
bte <- bte |>
  mutate(dom_sp = substr(GR_ESS,1, 2))

most_prev_dom <- bte |> 
  count(dom_sp) |> arrange(desc(n)) |> head(100)

#### Divide into classes as defined by Grouping.R ####

