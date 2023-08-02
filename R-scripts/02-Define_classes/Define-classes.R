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

#filter out species that have weird or no documentation
bte_filt <- bte |>
  filter(! dom_sp %in% 
           c("TR", "MS", "C", "M", "R", "FS", "MH", "SR", "MX",
             "PU", "RR", "P", "BG", "PS", "F"))

# Create the numbered class vector
bte_filt$sp_class <- determine_sp_class(bte_filt$dom_sp, msm.class.normal)

#table(bte_filt$dom_sp[test==10])

write.csv(bte_filt , "../../../Data/BTE/bte_sp_class.csv")
