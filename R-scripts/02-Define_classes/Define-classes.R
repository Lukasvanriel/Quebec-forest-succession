#### Script to group species classes ####

#### Load data ####
bte <- read.csv(here("Data", "BTE", "bte_all_cov.csv"))[,-1]

#### Extract the most prevalent classes ####
most_prev <- bte |> 
  count(GR_ESS) |> arrange(desc(n)) |> head(100)

#Also considering only the dominant species
bte <- bte |>
  mutate(dom_sp = substr(GR_ESS,1, 2))

most_prev_dom <- bte |> 
  count(dom_sp) |> arrange(desc(n)) |> head(100)

#filter out species that have weird or no documentation
#TODO: or rename
docu <- read.csv(here("Data", "Documentation", "Groupement-essences.csv")) |>
  rename(GR_ESS = Code)

docu_filt <- docu
docu_filt$Description <- sub('\\(BFEC 3-4\\) ', '', docu_filt$Description)
docu_filt <- docu_filt[! duplicated(docu_filt[,c("GR_ESS", "Description")]),]


table(bte$GR_ESS[bte$dom_sp %in% c("TR", "MS", "C", "M", "R", "FS", "MH", "SR", "MX",
                                  "RR", "P", "PS", "F", "EV")])

for(str in c("TR", "MS", "C", "M", "R", "FS", "MH", "SR", "MX",
             "RR", "P", "PS", "F", "EV")) {
  print(sapply(names(table(bte$GR_ESS[bte$dom_sp == str])),
         FUN=function(x) docu$Description[docu$GR_ESS == x]))
}

#Filter out the remaining:
bte_filt <- bte |>
  filter(! dom_sp %in% 
           c("TR", "MS", "C", "M", "R", "FS", "MH", "SR", "MX",
             "RR", "P", "PS", "F", "EV"))

#### Divide into classes as defined by Grouping.R ####

# Create the numbered class vector
#Load the script that produces the necessary dataframe
source(here("R-scripts", "02-Define_classes", "Grouping.R"))

bte_filt$sp_class <- determine_sp_class(bte_filt$dom_sp, msm.class.normal)

table(bte_filt$sp_class)
bte_filt[bte_filt$sp_class == 10,]$GR_ESS

### Check out the different classes ###
for(i in 1:max(bte_filt$sp_class)) {
  sp_gr <- bte_filt %>% 
    filter(sp_class==i) %>% 
    group_by(GR_ESS) %>% 
    summarise(n=n()) %>% 
    arrange(desc(n)) %>% 
    left_join(., select(docu_filt, GR_ESS, Description), by="GR_ESS")
  
  write.csv(sp_gr , here("Data", "Documentation", paste0("Docu_Group", as.character(i),
                                                         ".csv")))
}

## Correct some of the mistakes:
#TODO: Find way to check the classes without documentatin (NA)

bte_filt[bte_filt$GR_ESS == "RFH",]$sp_class <- 9
bte_filt[bte_filt$GR_ESS == "SPE",]$sp_class <- 2
bte_filt[bte_filt$GR_ESS == "RBB",]$sp_class <- 9
bte_filt[bte_filt$GR_ESS == "RBJ-",]$sp_class <- 9
bte_filt[bte_filt$GR_ESS == "RBJ+",]$sp_class <- 9
bte_filt[bte_filt$GR_ESS == "EBB",]$sp_class <- 7
  
# Write out results 
write.csv(bte_filt , here("Data", "BTE", "bte_cov_class.csv"))

### Inspect classes ###

plot.sifort.classes <- function(data, inv) {
  data.inv <- data |>
    filter(NO_PRG == inv) |>
    mutate(sp_class=factor(sp_class))
  
  my_colours <- c("black", "darkblue", "yellow", "red", "orange", "green", "darkgreen", "gray", "cyan")
  ggplot(data.inv) +
    geom_point(aes(x=LONGI, y=LATIT, col=sp_class), size=0.1) +
    scale_color_manual(values = my_colours) +
    theme_bw()
}

pdf("/Users/lukas/Desktop/classes1.pdf")
plot.sifort.classes(data, 1)
dev.off()
pdf("/Users/lukas/Desktop/classes2.pdf")
plot.sifort.classes(data, 2)
dev.off()
pdf("/Users/lukas/Desktop/classes3.pdf")
plot.sifort.classes(data, 3)
dev.off()
pdf("/Users/lukas/Desktop/classes4.pdf")
plot.sifort.classes(data, 4)
dev.off()
pdf("/Users/lukas/Desktop/classes5.pdf")
plot.sifort.classes(data, 5)
dev.off()



########### 
bte_filt <- read.csv(here("Data", "BTE", "bte_cov_class.csv"))
docu <- read.csv(here("Data", "Documentation", "Groupement-essences.csv")) |>
  rename(GR_ESS = Code)
docu_sub <- docu[c(4, 13, 39987, 3532,1918,39719, 39735, 1909,3515),]

#first get rid of the (BFEC 3-4) 
docu_sub_filt <- docu_sub
docu_sub_filt$Description <- sub('\\(BFEC 3-4\\) ', '', docu_sub_filt$Description)
docu_sub_filt <- docu_sub_filt[! duplicated(docu_sub_filt[,c("GR_ESS", "Description")]),]

docu[docu$GR_ESS == "FIPG",]
docu[docu$GR_ESS == "BPBP",][1,3] == docu[docu$GR_ESS == "BPBP",][2,"Description"]
duplicated(docu[docu$GR_ESS == "BPBP",c("GR_ESS", "Description")])

duplicated(docu_sub[,c("GR_ESS", "Description")])

docu[docu$GR_ESS == "PISE",]

table(bte_filt[bte_filt$GR_ESS == "PISE",]$NO_PRG)



