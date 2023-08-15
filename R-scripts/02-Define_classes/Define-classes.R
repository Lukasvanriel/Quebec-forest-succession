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
bte_filt <- bte |>
  filter(! dom_sp %in% 
           c("TR", "MS", "C", "M", "R", "FS", "MH", "SR", "MX",
             "PU", "RR", "P", "BG", "PS", "F", "EV"))

#### Divide into classes as defined by Grouping.R ####

# Create the numbered class vector
#Load the script that produces the necessary dataframe
source(here("R-scripts", "02-Define_classes", "Grouping.R"))

bte_filt$sp_class <- determine_sp_class(bte_filt$dom_sp, msm.class.normal)

table(bte_filt$sp_class)
bte_filt[bte_filt$sp_class == 10,]$GR_ESS


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
