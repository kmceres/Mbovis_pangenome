library(sf)
library(raster)
library(dplyr)
library(spData)
library(spDataLarge)
library(tmap)    
library(leaflet)
library(ggplot2) 
library(RColorBrewer)
library(ggpubr)

meta <- read.csv("data/full_metadata_quast.csv")

meta <- meta %>% mutate(Country_corrected = ifelse(
  Country_corrected == "United States of America","United States",Country_corrected
))
meta <- meta %>% distinct(Run,.keep_all=T)
meta$name_long <- meta$Country_corrected
meta <- meta %>% ungroup %>% group_by(Country_corrected) %>% add_tally(name="n_country")
world_meta <- left_join(world, meta, by ="name_long")


facets = "n_country"

world_plot <-  tm_shape(world_meta) + 
  tm_polygons(facets,
              pal=c("#481A6CFF", "#414487FF", "#31688EFF", "#23888EFF" ,
                    "#22A884FF", "#54C568FF", "#A5DB36FF" ,"#DCE318FF"),
                    style="log10")+
  tm_facets(nrow = 1, sync = TRUE)

tmap_save(world_plot, "world_plot.pdf")
