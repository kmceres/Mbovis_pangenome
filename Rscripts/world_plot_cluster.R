library(dplyr)
library(stringr)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(geodist)
library(rgeos)
library(sf)

metadata <- read.csv("data/full_SRA_metadata.csv")
metadata <- metadata %>% select(BioSample, geo_loc_name, geo_loc_name_country)
genepa <- read.csv("data/gene_presence_absence_roary.csv")
genomes <- colnames(genepa)[15:ncol(genepa)]
genomes <- str_remove(genomes,"_assembled")
metadata <- metadata %>% filter(BioSample %in% genomes) %>% 
  separate(geo_loc_name, into = c("Country", "State"), sep=":") 
metadata <- metadata %>% separate(State, into = c("loc1", "loc2"), sep=",")
states <- read.csv("data/locations.csv")
states <- states %>% filter(BioSample %in% genomes)

#create clusters based on distance
dist_tib <- states%>% select(latitude, longitude)
dist <- geodist(dist_tib)
colnames(dist) = states$id
rownames(dist) = states$id
Delta <- dist^2/(2*nrow(dist))
Delta <- as.dist(Delta)
tree <- hclust(Delta,method="ward.D")
p10 <- cutree(tree, k=10)
p10_df <- tibble(id = rownames(dist),cluster=p10)
#write out clusters
#write.csv(p10_df,"data/wards_clusters_10.csv")
p10_df <- read.csv("data/wards_clusters_10.csv")
states <- left_join(states, p10_df, by="id")
states <- states %>% group_by(cluster) %>% add_tally(name="n_cluster")
states$geo_loc_name_country <- ifelse(states$geo_loc_name_country == "USA",
                                      "United States of America",states$geo_loc_name_country)

# plot world
world <- ne_countries(scale = "medium", returnclass = "sf")
world_df <- as.data.frame(world)
world_df$geo_loc_name_country = world_df$sovereignt
world_df <- left_join(world_df,states,by="geo_loc_name_country")
world <- st_as_sf(world_df)
states_f <- states %>% distinct(cluster,.keep_all=T)
c<- viridis(10)
#write.table(c,"data/colors_world_plot.txt")
c <- read.table("data/colors_world_plot.txt")
c <- c(c$x)
# plot
world_plot_cluster <- ggplot() +
  geom_sf(data = world) +
  geom_point(data=states_f,aes(x = longitude, y = latitude, size = n_cluster, 
            color=as.factor(cluster)))+
  scale_color_manual(values=c)+
  guides(color = "none")+labs(size="n")+theme_bw()

ggsave("plots/Figure1_world_plot_cluster.svg", world_plot_cluster,width = 5,height=6)




