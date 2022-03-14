source("usda_snp_tab_to_vcf.R")
source("snps_to_boolean_pca.R")
library(ape)
library(ggplot2)
library(dplyr)
library(viridis)
library(seqinr)
library(factoextra)
library(abdiv)
library(ggpubr)
library(svglite)

### SNPS #####
snps<- seqinr::read.fasta("data/core_genome.filtered_polymorphic_sites.fasta")
snps2 <- snps[names(snps)!="h37rv.gb.fna"]
snps_table <- tibble(id = rep(NA, 1463), dna=rep(NA,1463))

for(s in 1:length(snps2)){
  seq = snps2[[s]]
  snps_table$id[s] = attr(seq, "name")
  snps_table$dna[s] = paste(unlist(as.character(toupper(seq[1:length(seq)]))),collapse="")
}

write.table(snps_table, "data/concat_manual.filtered.snps.txt",sep="\t", row.names = F, col.names = F, quote = F)
meta <- read.csv("data/full_metadata_quast.csv")
lineages <- meta %>% ungroup() %>% 
  dplyr::select(Bin,L2_corrected,Continent,Country_corrected) %>%
  mutate(BioSample = Bin) %>% select(-Bin)
#pca 
b <- snps_to_boolean_vector("data/concat_manual.filtered.snps.txt")
d <- prepare_pca(b)

pca <- prcomp(d, center = F, scale. = F) # prepare_pca already centers and scales

#eigenvalue screeplot: Show the percentage of variances explained by each principal component.
var <- fviz_eig(pca)
eig <- get_eig(pca)
#plot pcs 
pcs <- as.data.frame(pca$x)
pcsdf <- pcs
pcsdf$BioSample <- str_remove(rownames(pcs),".fasta")
## kmeans cluster pcs
opt_n_clusters <- fviz_nbclust(pcs,kmeans, "silhouette",k.max=20) # get optimal number of clusters = 10
seed=19
k <- kmeans(pcs, 10, nstart=25, iter.max=1000)
km_df <- k$cluster
km_df <- tibble(BioSample = str_remove(names(km_df),".fasta"),
                cluster=km_df)
colors <- c(viridis(5,option="plasma"),viridis(5,option="D"))

colors <- gplots::col2hex(colors)
color_df <- tibble(color = colors, cluster=1:10)
km_color <- left_join(km_df, color_df) %>% relocate(BioSample, color, cluster)
km_color <- left_join(km_color, lineages)

pcsdf <- left_join(pcsdf, km_color)
order <- pcsdf %>% arrange(cluster)



lineage_pal = c("#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#E6F598", "#ABDDA4" ,"#66C2A5", "#3288BD")
top_eig <- eig[1:10,]
top_eig$Dimension <- c(1:10)
var_plot <- ggplot(top_eig)+geom_col(aes(x=as.factor(Dimension),y=variance.percent))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="Dimensions", y="Percentage of explained variance")


 # scale_color_viridis(discrete=T,option="plasma")+theme_bw()+labs(color="Kmeans cluster")
ps.1.2_l<-ggplot(pcsdf) + geom_point(aes(x=PC1,y=PC2,col=as.factor(L2_corrected))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values=lineage_pal)+
  labs(x="PC1 17.8%", y="PC2 13.9%")

ps.1.3_l<-ggplot(pcsdf) + geom_point(aes(x=PC1,y=PC3,col=as.factor(L2_corrected))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values=lineage_pal)+
  labs(x="PC1 17.8%", y="PC3 8.32")

ps.1.4_l<-ggplot(pcsdf) + geom_point(aes(x=PC1,y=PC4,col=as.factor(L2_corrected)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values=lineage_pal)+
  labs(x="PC1 17.8%", y="PC4 6.01%")

ps.1.5_l<-ggplot(pcsdf) + geom_point(aes(x=PC1,y=PC5,col=as.factor(L2_corrected)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values=lineage_pal)+
  labs(x="PC1 17.8%", y="PC5 5.92%")

p1<-ggarrange(ps.1.2_l,ps.1.3_l,ps.1.4_l,ps.1.5_l,common.legend = T,nrow=1)

# lineage_cluster <- ggplot(km_color)+geom_bar(aes(x=as.factor(cluster),fill=L2_corrected))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
#   scale_fill_manual(values=lineage_pal)+
#   labs(x="Cluster",y="Number of genomes",fill="Lineage")

### acc genes #########
## pca on gene presence absence matrix 
genepa <- read.csv("data/gene_presence_absence_roary_01.csv")
genes <- genepa$Gene
genepa<- genepa[,16:ncol(genepa)]
jaccard_dist <- matrix(NA, ncol(genepa), ncol(genepa))
rownames(jaccard_dist) = colnames(jaccard_dist) = colnames(genepa)
for (c1 in 1:ncol(genepa)){
  for (c2 in 1:ncol(genepa)){
    jaccard_dist[c1,c2] <- jaccard(genepa[,c1], genepa[,c2])
  }
}

pca_genepa <- prcomp(jaccard_dist, center = T, scale. = T) 
pca_genepa_df <- as.data.frame(pca_genepa$x) 
pca_genepa_df$BioSample <- rownames(pca_genepa_df)

#plot pcs
#eigenvalue screeplot: Show the percentage of variances explained by each principal component.
fviz_eig(pca_genepa)
eig_acc <- get_eig(pca_genepa)
pcs_acc <- left_join(pca_genepa_df, km_color)
ps.2.1.2_l<-ggplot(pcs_acc) + geom_point(aes(x=PC1,y=PC2,col=as.factor(L2_corrected))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values=lineage_pal)+
  labs(x="PC1 50.9%", y="PC2 25.0%")

ps.2.1.3_l<-ggplot(pcs_acc) + geom_point(aes(x=PC1,y=PC3,col=as.factor(L2_corrected))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values=lineage_pal)+
  labs(x="PC1 50.9%", y="PC3 8.24%")

ps.2.1.4_l<-ggplot(pcs_acc) + geom_point(aes(x=PC1,y=PC4,col=as.factor(L2_corrected))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values=lineage_pal)+
  labs(x="PC1 50.9%", y="PC4 5.61%")

ps.2.1.5_l<-ggplot(pcs_acc) + geom_point(aes(x=PC1,y=PC5,col=as.factor(L2_corrected))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values=lineage_pal)+
  labs(x="PC1 50.9%", y="PC5 3.57%")

p2<-ggarrange(ps.2.1.2_l,ps.2.1.3_l,ps.2.1.4_l,ps.2.1.5_l,common.legend = T,nrow=1)
pca_plot <- ggarrange(p1,p2,ncol=1,common.legend=T)

hc <- hclust(as.dist(jaccard_dist), method = "ward.D2")
tree<-as.phylo(hc)

plot(tree, cex = 0.6, label.offset = 0.5)
write.tree(tree,"tree_files/unfiltered_acc_gene_dendrogram.nwk")

## filtered acc genes
filt <- read.table("data/filtered_accessory_genes.txt")
jaccard_dist_filt <- matrix(NA, nrow(filt), nrow(filt))
rownames(jaccard_dist_filt) = colnames(jaccard_dist_filt) = rownames(filt)
for (r1 in 1:nrow(filt)){
  for (r2 in 1:nrow(filt)){
    jaccard_dist_filt[r1,r2] <- jaccard(filt[r1,], filt[r2,])
  }
}

genepa_01 <- genepa %>% mutate(Gene = genes)
phage_l <- c("group_1979", "group_616","group_2566", "group_2533", "group_2405",
             "group_2345", "group_601",  "group_373",  "group_2670","group_1978")

km_color_a <-left_join(km_color, phage_genes) 

pca_genepa_filt <- prcomp(jaccard_dist_filt, center = T, scale. = T) 
pca_genepa_filt_df <- as.data.frame(pca_genepa_filt$x) 
pca_genepa_filt_df$BioSample <- rownames(pca_genepa_filt_df)

fviz_eig(pca_genepa_filt)
eig_filt <- get_eig(pca_genepa_filt)
eig_acc_filt <- get_eig(pca_genepa_filt)
pcs_acc_filt <- left_join(pca_genepa_filt_df, km_color_a)

filtered_lin_12_f<-ggplot(pcs_acc_filt) + geom_point(aes(x=PC1,y=PC2,col=as.factor(L2_corrected))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_color_manual(values=lineage_pal)+
  labs(x="PC1 72.3", y="PC2 15.5")


phage <- ggplot(pcs_acc_filt) + geom_point(aes(x=PC1,y=PC2,col=as.factor(Phage_present))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(x="PC1 72.3", y="PC2 15.5")

