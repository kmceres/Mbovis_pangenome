source("Rscripts/usda_snp_tab_to_vcf.R")
source("Rscripts/snps_to_boolean_pca.R")
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
snps<- seqinr::read.fasta("data/concat_gblocks_aln.filtered_polymorphic_sites.fasta")
snps_table <- tibble(id = rep(NA, 851), dna=rep(NA,851))

for(s in 1:length(snps)){
  seq = snps[[s]]
  snps_table$id[s] = attr(seq, "name")
  snps_table$dna[s] = paste(unlist(as.character(toupper(seq[1:length(seq)]))),collapse="")
}

write.table(snps_table, "data/concat_manual.filtered.snps.txt",sep="\t", row.names = F, col.names = F, quote = F)


#pca 
b <- snps_to_boolean_vector("data/concat_manual.filtered.snps.txt")
d <- prepare_pca(b)

pca <- prcomp(d, center = F, scale. = F) # prepare_pca already centers and scales

#eigenvalue screeplot: Show the percentage of variances explained by each principal component.
fviz_eig(pca)

#plot pcs 
pcs <- as.data.frame(pca$x)
pcsdf <- pcs
pcsdf$BioSample <- rownames(pcs)

## kmeans cluster pcs
fviz_nbclust(pcs,kmeans) # get optimal numebr of clusters = 9
k <- kmeans(pcs, 9, nstart=25, iter.max=1000)
km_df <- k$cluster
km_df <- tibble(BioSample = names(km_df),
                cluster=km_df)
colors <- viridis_pal(option="plasma")(9)
colors <- gplots::col2hex(colors)
color_df <- tibble(color = colors, cluster=1:9)
km_color <- left_join(km_df, color_df) %>% relocate(BioSample, color, cluster)
write.table(km_color,"data/km_pca.csv",sep=",",quote=F,row.names = F)
km_color <- read.csv("data/km_pca.csv")
pcsdf <- left_join(pcsdf, km_color)

p1.1<-ggplot(pcsdf) + geom_point(aes(x=PC1,y=PC2,col=as.factor(cluster)))+ 
  scale_color_viridis(discrete=T,option="plasma")+theme_bw()+labs(color="Kmeans cluster")
p1.2<-ggplot(pcsdf) + geom_point(aes(x=PC1,y=PC3,col=as.factor(cluster)))+ 
  scale_color_viridis(discrete=T,option="plasma")+theme_bw()+labs(color="Kmeans cluster")
p1.3<-ggplot(pcsdf) + geom_point(aes(x=PC1,y=PC4,col=as.factor(cluster)))+ 
  scale_color_viridis(discrete=T,option="plasma")+theme_bw()+labs(color="Kmeans cluster")

### acc genes #########
## pca on gene presence absence matrix 
genepa <- read.csv("data/gene_presence_absence_roary_10.csv")
genes <- genepa$Gene
genepa<- genepa[,16:ncol(genepa)]
jaccard_dist <- matrix(NA, ncol(genepa), ncol(genepa))
rownames(jaccard_dist) = colnames(jaccard_dist) = colnames(genepa)
for (c1 in 1:ncol(genepa)){
  for (c2 in 1:ncol(genepa)){
    jaccard_dist[c1,c2] <- jaccard(genepa[,c1], genepa[,c2])
  }
}
write.table(jaccard_dist,"acc_genes_jaccard_dist.txt", sep=" ", quote=F)
jaccard_dist <- read.table("data/acc_genes_jaccard_dist.txt")

pca_genepa <- prcomp(jaccard_dist, center = T, scale. = T) 
pca_genepa_df <- as.data.frame(pca_genepa$x) 
pca_genepa_df$BioSample <- rownames(pca_genepa_df)

#plot pcs
#eigenvalue screeplot: Show the percentage of variances explained by each principal component.
fviz_eig(pca_genepa)
pcs_acc <- left_join(pca_genepa_df, km_color)
p2.1<- ggplot(pcs_acc) + geom_point(aes(x=PC1,y=PC2,col=as.factor(cluster)))+ 
  scale_color_viridis(discrete=T,option="plasma")+theme_bw()+labs(color="Kmeans cluster")
p2.2 <-ggplot(pcs_acc) + geom_point(aes(x=PC1,y=PC3,col=as.factor(cluster)))+ 
  scale_color_viridis(discrete=T,option="plasma")+theme_bw()+labs(color="Kmeans cluster")
p2.3<-ggplot(pcs_acc) + geom_point(aes(x=PC2,y=PC4,col=as.factor(cluster)))+ 
  scale_color_viridis(discrete=T,option="plasma")+theme_bw()+labs(color="Kmeans cluster")

# combine plots 
pca_plot <- ggarrange(p1.1,p1.2,p2.1,p2.2,nrow=2,ncol=2,common.legend = T)
ggsave("plots/Figure2B.svg",pca_plot)

hc <- hclust(as.dist(jaccard_dist), method = "ward.D2")
tree<-as.phylo(hc)

plot(tree, cex = 0.6, label.offset = 0.5)
write.tree(tree,"data/iTOL/unfiltered_acc_gene_dendrogram.nwk")
