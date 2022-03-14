## assessing genome quality
library(ggplot2)
library(dplyr)
library(ggpubr)
library(seqinr)
library(tidyr)

## quast results
metadata <- read.csv("data/full_metadata_quast.csv")
#plot results
contigs <- ggplot(metadata)+geom_histogram(aes(x=X..contigs))+labs(x= "Number of contigs")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
N50 <- ggplot(metadata)+geom_histogram(aes(x=N50))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
NG50 <- ggplot(metadata)+geom_histogram(aes(x=NG50)) 

mean(metadata$X..contigs)
median(metadata$X..contigs)
range(metadata$X..contigs)
median(metadata$N50)
range(metadata$N50)

## checkM results
checkM_p <- ggplot(metadata) + geom_point(aes(x=Contamination,y=Completeness)) + 
  labs(x="Contamination %", y="Completeness %")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

## save quality control plot
q_out <- ggarrange(contigs,N50,checkM_p,nrow=1, labels = c("A","B","C"))
ggsave("~/PhD/complete_genome_set_bovis_pangenome/plots/quality_results.png",width=8,height=3)


