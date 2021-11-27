## assessing genome quality
library(ggplot2)
library(dplyr)
library(ggpubr)
library(seqinr)
library(tidyr)


## quast results
quast <- read.csv("data/quast_cows_150.csv")
#plot results
contigs <- ggplot(quast)+geom_histogram(aes(x=X..contigs))+labs(x= "Number of contigs")
N50 <- ggplot(quast)+geom_histogram(aes(x=N50))
NG50 <- ggplot(quast)+geom_histogram(aes(x=NG50)) 

mean(quast$X..contigs)
median(quast$X..contigs)
range(quast$X..contigs)
median(quast$N50)
range(quast$N50)

## checkM results
checkM <- read.csv("data/checkm_out.csv")
checkM_p <- ggplot(checkM) + geom_point(aes(x=Contamination,y=Completeness)) + 
  labs(x="Contamination %", y="Completeness %")

## save quality control plot
q_out <- ggarrange(contigs,N50,checkM_p,nrow=1, labels = c("A","B","C"))
ggsave("plots/quality_results.png",width=8,height=3)


