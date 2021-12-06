library(dplyr)
library(stringr)
library(tidyr)
library(Biostrings)

ref_matches <- read.csv("data/ref_matches.csv")
outbreaks <- read.csv("data/outbreak_color_label.csv")
filt_genes <- read.csv("~/PhD/Mbovis_pangenome/data/filtered_genes.csv")
phage_l <- c( "group_1898", "group_1296", "group_1215", "group_1173", 
              "group_1033", "group_979", "group_955",  "group_902",  
              "group_742",  "group_350",  "group_28",   "group_16", "group_6" )
filt_genes <- filt_genes %>% filter(No..isolates < 851)#.95*851)
genes <- read.csv("data/gene_presence_absence_roary_10.csv")
genes <- genes %>% filter(Gene %in% filt_genes$Gene)
outbreaks <- outbreaks %>% 
  filter(Biosample %in% str_remove(colnames(genes)[16:ncol(genes)],"_assembled"))

long <- genes %>% select(-colnames(genes)[3:15]) %>% select(-X) %>%
  gather("Biosample", "Present", -Gene)
long$Biosample <- str_remove(long$Biosample, "_assembled")

long_outbreak <- long %>% filter(Biosample %in% outbreaks$Biosample)
long_outbreak <- left_join(long_outbreak, outbreaks)
long_outbreak <- long_outbreak %>%  mutate(group_n = as.numeric(as.factor(Outbreak))) %>% 
  group_by(Gene) %>% mutate(n_gene = sum(Present)) %>% filter(n_gene > 0) %>%
  filter(n_gene < nrow(outbreaks))
ref_matches_f <- ref_matches %>% mutate(Gene = first.names) %>% select(Gene, Name)
long_outbreak <- left_join(long_outbreak, ref_matches_f)
long_outbreak <- long_outbreak %>% 
  mutate(Name = ifelse(Gene %in% phage_l, "PhiRv1",Name))  %>% 
  filter(is.na(Name)==F)

order <- long_outbreak %>% arrange(Name) %>% distinct(Gene, .keep_all=T)
long_outbreak$Gene <- factor(long_outbreak$Gene, levels=order$Gene)
order2 <- long_outbreak %>% arrange(Outbreak) %>% ungroup() %>% distinct(Biosample, .keep_all=T)
long_outbreak$Biosample <- factor(long_outbreak$Biosample, levels=order2$Biosample)
long_outbreak <- long_outbreak %>% 
  mutate(has_PE_PPE = ifelse((str_detect(Name, "PPE")==T | 
                                str_detect(Name, "PE")==T), 1,0))

outbreak_w_PE_PPE <- ggplot(long_outbreak) + 
  geom_tile(aes(x=Gene, y=Biosample, fill=as.factor(Present),color=Name, group=Outbreak))+
  scale_fill_manual(values=c("#EEEEEE","#020247"))+
  theme(axis.text.x=element_text(size=6, angle=90,hjust=0.95,vjust=0.2))
no_PE_PPE <- long_outbreak %>% filter(has_PE_PPE == 0)


# plot outbreak trees (hamming distance == SNP distance)
aln <- ape::read.FASTA("data/concat_gblocks_aln.filtered_polymorphic_sites.fasta")
aln = readDNAStringSet("data/concat_gblocks_aln.filtered_polymorphic_sites.fasta")

combined <- aln[str_remove(names(aln),"_assembled")%in%outbreaks$Biosample[outbreaks$Outbreak=="MN_Beef"]|
                str_remove(names(aln),"_assembled")%in%outbreaks$Biosample[outbreaks$Outbreak=="NM_Dairy_07-B"]|
                str_remove(names(aln),"_assembled")%in%outbreaks$Biosample[outbreaks$Outbreak=="CA_Dairy_11-B"]|
                str_remove(names(aln),"_assembled")%in%outbreaks$Biosample[outbreaks$Outbreak=="CA_Dairy_13-A"]
                ]

combined_dist <- stringDist(combined, method="hamming")
clust <- hclust(combined_dist)
combined_tree <- as.phylo(clust)
plot(combined_tree)
write.tree(combined_tree,"data/hamming_tree_outbreaks.tre")


# order of tree nodes
l<-combined_tree$tip.label[combined_tree$edge]
l<-l[!is.na(l)]
l <- str_remove(l, "_assembled")
no_PE_PPE$Biosample <- factor(no_PE_PPE$Biosample, level=rev(l))
no_PE_PPE <- no_PE_PPE %>% ungroup() 
outbreak_wo_PE_PPE <- ggplot(no_PE_PPE) + 
  geom_tile(aes(x=Gene, y=Biosample, fill=as.factor(Present)))+
  scale_fill_manual(values=c("#EEEEEE","#020247"))+
  theme(axis.text.x=element_text(size=6, angle=90,hjust=0.95,vjust=0.2))+
  geom_text(aes(x=Gene, y=65,label=Name,angle = 90),size=2,)
outbreak_wo_PE_PPE
ggsave("plots/Figure5.svg",outbreak_wo_PE_PPE,width=10, height=12)

