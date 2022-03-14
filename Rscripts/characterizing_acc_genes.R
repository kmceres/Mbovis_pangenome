library(ggplot2)
library(dplyr)
library(viridis)
library(svglite)
library(tidyr)
library(stringr)
genepa<- read.csv("data/gene_presence_absence_roary.csv")
ngenomes = max(genepa$No..isolates)
phage_l <- c("group_1979", "group_616","group_2566", "group_2533", "group_2405",
             "group_2345", "group_601",  "group_373",  "group_2670","group_1978")

ends <- read.csv("data/genes_at_end_contig.csv")
ref_matches <- read.csv("data/ref_matches.csv")
ref_matches <- ref_matches %>% mutate(Gene = first.names) %>% group_by(first.names) %>%
  mutate(max_percent_match=max(overlap_percent) )%>%
  filter(overlap_percent==max_percent_match) %>% ungroup() %>%
  distinct(first.names,.keep_all=T)


gpa_long <- genepa %>% gather("BioSample","Present",colnames(genepa)[15:ncol(genepa)])
gpa_long <- gpa_long %>% dplyr::select(Gene, No..isolates, Order.within.Fragment, BioSample, Present)
matches <- ref_matches %>% group_by(second.ref_name) %>%
  filter(n_ref_corr >= 2)
gpa <- genepa %>% dplyr::select(Gene, No..isolates)
matches <- left_join(matches, gpa)
matches <- matches %>% filter(No..isolates < ngenomes) %>% 
  ungroup() %>% 
  group_by(second.ref_name) %>%
  add_tally(name="n_acc_gene_match") %>% filter(n_acc_gene_match > 1)
ref_l <- unique(matches$second.ref_name)
ref_df <- matches %>% ungroup() %>% mutate(ref = second.ref_name) %>% 
  dplyr::select(ref, Name) %>% distinct(ref,.keep_all=T)


gene_match_df <- tibble(ref = ref_l, 
                        n_gene_match = numeric(length=length(ref_l)),
                        avg_match_gene_per_sample = numeric(length=length(ref_l)),
                        n_genes_at_end = numeric(length=length(ref_l)))

for (g in 1:length(ref_l)){
  genes <- matches$Gene[matches$second.ref_name==ref_l[g]]
  samp <- gpa_long %>% filter(Gene %in% genes)
  samp <- samp %>% filter(Present != "") %>% group_by(BioSample) %>% 
    add_tally() %>%
    mutate(gene_per_samp = n / length(genes))
  gene_match_df$ref[g] = ref_l[g]
  gene_match_df$n_gene_match[g] = length(genes)
  gene_match_df$avg_match_gene_per_sample[g] = mean(samp$gene_per_samp)
  gene_match_df$n_genes_at_end[g] = sum(genes %in% ends$Gene == T)
}


gene_match_df <- gene_match_df %>% 
  mutate(percent_genes_at_end = n_genes_at_end / n_gene_match)
gene_match_df <- left_join(gene_match_df, ref_df)
gene_match_df <- gene_match_df %>% filter(ref != "1BQ2027_IS1534-2") 
order <- gene_match_df %>% arrange(n_gene_match)
gene_match_df$Name <- factor(gene_match_df$Name, levels=order$Name)

gene_match_df <- gene_match_df %>% mutate(has_end = ifelse(n_genes_at_end > 0, 1, 0))
gene_match_df_f <- gene_match_df %>% filter(n_gene_match > 2)
order <- gene_match_df_f %>% arrange(n_gene_match)
gene_match_df_f$Name <- factor(gene_match_df_f$Name, levels=order$Name)

p <- ggplot(gene_match_df_f) + geom_col(aes(x=Name,y=n_gene_match, 
                                            fill=avg_match_gene_per_sample))+
  scale_fill_viridis()+#coord_flip()+
  theme_bw()+theme(legend.position="bottom", text = element_text(size=8), 
                   axis.text.x = element_text(angle = 90,hjust=1),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_point(aes(x=Name,y=n_genes_at_end))+
  labs(y="Number of genes mapped to reference gene", x=NULL,
       fill = "Average fraction of mapped genes contained in one genome")

ggsave("Figure3A.svg",p, width = 8, height= 4)

ref <- ref_matches %>% mutate(ref = second.ref_name) %>% select(ref, Name, first.names)
high_n <- gene_match_df_f %>% filter(n_gene_match > 7) %>% select(ref, Name)
high_n <- matches %>% filter(second.ref_name %in% high_n$ref)
unique_names <- unique(high_n$Name)
# Write out tables for ITOL figure (Supp figure 3)

for (g in 1:length(unique_names)){
  out <- genepa_10 %>% filter(Gene %in% high_n$first.names[high_n$Name==unique_names[g]]) %>%
    select(-colnames(genepa_10)[c(1,3:15)])
  out_t <- t(out)
  outname = paste(unique_names[g],"csv",sep=".")
  write.csv(out_t, outname,row.names = T,col.names = F,quote=F)
}

## get barplot comparing gene counts
# unfiltered: counts without filtering redundant genes
unfilt_genes <- genepa %>% dplyr::select(Gene, No..isolates)
unfiltered <- tibble(c = nrow(unfilt_genes %>% 
                                filter(No..isolates < ngenomes &
                                         No..isolates >= .95*ngenomes)),
                     m = nrow(unfilt_genes %>% 
                                filter(No..isolates >= .15 * ngenomes & 
                                         No..isolates < .95 * ngenomes)),
                     r = nrow(unfilt_genes %>% filter(No..isolates < .15 * ngenomes)),
                     group = "unfiltered")

# filtered: remove genes at ends of contigs, and those linked to genes at ends
ref_genes_to_remove <- gene_match_df %>% 
  filter(n_gene_match > 1 | n_genes_at_end > 0)

genes_to_remove <- ref_matches %>% filter(second.ref_name %in% ref_genes_to_remove$ref |
                                            str_detect(Name,"PE_")==T | str_detect(Name,"PPE")==T)


palindromic_contig_genes <-read.csv("data/palindromic_genes.csv")
non_blast_matches <- read.csv("data/genes_not_blasted.csv")
non_blast_matches <- non_blast_matches %>% filter(Gene != "esxM" & Gene != "espI_2")
filt_genes_no_matches_ends <- genepa %>% dplyr::select(Gene, No..isolates) %>% 
  filter(!(Gene %in% ends$Gene)) %>%
  filter(!(Gene %in% non_blast_matches$Gene))%>%
  filter(!(Gene %in% genes_to_remove$Gene)) 
filt_genes_no_matches_ends_ppe <- filt_genes_no_matches_ends %>% 
  filter(str_detect(Gene,"PE_")==F & str_detect(Gene,"PPE")==F)
filt_genes_no_matches_ends_ppe_paralogs  <- filt_genes_no_matches_ends_ppe %>% 
  filter(!(Gene %in% palindromic_contig_genes$Gene))


filt_genes <- genepa %>% dplyr::select(Gene, No..isolates) %>%
  filter(!(Gene %in% ends$Gene)) %>%
  filter(!(Gene %in% genes_to_remove$Gene)) %>%
  filter(!(Gene %in% palindromic_contig_genes$Gene))%>%
  filter(!(Gene %in% non_blast_matches$Gene))%>%
  filter(str_detect(Gene,"PE_")==F & str_detect(Gene,"PPE")==F)


filtered<- tibble(c = nrow(filt_genes %>%
                                filter(No..isolates >= .95*ngenomes &
                                         No..isolates < ngenomes)),
                     m = nrow(filt_genes %>%
                                filter(No..isolates >= .15 * ngenomes &
                                         No..isolates < .95 * ngenomes)),
                     r = nrow(filt_genes %>% filter(No..isolates < .15 * ngenomes)),
                     group = "Multiple matches non-blast matches filtered")

                                           
combined <- bind_rows(unfiltered,filtered)
combined <- combined %>% gather("category","count",-group)
combined$group <- factor(combined$group, levels=c("unfiltered","filtered"))

p2 <- ggplot(combined,aes(x=category, y= count,fill=group))+
  geom_col(position = "dodge")+
  geom_text(aes(label=count), 
            position = position_dodge(width=.9),
            vjust=1.5,color="white") +
  scale_fill_manual(values=c("gray20","gray60"))+theme_bw()+
  labs(y="Number of genes in category", x = NULL, fill=NULL)+
  theme(legend.position = "bottom",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("Figure3D.svg",p2, width=4,height=4)


## output matrix of high N genes for ITOL
## write out acc genes with no PE/PPE and no ends
genepa_10 <- read.csv("data/gene_presence_absence_roary_01.csv")
filtered <- filt_genes %>% mutate(first.names = Gene)
filtered <- left_join(filtered, ref_matches, by="first.names")
filtered_gpa <- genepa_10 %>% filter(Gene %in% filtered$first.names)
filtered_gpa <- filtered_gpa %>% filter(No..isolates < .95*ngenomes) 
mat <- filtered_gpa%>%
  select(colnames(filtered_gpa)[16:ncol(filtered_gpa)])
mat <- t(mat)
colnames(mat) <- filtered_gpa$Gene
colSums(mat)
write.table(mat, "data/filtered_accessory_genes.txt",sep=" ",quote=F)

cyto <- read.csv("data/cytoscape.csv")
paralog_l <- cyto %>% distinct(name,.keep_all=T) %>% filter(name %in% filt_genes$Gene) %>%
  filter(paralog == 1)

## lists for PANTHER overrepresentation analysis
essential <- read.csv("data/ref_matches_essential_mtb.csv")
essential$Name <- essential$synonym
essential$Gene <- essential$first.names
essential <- essential %>% select(Name, Gene, Mtb_ID, category)
filt_ref <- left_join(filt_genes, essential)
core <- filt_ref %>% filter(No..isolates == ngenomes)
acc <- filt_ref %>% filter(No..isolates < .95*ngenomes) %>%
  filter(!(Gene %in% paralog_l$name)) %>% 
  filter(!(Gene %in% phage_l)) 

core_l <- core %>% filter(is.na(Name)==F) %>% select(Mtb_ID)
acc_l <- acc %>% filter(is.na(Name)==F) %>% select(Mtb_ID)
core_l$Mtb_ID <- sub("and.*","",core_l$Mtb_ID)
acc_l$Mtb_ID <- sub("and.*","",acc_l$Mtb_ID)

core_l <- bind_rows(core_l,acc_l)
write.table(core_l, "data/core_list_for_PANTHER.txt",row.names = F,col.names = F,quote=F)
write.table(acc_l, "~/PhD/complete_genome_set_bovis_pangenome/acc_list_for_PANTHER.txt",row.names = F,col.names = F,quote=F)

# PANTHER overrepresentation test using acc_l as test set and core_l as ref set
# and using M. tuberculosis as the species
# sig determined by FDR p < .05, fishers exact, one tailed (only looking at over-represented genes)



