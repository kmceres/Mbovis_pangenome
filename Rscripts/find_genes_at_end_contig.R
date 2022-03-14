setwd("~/PhD/complete_genome_set_bovis_pangenome/gffs_panaroo/panaroo_results/")
library(dplyr)
library(tidyr)
library(ape)
gene_df <- read.csv("gene_data.csv")
genepa <- read.csv("gene_presence_absence_roary.csv")
genepa_long <- genepa %>% gather("Genome","annotation_id",-colnames(genepa)[1:15])
genepa_long <- genepa_long %>% dplyr::select(Gene, No..isolates,Genome,annotation_id) %>% 
  filter(annotation_id != "")
gff_df <- tibble(scaffold_name=character(),start=numeric(),
                 end=numeric(),annotation_id=character())
list <- read.table("../full_gff_list.txt")
list <- c(list$V1)
for (g in 1:length(list)){
  filename=paste("../",list[g],sep="")
  filename=paste(filename,"_no_fasta.gff",sep="")
  gff <- read.gff(filename)
  gff$annotation_id <- sub("ID=","",gff$attributes)
  gff$annotation_id <- sub(";.*","",gff$annotation_id)
  gff <- gff %>% mutate(scaffold_name = seqid) %>%
    dplyr::select(scaffold_name,start,end,annotation_id)
  gff_df <- bind_rows(gff_df, gff)
}

gene_data_gff <- left_join(gene_df, gff_df)
gene_data_gff <- left_join(gene_data_gff, genepa_long)
gene_data_gff <- gene_data_gff %>% filter(is.na(Gene)==F)
gene_data_gff$contig_length <- sub(".*length_","",gene_data_gff$scaffold_name)
gene_data_gff$contig_length <- as.numeric(sub("_.*","",gene_data_gff$contig_length))
gene_data_gff <- gene_data_gff %>% mutate(right_tail = contig_length - end) %>% 
  mutate(left_tail = start) %>% 
  mutate(at_end = ifelse(right_tail < 500 | left_tail < 500,1,0))
gene_data_gff_f <- gene_data_gff %>% filter(at_end == 1) %>%
  group_by(Gene,Genome)  %>% sample_n(1) %>% 
  ungroup() %>% group_by(Gene) %>% add_tally()
write.csv(gene_data_gff, "gene_data_gffs.csv")
always_at_end_500 <- gene_data_gff_f %>% filter(n == No..isolates) %>% ungroup() %>% distinct(Gene,.keep_all=T)
always_at_end_500 <- always_at_end_500 %>% dplyr::select(Gene, No..isolates)
write.csv(always_at_end_500,"genes_at_end_contig.csv")
