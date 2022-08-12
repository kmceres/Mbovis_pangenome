library(ggplot2)
library(dplyr)
library(viridis)
library(svglite)
library(tidyr)
library(stringr)
library(ape)
genepa<- read.csv("additional_data/gene_presence_absence_roary.csv")
ngenomes = max(genepa$No..isolates)
gene_df <- read.csv("additional_data/gene_data.csv") ## too large to upload to github; email kc649@cornell.edu if interested
genepa_long <- genepa %>% gather("Genome","annotation_id",-colnames(genepa)[1:15])
genepa_long <- genepa_long %>% dplyr::select(Gene, No..isolates,Genome,annotation_id) %>% 
  filter(annotation_id != "")
gff_df <- tibble(scaffold_name=character(),start=numeric(),
                 end=numeric(),annotation_id=character())
list <- colnames(genepa)[15:ncol(genepa)]
# combine gene presence/absence matrix with location about position of gene within contig
for (g in 1:length(list)){
  filename = paste("./gffs",list[g],sep="/")
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
write.csv(gene_data_gff, "additional_data/gene_data_gffs.csv")
always_at_end_500 <- gene_data_gff_f %>% filter(n == No..isolates) %>% ungroup() %>% distinct(Gene,.keep_all=T)
always_at_end_500 <- always_at_end_500 %>% dplyr::select(Gene, No..isolates)
write.csv(always_at_end_500,"additional_data/genes_at_end_contig.csv")



## ref matches
library(dplyr)
library(IRanges)
library(stringr)
library(igraph)
library(ape)
# blasted pan genome reference genes to AF2122/97 (NC002945)
blast_results <- read.csv("additional_data/blast_results.csv")
genepa <- read.csv("additional_data/gene_presence_absence_roary.csv")
genes_not_blasted <- genepa %>% filter(!(Gene %in% blast_results$qseqid))
ends <- read.csv("additional_data/genes_at_end_contig.csv")
ends$End = 1
genes_not_blasted <- left_join(genes_not_blasted,ends,by=c("Gene"))

max_full_blast <- blast_results %>% group_by(qseqid) %>% 
  mutate(max_bit = max(bitscore)) %>%
  filter(bitscore == max_bit)
max_full_blast <- max_full_blast %>% ungroup() %>% 
  mutate(dif=send-sstart) %>% 
  mutate(start=ifelse(dif<0, send, sstart)) %>% 
  mutate(end=ifelse(dif<0, sstart, send))
max_full_blast <- max_full_blast %>% filter(is.na(qseqid)==F)

# read in annotation for af2122.97
an <- read.gff("data/GCA_000195835.3_ASM19583v2_genomic.gff")
an$Gene <- sub(";.*","",an$attributes)
an$Gene <- sub("ID=gene-","",an$Gene)
an$Name <- sub(";gbkey.*","",an$attributes)
an$Name <- sub(".*Name=","",an$Name)
an <- an %>% filter(type=="gene")


# find overlapping ranges between blast matches and annotations
ann_R<- IRanges(start= as.numeric(an$start),
                end = as.numeric(an$end),
                names= an$Name,
                ref_name = an$Gene)

blast_R <- IRanges(start=as.numeric(max_full_blast$start),
                   end=as.numeric(max_full_blast$end),
                   names=max_full_blast$qseqid)


o<- findOverlapPairs(blast_R, ann_R)
o_df <-as.data.frame(o)
o_df <- o_df %>% dplyr::select(first.start, first.end, second.X.start, second.X.end, first.names, second.ref_name) %>%
  group_by(second.ref_name) %>% add_tally(name="n_ref")

## for remaining genes,
## find any remaining pileup genes that weren't covered in self match
##    i.e. those genes that overlap the same ref but not one another
## filter out short overlaps (not real gene overlaps)
## and categorize:
##  core = genes that overlap ref genes completely + perfectly
##  pseudogenes = genes that are within ref genes, but shorter
##  pseudogenes + full = genes that sometimes are core, but sometimes are pseudogenes
ref_matches <- o_df %>%
  # filter(!(first.names %in% pileup_genes_to_remove)) %>%
  mutate(gene_within_ref = ifelse(second.X.start <= first.start & second.X.end >= first.end, 1, 0)) %>%
  mutate(query.width = first.end - first.start + 1, ref.width = second.X.end - second.X.start + 1) %>%
  mutate(diff_starts = second.X.start - first.start) %>%
  mutate(diff_ends = second.X.end - first.start) %>%
  mutate(overlap = case_when(
    gene_within_ref == 1 ~ "within",
    diff_starts > 0 ~ "left",
    diff_ends > 0 ~ "right"
  )) %>% 
  mutate(overlap_len = case_when(
    overlap == "left" ~ first.end - second.X.start,
    overlap == "right" ~  second.X.end - first.start,
    overlap != "within" & ( diff_starts == 0 | diff_ends == 0) ~ as.integer(0)
  )) %>% 
  mutate(overlap_percent = ifelse(overlap == "within", query.width/query.width,
                                  overlap_len/query.width)) %>%
  filter(overlap_percent > .75) %>%
  group_by(first.names) %>% 
  add_tally(name="n_gene") %>% ungroup() %>%
  group_by(second.ref_name)%>% add_tally(name="n_ref_corr") %>%
  mutate(len_same = ifelse(query.width == ref.width, 1,0))%>%
  mutate(same_start = ifelse(first.start == second.X.start, 1, 0)) %>%
  mutate(is_core = ifelse(len_same == 1 & same_start == 1, 1,0)) %>%
  mutate(left_overlap = first.start - second.X.start) %>%
  mutate(right_overlap = first.end - second.X.end) %>% 
  mutate(gene_within_ref = ifelse(first.start >= second.X.start & first.end <= second.X.end, 1, 0)) %>%
  group_by(second.ref_name) %>%
  mutate(ref_gene_has_core = ifelse(sum(is_core > 0),1,0)) %>%
  mutate(pseudo_status = case_when(
    gene_within_ref == 1 & n_ref_corr == 1 & len_same == 0~ "pseudogene",
    gene_within_ref == 1 & ref_gene_has_core == 1 & n_gene > 1~"pseudogene + full",
    is_core == T ~ "core",
    left_overlap <=0 & right_overlap >= 0 ~ "gene spans ref",
    gene_within_ref == 1 & n_ref_corr > 1 ~ "gene within, multiple genes",
    left_overlap <=0 & right_overlap <= 0 | left_overlap >= 0 & right_overlap >=0 ~"gene on one side ref",
    TRUE ~ "0")
  ) 
an$second.ref_name <- an$Gene
ref_matches <- left_join(ref_matches, an, by="second.ref_name")
write.csv(ref_matches, "additional_data/ref_matches.csv")


##
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



## get barplot comparing gene counts
# unfiltered: counts without filtering redundant genes
unfilt_genes <- genepa %>% dplyr::select(Gene, No..isolates)
unfiltered <- tibble(soft_core = nrow(unfilt_genes %>% 
                                filter(No..isolates >= .95*ngenomes)),
                     shell = nrow(unfilt_genes %>% 
                                filter(No..isolates >= .15 * ngenomes & 
                                         No..isolates < .95 * ngenomes)),
                     cloud = nrow(unfilt_genes %>% filter(No..isolates < 3)),
                     group = "unfiltered")

# filtered: remove genes at ends of contigs, and those linked to genes at ends
ref_genes_to_remove <- gene_match_df %>% 
  filter(n_gene_match > 1 | n_genes_at_end > 0)

genes_to_remove <- ref_matches %>% filter(second.ref_name %in% ref_genes_to_remove$ref |
                                            str_detect(Name,"PE_")==T | str_detect(Name,"PPE")==T)




filt_genes <- genepa %>% dplyr::select(Gene, No..isolates) %>%
  filter(!(Gene %in% ends$Gene)) %>%
  filter(!(Gene %in% genes_to_remove$Gene)) %>%
  filter(str_detect(Gene,"PE_")==F & str_detect(Gene,"PPE")==F)
filtered <- tibble(soft_core = nrow(filt_genes %>% 
                                filter(No..isolates >= .95*ngenomes)),
                     shell = nrow(filt_genes %>% 
                                filter(No..isolates >= .15 * ngenomes & 
                                         No..isolates < .95 * ngenomes)),
                     cloud = nrow(filt_genes %>% filter(No..isolates < 3)),
                     group = "filtered")
reis <- tibble(soft_core = 3708,
                 shell = 1341,
                 cloud =2656 ,
                 group="Reis & Cunha 2021" )

combined <- bind_rows(unfiltered,filtered,reis)
combined <- combined %>% gather("category","count",-group)
combined$group <- factor(combined$group, levels=c("unfiltered","filtered", "Reis & Cunha 2021"))
combined$category <- factor(combined$category, levels=c("soft_core","shell","cloud"))



p2 <- ggplot(combined,aes(x=category, y= count,fill=group))+
  geom_col(position = "dodge")+
  geom_text(aes(label=count), 
            position = position_dodge(width=.9),
            vjust=0,color="black") +
  scale_fill_viridis(discrete=T)+
  labs(y="Number of genes in category", x = NULL, fill=NULL)+
  theme(legend.position = "bottom",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


ggsave("Pangenome_compare.svg",p2, width=5,height=5)


# reis and cunha estimate
# core = all = 2736
# soft core  > 95%  =  3708
# acc = shell + cloud = rare + only one or two 
# shell = 1341
# cloud = 2656 (singletons/doubletons)
