library(GenomicRanges)
library(plyranges)
library(dplyr)
library(IRanges)
library(stringr)
library(igraph)
library(ape)
# blasted pan genome reference genes to AF2122/97 (NC002945)

blast_results <- read.csv("./data/blast_results.csv")
genepa <- read.csv("gene_presence_absence_roary.csv")
genes_not_blasted <- genepa %>% filter(!(Gene %in% blast_results$qseqid))

max_full_blast <- blast_results %>% group_by(qseqid) %>% 
  mutate(max_bit = max(bitscore)) %>%
  filter(bitscore == max_bit)
max_full_blast <- max_full_blast %>% ungroup() %>% 
  mutate(dif=send-sstart) %>% 
  mutate(start=ifelse(dif<0, send, sstart)) %>% 
  mutate(end=ifelse(dif<0, sstart, send))
max_full_blast <- max_full_blast %>% filter(is.na(qseqid)==F)
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

ov <- join_overlap_left(blast_R, ann_R)
ov_df <- as.data.frame(ov)
o<- findOverlapPairs(blast_R, ann_R)
o_df <-as.data.frame(o)
o_df <- o_df %>% select(first.start, first.end, second.X.start, second.X.end, first.names, second.ref_name) %>%
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
write.csv(ref_matches, "ref_matches.csv")
