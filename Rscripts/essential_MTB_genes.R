library(dplyr)
library(tidyr)
library(stringr)
library(genoPlotR)
mtb_ess <- read.csv("data/Mtb_essential_genes.csv")
mtb_ess <- mtb_ess %>% filter(Final.Call == "ES")
mtb_to_mbov <- read_dna_seg_from_embl("~/PhD/Mbovis_pangenome/First_draft/data/LT708304.1.txt",
                                      tagsToParse ="CDS",extra_fields="note")
mtb_to_mbov <- as.data.frame(mtb_to_mbov)
mtb_to_mbov_n <- mtb_to_mbov %>%  
  separate(note, into=c("MbovID","Gene","Mtb_ID"),sep=",")
n <- mtb_to_mbov_n %>% separate(Mtb_ID, into=c("len","Mtb_ID"),
                                sep="Equivalent to ")
n$Mtb_ID <-  sub("5' end of ","",n$Mtb_ID)
n$Mtb_ID <-  sub("3' end of ","",n$Mtb_ID)
n$Mtb_ID <-  sub("the ","",n$Mtb_ID)
n$Mtb_ID <- sub("5'end of ","", n$Mtb_ID)
n <- n %>% separate(Gene, into=c("Gene","id"),sep="Equivalent to")
n$id <- sub(" len.*","",n$id)
n <- n %>% mutate(Mtb_ID = ifelse(is.na(Mtb_ID)==T,id,Mtb_ID))
n <- n %>% separate(len, into=c("len","extra_id"),sep=" Similar to ")
n$extra_id <- sub("3' end of ","",n$extra_id)
n$extra_id <- sub("5' end of ","",n$extra_id)
n$extra_id <- sub("the ","",n$extra_id)
n <- n %>% mutate(Mtb_ID = ifelse(is.na(Mtb_ID)==T,extra_id,Mtb_ID))
mtb_ess$Mtb_ID <- mtb_ess$ORF.ID

mtb_ess_in_bov <- left_join(n,mtb_ess)
# all but rrs is in bovis
mtb_ess_in_bov <- mtb_ess_in_bov %>% mutate(Mtb_name = Name, 
                                            Mtb_description=Description) %>%
  select(name, start,end, length,gene,synonym,product,proteinid,MbovID,Mtb_ID,
         Mtb_name,Mtb_description,Final.Call)

paralog_l <- read.csv("./data/cytoscape.csv")
paralog_l <- paralog_l %>% filter(paralog == 1)

filtered_genes <- read.table("data/filtered_accessory_genes.txt")
filtered_genes <- colnames(filtered_genes) ; filtered_genes = gsub("\\.","~",filtered_genes)
ref_matches <- read.csv("data/ref_matches.csv")
ref_matches$ref_id <- ref_matches$second.ref_name
ref_matches <- ref_matches %>% filter(first.names %in% filtered_genes) %>% 
  filter(n_ref_corr > 1) 
  
mtb_ess_in_bov$ref_id <- mtb_ess_in_bov$synonym

ref_matches_ess <- left_join(ref_matches, mtb_ess_in_bov, by="ref_id")
ref_matches_ess <- ref_matches_ess %>% 
  mutate(category=ifelse(Final.Call=="ES","ES","not_ES")) %>% 
  mutate(category=ifelse(is.na(category == T),"not_ES",category))  %>% 
  filter(!(Name %in% paralog_l$shared.name)) 

ess <- ref_matches_ess %>% filter(category == "ES") %>% select(first.names)
l <- genepa_long %>% filter(Gene %in% ess$first.names) %>% filter(annotation_id != "")

## evaluating each potential pseudogene essential gene

# essential ref genes that have more than 1 pangenome match 
# group_3101 has a small insertion in a repetitive region of glmS
# group_2806 is near the end of contigs but not within 500 bp
# group_3575 is near the end of contigs, but not within 500 bp
# singletons:
# group_3132 is near the end of its contig
# group_774 is near the end of its contig
# group_1783 is a paralog
# group_367 rearrangement -- not pseudogene 
# group_3295 is a paralog

# in total, the only essential gene that were detected as pseudogenes 
# was group_3101 which has a small insertion in a repetitive region of glmS


