library(dplyr)
library(tidyr)
library(stringr)
library(genoPlotR)
mtb_ess <- read.csv("data/Mtb_essential_genes.csv")
mtb_ess <- mtb_ess %>% filter(Final.Call == "ES")
mtb_to_mbov <- read_dna_seg_from_embl("data/LT708304.1.txt",
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


ref_matches <- read.csv("data/ref_matches.csv")
ref_matches$ref_id <- ref_matches$second.ref_name

mtb_ess_in_bov$ref_id <- mtb_ess_in_bov$synonym


ref_matches_ess <- left_join(ref_matches, mtb_ess_in_bov, by="ref_id")
ref_matches_ess <- ref_matches_ess %>% 
  mutate(category=ifelse(Final.Call=="ES","ES","not_ES")) %>% 
  mutate(category=ifelse(is.na(category == T),"not_ES",category)) 

tab <- table(ref_matches_ess$category, ref_matches_ess$n_ref_corr)
# 4 paralogs and 3 core-rare matches (2 or 1 copies)
# of the 3 core-rare matches, 2 of the rare matches are ends
# glgB (851 samples) matches with group_451 (2 samples)
# eccD3 (851 samples) matches with group_2377 (2 samples)

rare_matches <- genepa_10 %>% filter(Gene == "eccD3"|Gene=="group_2377"|Gene == 
                       "glgB"| Gene == "group_451") 
rare_matches_long <- rare_matches %>% 
  gather("BioSample","Present",-c(colnames(rare_matches)[1:15]))

group_2377 <- rare_matches_long %>% filter(Gene == "group_2377")%>%
  filter(Present ==1 )
group_451 <- rare_matches_long %>% filter(Gene == "group_451")%>%
  filter(Present ==1 )

