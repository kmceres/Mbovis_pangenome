library(dplyr)
library(readxl)
t4 <- read.table("./data/Table4_Zwyer_et_al.txt",header=T,sep="\t")
t4$Marker_gene <- t4$Gene
t4 <- t4 %>% filter(KvarQ_informative=="True")
n <- read.csv("./data/mtb_to_mbov.csv")
n$Marker_gene <- n$Mtb_ID
marker_genes <- left_join(t4, n, by="Marker_gene")
marker_genes$Mtb_len = marker_genes$End-marker_genes$Start + 1

out <- marker_genes %>% mutate(CHROM = "NC_000962.3") %>%
  mutate(POS = Position_ref) %>% 
  mutate(REF = ancestral) %>%
  mutate(ALT = derived) %>%
  dplyr::select(CHROM, POS, REF, ALT, PhylogeneticSNP)

out <- out %>% mutate(lineage_group = 
                        ifelse(PhylogeneticSNP == "La1" | 
                               PhylogeneticSNP == "La2"|
                               PhylogeneticSNP == "La3" |
                               PhylogeneticSNP == "La1_La2"|
                               PhylogeneticSNP == "unk4" | 
                               PhylogeneticSNP == "unk5"|
                               PhylogeneticSNP == "unk6", "group1","group2"))
  
write.table(out,"./data/lineage_positions.txt",sep=",",quote=F,row.names=F)

get_lineage_table <- function(vcf_list){
  lineage_table <- tibble(ID = character(), 
                          L1 = character(),
                          L2 = character())
  
  for (v in 1:length(vcf_list)){
    print(vcf_list[v])
    vcf_name <- paste(vcf_list[v],"filtered_hapall.vcf",sep="_")
    tmp_vcf<-readLines(vcf_name)
    tmp_vcf_data<-read.table(vcf_name, stringsAsFactors = FALSE)
    
    # filter for the columns names
    tmp_vcf<-tmp_vcf[-(grep("#CHROM",tmp_vcf)+1):-(length(tmp_vcf))]
    vcf_names<-unlist(strsplit(tmp_vcf[length(tmp_vcf)],"\t"))
    names(tmp_vcf_data)<-vcf_names
    all_vars <- tmp_vcf_data %>% gather("BioSample","Var",-colnames(tmp_vcf_data)[1:9])
    #all_vars <- all_vars %>% mutate(Variant_present = ifelse(Var == "./.:.:.:.:.:.:.:.", 0, 1))
    
    all_vars <- all_vars %>% mutate(CHROM = `#CHROM`) %>% 
      select(CHROM, POS, REF, ALT)
    
    all_vars_f <- left_join(all_vars, out) %>% filter(is.na(PhylogeneticSNP)==F)
    all_vars_f_group1 <- all_vars_f %>% 
      filter(lineage_group == "group1") %>% distinct(PhylogeneticSNP) %>% 
      mutate(group1_n = paste(PhylogeneticSNP, collapse = "~~")) %>% distinct(group1_n)
    all_vars_f_group2 <- all_vars_f %>% 
      filter(lineage_group == "group2") %>% distinct(PhylogeneticSNP)%>% 
      mutate(group2_n = paste(PhylogeneticSNP, collapse = "~~")) %>% distinct(group2_n)
    
   
    new_row <- tibble(ID = vcf_list[v], 
                      L1 = all_vars_f_group1$group1_n, 
                      L2 = all_vars_f_group2$group2_n)
    
    lineage_table <- bind_rows(lineage_table, new_row)
  }
  return(lineage_table)
}

vcf_list <- read.table("path to list of vcf files")
vcf_list <- c(vcf_list$V1)
lineage_table <- get_lineage_table(vcf_list)

write.csv(lineage_table, "data/lineage_table.csv")
