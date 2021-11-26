library(ape)
library(IRanges)
library(dplyr)
library(tidyr)
library(stringr)

recomb <- read.gff("data/concat_gblocks_aln.recombination_predictions.gff")
header<- read.table("data/header.tsv",sep="\t", header=T)
header<- header %>% mutate(end = cumsum(Sequence_Length)) %>% 
  mutate(start = end - Sequence_Length + 1)
header$Gene <- str_remove(header$Name, ".aln")
header <- header %>% select(Gene, start, end, Percent_GC)

recomb_R <- IRanges(start=as.numeric(recomb$start),
                    end=as.numeric(recomb$end),
                    names=recomb$attributes)
header_R <- IRanges(start=as.numeric(header$start),
                    end = as.numeric(header$end),
                    names=header$Gene)

pairs <- findOverlapPairs(header_R, recomb_R)
recomb_df <- as.data.frame(pairs)
recomb_df$neg_lik <- sub(".*neg_log_likelihood=","",recomb_df$second.names)
recomb_df$node <- sub(".*node=","",recomb_df$second.names)
recomb_df$node <- str_remove(recomb_df$node,";.*")
recomb_df<-recomb_df %>% mutate(start=second.start, end=second.end) %>%
  select(first.names, start, end, node)
write.csv(recomb_df, "data/recombination_df.csv")

nrow(recomb_df %>% distinct(first.names)) # 100 genes have recomb events
recomb_df <- recomb_df %>% 
  mutate(internal = ifelse(str_detect(recomb_df$node,"->internal"),1,0))

internal <- recomb_df %>% filter(internal == 1) %>% distinct(first.names,.keep_all=T)
nrow(internal %>% distinct(first.names))
# 40 have recombination events on internal nodes (more likely true recombination events)
write.csv(internal,"data/recomb_events_internal_nodes.csv")
