# Scripts and accessory data files for: A critical evaluation of *Mycobacterium bovis* pangenomics, with reference to its utility in outbreak investigation
Maintained by Kristina Ceres: kc649@cornell.edu

## Rscripts
* world_plot_cluster.R: scripts to create geographic clusters and map plot for **Figure 1**
* PCA.R: Core and accessory genome clustering and plots for **Figure 2**
* ref_blast_matches.R: analyzing where M. bovis genes identifed by Panaroo map to the M. bovis reference af2122/97
* recombination_analysis.R: organizing results from Gubbins for manual contig alignment inspection to determine whether recombination event is caused by error
* characterizing_acc_genes.R: plots for **Figure 3** related to differnet types of accessory genes
* essential_MTB_genes.R: determining if any pseudogenes were essential in Mtb
* outbreak_analysis.R: looking at variable genes over the four outbreaks included in the study, and plot for **Figure 5**

## iTOL annotation files
* ITOL_colorstrip_kmeans_cluster.txt: Kmeans labels for Jaccard dendrogram
* ITOL_mbovis_group_colorstrip.txt: M. bovis phylogenetic labels 
* ITOL_glpK.txt: heatmap showing presence and absence of glpK 7C-->8C insertion used in **Figure 4**

## tree files
* unfiltered_acc_gene_dendrogram.nwk: Jaccard distance dendrogram used in **Figure2A**
* concat_gblocks_aln.node_labelled.final_tree.tre: core phylogeny used in all figures that display core trees
* hamming_tree_outbreaks.tre: SNP tree for **Figure 5** 




