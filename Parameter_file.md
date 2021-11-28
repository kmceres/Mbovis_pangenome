# Parameters for bioinformatic analyses

## Genome assembly
* SPAdes:
	- spades.py -1 \*.\_1P.fq.gz -2 \*.\_2P.fq.gz -o /path/to/contig -t 16 --careful 

* Pilon polishing:
	- bwa index /path/to/contig/contigs.fasta 
	- bwa mem -t 14 /path/to/contig/contigs.fasta -r contig_1P.fq.gz contig_2P.fq.gz | samtools view -Sb | samtools sort -@14 -o /path/to/contig/contig.sorted.bam 
	- samtools index /path/to/contig/contig.sorted.bam
	- pilon --genome /path/to/contig/contigs.fasta --bam /path/to/contig/contig.sorted.bam

* annotation: 
	- prokka /path/to/fasta --proteins af2122.97.gb --outdir /path/to/fasta/prokka

## Assembly quality control
* Quast:
	- quast.py path/to/spades/contigs.fasta -r af2122.97.fasta -g af2122.97.gff3 -o path/to/contig/quast_results/
* CheckM:
	- checkm lineage_wf -x fa --threads 8 --tmdir /path/to/fastas/ /checkm_out/

## Pangenome analysis
* Panaroo:
	- panaroo -i \*.gff -o panaroo_results --clean-mode strict --remove-invalid-genes -t 4

## Core genome alignment
* Panaroo core gene alignment:
	- panaroo-msa --aligner mafft -a core -t 8 --verbose -o . 

* Gblocks: ran on each core gene alignment output by panaroo
	- Gblocks_0.91b/Gblocks gene.aln -t=d 

## Recombination analysis
* IQTree (starting tree for Gubbins):
	- iqtree-2.1.2-MacOSX/bin/iqtree2 -s concat_aln.fasta -m GTR -nt AUTO

* Gubbins
	- run_gubbins.py concat_aln.fasta -s concat_gblocks_aln.fasta.treefile -p concat_gblocks_aln

## BLAST-n analysis
* make local blast database:
	- makeblastdb -in GCF_000195835.3_ASM19583v2_genomic.gbff -parse_seqids -blastdb_version 5 -dbype nucl -out bovis_pangenome_ref

* blast accessory genes:
	- blastn -db bovis_pangenome_ref -query pan_genome_reference.fa -out blast_results.out -outfmt 6 

## SNP calling and annotation
* vSNP
	- vSNP_step1.py -r1 \*\_1.fastq.gz -r2 \*\_2.fastq.gz -r Mycobacterium_AF2122 
	- vSNP_step2.py -n -a 
* filtering SNPs, merge vcf files into one for annotation
	- vcffilter -f "MQ > 30 & DP > 10" file.vcf > file.filtered.vcf
	- bgzip file.filtered.vcf
	- tabix file.filtered.vcf.gz 
	- bcftools merge \*.vcf.gz -Oz -o Merged.vcf.gz 

* snpEFF
	- java -Xmx4g -jar snpEff af2122v4 Merged.vcf.vcf.gz -no-downstream -no-upstream > Merged_annotated.vcf 

## PANTHER over-representation test
* gene labels for Mtb were used because PANTHER has an Mtb datatabse
* input files:
	- test set: acc_list_for_PANTHER.txt
	- ref set: core_list_for_PANTHER.txt
	- annotation dataset: PANTHER Protein Class
* significance: 
	- FDR p < .05, Fisher's exact test
	- only interested in gene categories over-represented in test set 

