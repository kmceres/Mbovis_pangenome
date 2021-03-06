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
* Panaroo individual gene alignments:
	- panaroo-msa --aligner mafft -a core -t 8 --verbose -o . 

* Parsnp core alignment:
* 	- parsnp -g path/h37rv.gb -d ./fasta_dir/ -c
## Recombination analysis
* Gubbins
	- run_gubbins.py parsnp.fasta -s parsnp.treefile -p bovis_pangenome
* FastGEAR
	- run_fastGEAR.sh /path/to/MATLAB_Runtime/v901 /path/to/alignment /path/to/outdir/out.mat ./fG_input_specs.txt
	- (default specs)

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
	- vcffilter -f "MQ > 30 & DP > 10 & AF > 0.9 " file.vcf > file.filtered.vcf
	- bgzip file.filtered.vcf
	- tabix file.filtered.vcf.gz 
	- bcftools merge \*.vcf.gz -Oz -o Merged.vcf.gz 

* snpEFF
	- java -Xmx4g -jar snpEff af2122v4 Merged.vcf.vcf.gz -no-downstream -no-upstream > Merged_annotated.vcf

## indel calling 
* create pileup
	- samtools mpileup -f ../alignment/NC_002945v4.fasta -s $line.bam  -o aln.sorted.bam.mpileup  #ScarTrek has specific naming requirements
* ScarTrek
 	-python find_scars.py -i ./dir_with_bam_files -g AF2122_genes.txt -p AF2122_proteins.txt 


## PANTHER over-representation test
* gene labels for Mtb were used because PANTHER has an Mtb datatabse
* input files:
	- test set: acc_list_for_PANTHER.txt
	- ref set: core_list_for_PANTHER.txt
	- annotation dataset: PANTHER Protein Class
* significance: 
	- FDR p < .05, Fisher's exact test
	- only interested in gene categories over-represented in test set 

