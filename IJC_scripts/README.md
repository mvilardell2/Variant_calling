This directory contains the scripts used for calling variants with GATK. 

Analysis are managed using SLURM job submission scripts. 


### MAPPING

DATE: 16/05/2024

- Alignment is performed with STAR (2.7.6a STAR/2.7.6a-GCC-11.2.0). It aligns RNA-seq reads to a reference genome.

- Need to prepare 'sample_sheet.csv' file, that contains the sample ID,the path to fastq files,and strandess for each of the samples to analyse.

- For each sample, a directory with the alignment outputs was created in /mnt/beegfs/mvilardell/manual_RNAseq/mapping

Code script: run_Array_star.sh  

Bam files generated have been indexed and collected statistics with SAMtools. Code script: run_samtools_idx.sh

### MARKDUPLICATES

Duplicate reads are taged in the bam file. Also generated metrics.
Code script: run_markd.sh


### SPLITNCIGARREADS.

Code script: run_splitncigarreads.sh

### RECALIBRATION.

Code script: run_recallibrate.sh


### VARIANT CALLING 

DATE 18/05/2024

Done with GATK4 Haplotypecaller. For each sample there is a VCF file created within its directory. The output is a VCF.

Code script: run_varcalling.sh.




DATE 22/05/2024

### ANOTATION

Variant annotation was performed using ANNOVAR software in each sample-vcf file.

ANNOVAR consists of two command steps: 

- 1) annotate_variation.pl: Downloads the databases files that are used for the annotation: refGene,clinvar_20220320,exac03,avsnp147,dbnsfp30a.

- 2) table_annovar.pl: Variant annotation   

Code script: run_annovar.sh


### MERGE VCF

BCFtools merge is used to join all vcf single-sample files into a multi-sample vcf file. 

Output: merged_all_files_annotation.vcf.gz

Then, this file is indexed with bcftools

Code script: run_mergevcf.sh



### VARIANT CALLING

DATE 27/05/24

Variant calling was re-done with GATK HaplotypeCaller, obtaining a GVCF file for each sample. 

### COMBINE GVCF

All 18 GVCF files were combined into a single GVCF file using GATK CombineGVCFs. Output file: join.g.vcf.gz

Code script: run_combine.sh


### GENOTYPE GVCF

Join genotyping (of all 18 samples) was performed. The output is a VCF file (no annotation yet). Output: join_genotype.vcf.gz

Code script: run_genotype.sh


### ANNOTATION: 

Variant annotation was performed with ANNOVAR in the multi-sample vcf file obtained in the previous step. Output: join_anno.hg38_multianno.vcf

Code script: run_annovar_oneinput.sh
