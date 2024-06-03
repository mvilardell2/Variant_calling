# NGS - Variant calling

This repository contain scripts and tools for Next-Generation Sequencing (NGS) variant calling. Variant calling is a crucial step in the analysis of NGS data, where genetic variations such as single nucleotide polymorphisms (SNPs), insertions, deletions, and structural variants are identified and annotated.

## Bioinformatic steps summary:

To start: Reads generated with NGS (FASTQ format).

1. QUALITY CONTROL: FASTQC
2. ALIGNMENT:
      
   Genome of reference---> Create a Genome INDEX
3. SORT ALIGNMENTS and generate metrics.
4. MARK DUPLICATE READS: Flag duplicated reads so in downstream analysis, GATK can ignore these reads. 
5. SPLIT N CIGAR READS
6. BASE QUALITY SCORE RECALIBRATION
7. VARIANT CALLING
   
   - Obtain a GVCF file --> Combine multiple samples into a single file (if any) with COMBINEGVCFs. Then, join genotyping (GENOTYPEGVCFs). The output now will be a VCF file.
   
   - Obtain a VCF file --> Merge multiple samples with bcftools.

8. SELECT SNPS or INDELS

9. VARIANT FILTRATION

   - Variant Quality Score Recalibration and ApplyVQSR
   
   - Hard Filtering using GATK VariantFiltration: filtering based on custom quality criteria such as sequencing depth, mapping quality etc. 
  
10. VARIANT ANNOTATION
  
   


![image](https://github.com/mvilardell2/Variant_calling/assets/154689619/e6b548f8-f166-4b1c-b460-4350bc0bacab)

## TOOLS: 

FASTQC: quality control.

SAMTools: managment of SAM files.

GATK: Tool for the analysis of variant calling. Covers from mark duplicate reads untill variant annotation steps. 

ANNOVAR: software for annotation variants.

BCFtools: Managment of VCF files.
