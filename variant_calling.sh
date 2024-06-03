#!/bin/bash

# Job description:
#SBATCH --job-name="variant_calling"

# ******** Main parameters: *********
#SBATCH --mem=100gb          # Total amount of memory
# ---- distributed ----> (doParallel,foreach,mclapply) "usually single task multiple workers"
#SBATCH --cpus-per-task=2 # number of cores per task <= 128
#SBATCH --ntasks=2          # number of tasks, (1cpu x task is default)unles MPI comunication
######SBATCH --nodes=1           # number of nodes ( <= ntasks since no communication between tasks)
#SBATCH --time=0-23:10:00   # walltime dd-hh:mm:ss. current max  = 07-00:00:00 (7 days)
#SBATCH -o var_calling.out # STDOUT
#SBATCH -e var_calling.err # STDERR




echo "Starting at $(date)"
module load BWA/0.7.17-foss-2021b
module load fastqc-0.11.9-gcc-11.2.0-dd2vd2m
module load SAMtools/1.19.2-foss-2021b
module load GATK/4.4.0.0-foss-2021b-Java-17.0.2




# Script to call germline variants in a human WGS paired end reads 2 X 100bp
# Following GATK4 best practices workflow - https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-
# This script is for demonstration purposes only


# download data
#wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
#wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz


echo "Run Prep files..."

################################################### Prep files (TO BE GENERATED ONLY ONCE) ##########################################################



# download reference files
#wget -P /ijc/LABS/MERKEL/DATA/PROJECTS/mvilardell/other/variant_calling/WGS_var_calling_tutorial/reference https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
#gunzip /ijc/LABS/MERKEL/DATA/PROJECTS/mvilardell/other/variant_calling/WGS_var_calling_tutorial/reference/hg38.fa.gz

# index ref - .fai file before running haplotype caller
#samtools faidx /ijc/LABS/MERKEL/DATA/PROJECTS/mvilardell/other/variant_calling/WGS_var_calling_tutorial/reference


# ref dict - .dict file before running haplotype caller
#gatk CreateSequenceDictionary R=/ijc/LABS/MERKEL/DATA/PROJECTS/mvilardell/other/variant_calling/WGS_var_calling_tutorial/reference/hg38.fa O=/ijc/LABS/MERKEL/DATA/PROJECTS/mvilardell/other/variant_calling/WGS_var_calling_tutorial/reference/hg38.dict


# download known sites files for BQSR from GATK resource bundle
#wget -P /ijc/LABS/MERKEL/DATA/PROJECTS/mvilardell/other/variant_calling/WGS_var_calling_tutorial/reference/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
#wget -P /ijc/LABS/MERKEL/DATA/PROJECTS/mvilardell/other/variant_calling/WGS_var_calling_tutorial/reference/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx



###################################################### VARIANT CALLING STEPS ####################################################################


# directories
ref="/ijc/LABS/MERKEL/DATA/PROJECTS/mvilardell/other/variant_calling/WGS_var_calling_tutorial/reference/hg38.fa"
known_sites="/ijc/LABS/MERKEL/DATA/PROJECTS/mvilardell/other/variant_calling/WGS_var_calling_tutorial/reference/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/ijc/LABS/MERKEL/DATA/PROJECTS/mvilardell/other/variant_calling/WGS_var_calling_tutorial/aligned_reads"
reads="/ijc/LABS/MERKEL/DATA/PROJECTS/mvilardell/other/variant_calling/WGS_var_calling_tutorial/reads"
results="/ijc/LABS/MERKEL/DATA/PROJECTS/mvilardell/other/variant_calling/WGS_var_calling_tutorial/results"
data="/ijc/LABS/MERKEL/DATA/PROJECTS/mvilardell/other/variant_calling/WGS_var_calling_tutorial/data"




# -------------------
# STEP 1: QC - Run fastqc 
# -------------------

echo "STEP 1: QC - Run fastqc"

fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/

# No trimming required, quality looks okay.


# --------------------------------------
# STEP 2: Map to reference using BWA-MEM
# --------------------------------------

echo "STEP 2: Map to reference using BWA-MEM"

# BWA index reference 
bwa index ${ref}


# BWA alignment
bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634.paired.sam




# -----------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# -----------------------------------------

echo "STEP 3: Mark Duplicates and Sort - GATK4"

gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam



# ----------------------------------
# STEP 4: Base quality recalibration
# ----------------------------------


echo "STEP 4: Base quality recalibration"

# 1. build the model
gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table


# 2. Apply the model to adjust the base quality scores
gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file {$data}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam 



# -----------------------------------------------
# STEP 5: Collect Alignment & Insert Size Metrics
# -----------------------------------------------


echo "STEP 5: Collect Alignment & Insert Size Metrics"

gatk CollectAlignmentSummaryMetrics R=${ref} I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam O=${aligned_reads}/alignment_metrics.txt
gatk CollectInsertSizeMetrics INPUT=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam OUTPUT=${aligned_reads}/insert_size_metrics.txt HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf



# ----------------------------------------------
# STEP 6: Call Variants - gatk haplotype caller
# ----------------------------------------------

echo "STEP 6: Call Variants - gatk haplotype caller"

gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf



# extract SNPs & INDELS

gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf






