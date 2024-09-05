#!/bin/bash
####### description:

#SBATCH --job-name="Fitler variants"

# ******** Main parameters: *********
#SBATCH --mem=100gb          # Total amount of memory
# ---- distributed ----> (doParallel,foreach,mclapply) "usually single task multiple workers"
#SBATCH --cpus-per-task=1 # number of cores per task <= 128
#SBATCH --ntasks=2          # number of tasks, (1cpu x task is default)unles MPI comunication
######SBATCH --nodes=1           # number of nodes ( <= ntasks since no communication between tasks)
#SBATCH --time=0-23:10:00   # walltime dd-hh:mm:ss. current max  = 07-00:00:00 (7 days)
#######SBATCH --array=2-19
#SBATCH -o filter_variants.out # STDOUT
#SBATCH -e filter_variants.err # STDERR



#----- Info:
#BASEDIR=/mnt/beefgs/mvilardell/manual_RNAseq/other_analysis
ref=/ijc/LABS/MERKEL/DATA/PROJECTS/idevillasante/20221005_MGraupera_PIK3CA-Var/RNAseq/manual/references/GENCODE.v40.fasta


echo "Starting at $(date)"
module load jdk
module load GATK/4.4.0.0-foss-2021b-Java-17.0.2



gatk VariantFiltration \
	-R ${ref} \
	-V results/raw_snps.vcf \
	-O results/filtered_snps.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"


# Select variants that pass the filter

gatk SelectVariants \
	--exclude-filtered \
	-V results/filtered_snps.vcf \
	-O results/analysis-ready-snps.vcf


echo "Finished at $(date)"

