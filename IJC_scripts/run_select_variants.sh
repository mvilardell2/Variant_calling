#!/bin/bash
####### description:

#SBATCH --job-name="Select variants"

# ******** Main parameters: *********
#SBATCH --mem=100gb          # Total amount of memory
# ---- distributed ----> (doParallel,foreach,mclapply) "usually single task multiple workers"
#SBATCH --cpus-per-task=1 # number of cores per task <= 128
#SBATCH --ntasks=2          # number of tasks, (1cpu x task is default)unles MPI comunication
######SBATCH --nodes=1           # number of nodes ( <= ntasks since no communication between tasks)
#SBATCH --time=0-23:10:00   # walltime dd-hh:mm:ss. current max  = 07-00:00:00 (7 days)
#######SBATCH --array=2-19
#SBATCH -o select_variants.out # STDOUT
#SBATCH -e select_variants.err # STDERR



#----- Info:
#BASEDIR=/mnt/beefgs/mvilardell/manual_RNAseq/other_analysis
#samples=/mnt/beefgs/mvilardell/manual_RNASeq/sample_sheet.csv
ref=/ijc/LABS/MERKEL/DATA/PROJECTS/idevillasante/20221005_MGraupera_PIK3CA-Var/RNAseq/manual/references/GENCODE.v40.fasta


echo "Starting at $(date)"
module load jdk
module load GATK/4.4.0.0-foss-2021b-Java-17.0.2



#ID=$(awk "FNR == $SLURM_ARRAY_TASK_ID" ${samples} | cut -d, -f1)
#cd  "$BASEDIR/../mapping/$ID/"


gatk SelectVariants -R ${ref} -V /mnt/beegfs/mvilardell/manual_RNAseq/mapping/Multi-sample/Done_w_combinegvcf/join_genotype.vcf --select-type SNP -O results/raw_snps.vcf
gatk SelectVariants -R ${ref} -V /mnt/beegfs/mvilardell/manual_RNAseq/mapping/Multi-sample/Done_w_combinegvcf/join_genotype.vcf --select-type INDEL -O results/raw_indels.vcf


echo "Finished at $(date)"

