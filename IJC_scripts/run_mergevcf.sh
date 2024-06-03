#!/bin/bash

# Job description:
#SBATCH --job-name="Merge_vcf"

#SBATCH --output=logs/merge.out      # Standard output and error log
#SBATCH --error=logs/merge.err       # Error log
#SBATCH --time=01:00:00                 # Wall time limit
#SBATCH --mem=16G                       # Memory per node



module load BCFtools



echo "Starting at $(date)"

#BASEPATH=/mnt/beegfs/mvilardell/manual_RNAseq

#ID=$(awk "FNR == $SLURM_ARRAY_TASK_ID" sample_sheet.csv | cut -d, -f1)
#cd  "$SLURM_SUBMIT_DIR/mapping/$ID/"


# Compress the file: 
#bgzip -c ${ID}_raw_variants.vcf > ${ID}_raw_variants.vcf.gz

# Create an index: 
#bcftoos index -t ${ID}_raw_variants.vcf.gz 


#----------------------------

#BASEPATH=/mnt/beegfs/mvilardell/manual_RNAseq

#cd "$BASEPATH"

# Go to manual_RNAseq directory and run: 

bcftools merge --file-list all_list.txt -Oz -o mapping/merged_all_files_annotation.vcf.gz  




echo "Finished at $(date)"




