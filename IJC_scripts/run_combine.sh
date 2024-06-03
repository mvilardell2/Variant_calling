#!/bin/bash
#SBATCH --job-name=combine_gvcfs     # Job name
#SBATCH --output=combine_gvcfs.out   # Standard output and error log
#SBATCH --error=combine_gvcfs.err    # Error log
#SBATCH --ntasks=1                   # Run on a single CPU
#SBATCH --time=48:00:00              # Time limit hrs:min:sec
#SBATCH --mem=16G                    # Memory limit



echo "Starting at $(date)"
module load jdk 
module load GATK/4.4.0.0-foss-2021b-Java-17.0.2

gatk CombineGVCFs --output "mapping/join.g.vcf.gz" --variant "all.list" --reference /ijc/LABS/MERKEL/DATA/PROJECTS/idevillasante/20221005_MGraupera_PIK3CA-Var/RNAseq/manual/references/GENCODE.v40.fasta 




echo "Finished at $(date)"




