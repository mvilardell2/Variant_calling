#!/bin/bash
####### description:
# Template for distributed task, where we want to perform a single run
# of our program with multiple CPUs and shared memory.

# Job description:
#SBATCH --job-name="var_calling"
#SBATCH --comment="star SBATCH"

# ******** Main parameters: *********
#SBATCH --mem=100gb          # Total amount of memory
# ---- distributed ----> (doParallel,foreach,mclapply) "usually single task multiple workers"
#SBATCH --cpus-per-task=1 # number of cores per task <= 128
#SBATCH --ntasks=20          # number of tasks, (1cpu x task is default)unles MPI comunication
#SBATCH --nodes=1           # number of nodes ( <= ntasks since no communication between tasks)
#SBATCH --time=1-23:10:00   # walltime dd-hh:mm:ss. current max  = 07-00:00:00 (7 days)
#SBATCH --array=2-19
#SBATCH -o slurm.%N.%J.%u.%a.out # STDOUT
#SBATCH -e slurm.%N.%J.%u.%a.err # STDERR
#----- Info:


echo "Starting at $(date)"
module load jdk 
module load GATK/4.4.0.0-foss-2021b-Java-17.0.2


ID=$(awk "FNR == $SLURM_ARRAY_TASK_ID" sample_sheet.csv | cut -d, -f1)
cd  "$SLURM_SUBMIT_DIR/mapping/$ID/"


gatk --java-options "-Xmx5g" HaplotypeCaller  --input "clean_{$ID}.bam" --output "${ID}_haplotypecaller.g.vcf.gz" --reference /ijc/LABS/MERKEL/DATA/PROJECTS/idevillasante/20221005_MGraupera_PIK3CA-Var/RNAseq/manual/references/GENCODE.v40.fasta --dbsnp /ijc/LABS/MERKEL/DATA/PROJECTS/idevillasante/20221005_MGraupera_PIK3CA-Var/RNAseq/manual/references/Homo_sapiens_assembly38.known_indels.vcf.gz  --dont-use-soft-clipped-bases true --standard-min-confidence-threshold-for-calling 20.0 -ERC GVCF 




