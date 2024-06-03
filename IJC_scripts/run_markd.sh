#!/bin/bash
####### description:
# Template for distributed task, where we want to perform a single run
# of our program with multiple CPUs and shared memory.

# Job description:
#SBATCH --job-name="markdup"
#SBATCH --comment="star SBATCH"

# ******** Main parameters: *********
#SBATCH --mem=100gb          # Total amount of memory
# ---- distributed ----> (doParallel,foreach,mclapply) "usually single task multiple workers"
#SBATCH --cpus-per-task=1 # number of cores per task <= 128
#SBATCH --ntasks=20          # number of tasks, (1cpu x task is default)unles MPI comunication
#SBATCH --nodes=1           # number of nodes ( <= ntasks since no communication between tasks)
#SBATCH --time=1-23:10:00   # walltime dd-hh:mm:ss. current max  = 07-00:00:00 (7 days)
#SBATCH --array=1-19
#SBATCH -o slurm.%N.%J.%u.%a.out # STDOUT
#SBATCH -e slurm.%N.%J.%u.%a.err # STDERR
#----- Info:

echo "Starting at $(date)"
module load jdk 
module load GATK/4.4.0.0-foss-2021b-Java-17.0.2


ID=$(awk "FNR == $SLURM_ARRAY_TASK_ID" sample_sheet.csv | cut -d, -f1)
cd  "$SLURM_SUBMIT_DIR/mapping/$ID/"

input=Aligned.sortedByCoord.out.bam
gatk --java-options "-Xmx46g" MarkDuplicates I= $input O=marked_duplicates.bam M=marked_dup_metrics.txt





