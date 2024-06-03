#!/bin/bash


# Job description:
#SBATCH --job-name="multiqc"

# ******** Main parameters: *********
#SBATCH --mem=40G         # Total amount of memory 
# ---- distributed ----> (doParallel,foreach,mclapply) "usually single task multiple workers"
#SBATCH --cpus-per-task=1 # number of cores per task <= 128 
#SBATCH --ntasks=22          # number of tasks, (1cpu x task is default)unles MPI comunication  
#SBATCH --nodes=1           # number of nodes ( <= ntasks since no communication between tasks)
#SBATCH --time=1-23:10:00   # walltime dd-hh:mm:ss. current max  = 07-00:00:00 (7 days)
#SBATCH --array=2-19
#SBATCH -o slurm.%N.%J.%u.%a.out # STDOUT
#SBATCH -e slurm.%N.%J.%u.%a.err # STDERR



module load MultiQC

#Analysis are runned from:
#BASEDIR=/mnt/beegfs/mvilardell/manual_RNAseq/

ID=$(awk "FNR == $SLURM_ARRAY_TASK_ID" sample_sheet.csv | cut -d, -f1)
cd  "$SLURM_SUBMIT_DIR/mapping/$ID/"

multiqc .




