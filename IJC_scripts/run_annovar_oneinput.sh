#!/bin/bash
####### description:
# Template for distributed task, where we want to perform a single run
# of our program with multiple CPUs and shared memory.

# Job description:
#SBATCH --job-name="annotation"
#SBATCH --comment="star SBATCH"

# ******** Main parameters: *********
#SBATCH --mem=100gb          # Total amount of memory
# ---- distributed ----> (doParallel,foreach,mclapply) "usually single task multiple workers"
#SBATCH --cpus-per-task=1 # number of cores per task <= 128
#SBATCH --ntasks=2          # number of tasks, (1cpu x task is default)unles MPI comunication
######SBATCH --nodes=1           # number of nodes ( <= ntasks since no communication between tasks)
#SBATCH --time=0-23:10:00   # walltime dd-hh:mm:ss. current max  = 07-00:00:00 (7 days)
#SBATCH -o annotation.out # STDOUT
#SBATCH -e annotation.err # STDERR



#----- Info:
#BASEDIR=/mnt/beefgs/mvilardell/manual_RNAseq
humandb=/ijc/LABS/MERKEL/DATA/PROJECTS/idevillasante/idevillasante/apps/annovar/humandb/

echo "Starting at $(date)"
module load jdk

module load Perl/5.34.0-GCCcore-11.2.0

gunzip "join_genotype.vcf.gz"
perl /ijc/LABS/MERKEL/DATA/PROJECTS/idevillasante/idevillasante/apps/annovar/table_annovar.pl "join_genotype.vcf" ${humandb} -buildver hg38 -out "join_anno" -remove -nastring . -protocol refGene,clinvar_20220320,exac03,avsnp147,dbnsfp30a -operation gx,f,f,f,f --vcfinput --polish



echo "Finished at $(date)"

