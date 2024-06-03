#!/bin/bash
####### description:
# Template for distributed task, where we want to perform a single run
# of our program with multiple CPUs and shared memory.

# Job description:
#SBATCH --job-name="star"
#SBATCH --comment="star SBATCH"

# ******** Main parameters: *********
#SBATCH --mem=40G         # Total amount of memory 
# ---- distributed ----> (doParallel,foreach,mclapply) "usually single task multiple workers"
#SBATCH --cpus-per-task=1 # number of cores per task <= 128 
#SBATCH --ntasks=22          # number of tasks, (1cpu x task is default)unles MPI comunication  
#SBATCH --nodes=1           # number of nodes ( <= ntasks since no communication between tasks)
#SBATCH --time=1-23:10:00   # walltime dd-hh:mm:ss. current max  = 07-00:00:00 (7 days)
#SBATCH --array=1-18
#SBATCH -o slurm.%N.%J.%u.%a.out # STDOUT
#SBATCH -e slurm.%N.%J.%u.%a.err # STDERR
#----- Info: 

echo "Starting at $(date)"
# ------ Environment Configuration ------
cd $SLURM_SUBMIT_DIR

#Analysis are runned from:
#BASEDIR=/mnt/beegfs/mvilardell/manual_RNAseq/


#n=$(echo "scale=2 ; ${SLURM_ARRAY_TASK_ID}/100" | bc)
ID=$(awk "FNR == $SLURM_ARRAY_TASK_ID" sample_sheet.csv | cut -d, -f1)
R1=$(awk "FNR == $SLURM_ARRAY_TASK_ID" sample_sheet.csv | cut -d, -f2)
R2=$(awk "FNR == $SLURM_ARRAY_TASK_ID" sample_sheet.csv | cut -d, -f3)



module load STAR/2.7.6a-GCC-11.2.0 



STAR --runThreadN $SLURM_NTASKS --sjdbOverhang 99 --genomeDir /ijc/LABS/MERKEL/DATA/PROJECTS/idevillasante/20221005_MGraupera_PIK3CA-Var/RNAseq/manual/references/star_index --sjdbGTFfile /ijc/LABS/MERKEL/DATA/PROJECTS/idevillasante/20221005_MGraupera_PIK3CA-Var/RNAseq/manual/references/GENCODE.v40_anno.gtf --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --quantMode GeneCounts --outFileNamePrefix ./mapping/$ID/ --readFilesCommand zcat --readFilesIn ${R1} ${R2} 


#STAR --runThreadN 2\
#--sjdbOverhang 99 \
#--genomeDir /ijc/LABS/MERKEL/DATA/PROJECTS/idevillasante/20221005_MGraupera_PIK3CA-Var/RNAseq/manual/references/star_index \
#--sjdbGTFfile /ijc/LABS/MERKEL/DATA/PROJECTS/idevillasante/20221005_MGraupera_PIK3CA-Var/RNAseq/manual/references/GENCODE.v40_anno.gtf \
#--outSAMtype BAM SortedByCoordinate \
#--outReadsUnmapped Fastx \
#--quantMode GeneCounts \
#--outFileNamePrefix ./star/ \
#--readFilesCommand zcat \
#--outSAMattrRGline ID:VM141 , ID:VM211 \
#--readFilesIn /ijc/LABS/GRAUPERA/RAW/ANE_MARTINEZ_LARRINAGA/MG07_Illumina_TotalRNASeq_SANDRA/FASTQs/VM141_1.fastq.gz,/ijc/LABS/GRAUPERA/RAW/ANE_MARTINEZ_LARRINAGA/MG07_Illumina_TotalRNASeq_SANDRA/FASTQs/VM211_1.fastq.gz /ijc/LABS/GRAUPERA/RAW/ANE_MARTINEZ_LARRINAGA/MG07_Illumina_TotalRNASeq_SANDRA/FASTQs/VM141_2.fastq.gz,/ijc/LABS/GRAUPERA/RAW/ANE_MARTINEZ_LARRINAGA/MG07_Illumina_TotalRNASeq_SANDRA/FASTQs/VM211_2.fastq.gz 





