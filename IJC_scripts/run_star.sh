#!/bin/bash
#SBATCH --comment=STAR_smp
#SBATCH --nodes=1

#SBATCH --mem=40G
#SBATCH --output=%j_STAR_VM204.log
#SBATCH --time=3-22:30:00


module load STAR/2.7.6a-GCC-11.2.0
 

STAR --runThreadN 2\
--sjdbOverhang 99 \
--genomeDir /ijc/LABS/MERKEL/DATA/PROJECTS/idevillasante/20221005_MGraupera_PIK3CA-Var/RNAseq/manual/references/star_index \
--sjdbGTFfile /ijc/LABS/MERKEL/DATA/PROJECTS/idevillasante/20221005_MGraupera_PIK3CA-Var/RNAseq/manual/references/GENCODE.v40_anno.gtf \
--outSAMtype BAM SortedByCoordinate \
--outReadsUnmapped Fastx \
--quantMode GeneCounts \
--outFileNamePrefix ./VM204/ \
--readFilesCommand zcat \
--readFilesIn /ijc/LABS/GRAUPERA/RAW/ANE_MARTINEZ_LARRINAGA/MG07_Illumina_TotalRNASeq_SANDRA/FASTQs/VM204_1.fastq.gz /ijc/LABS/GRAUPERA/RAW/ANE_MARTINEZ_LARRINAGA/MG07_Illumina_TotalRNASeq_SANDRA/FASTQs/VM204_2.fastq.gz 

