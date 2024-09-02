######### SAMTOOLS COMMANDS ########
####################################

# Convert SAM to BAM file format:
samtools view -b sample.sam -o sample.bam 

# -b indicates that the output is a bam
# -S indicates the input is a sam


# Display alignments
samtoos view sample.bam | head -n5


# Sort a BAM file (by genomic location)
samtools sort sample.bam -o sample.sorted.bam 


# Index a BAM file
samtools index sample.sorted.bam


# Extract alignments from a particular region (chromosome):
samtools view sample.sorted.bam 20:1.4e6-1.5e6 | head -n5 


# To count the number of aligments given a specific position:
samtools view -c sample.sorted.bam 20:1.4e6-1.5e6 
# -c is used to count 


# See the header (records information regarding the ref genom and how the BAM was processed): 
samtools view -H sample.sorted.bam


#  Display alignments that match specific filtering criteria
samtools view -f INT sample.bam | head
# can use -F to exclude


# Summary of the flags in the BAM file: 
samtools flagstats sample.sorted.bam

# Coverage: 
samtools coverage -r chr1:10000-15000 -m  Aligned.sortedByCoord.out.bam
# m: to display the histogram 
