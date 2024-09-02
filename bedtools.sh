###################################################################
################          BEDTOOLS          #######################
###################################################################

# Tools used to manipulate genomic interval files.


##################################################################
#################################################################

###  BEDTOOLS INTERSECT

# Identifies all the regions in the genome where the features in the two files overlap.

# Options: 
# -wa: returns the original records from A file that overlap
# -wb: returns the original records from B file that overlap
# -wo: report the number of base pairs that overlap
# -c: count the number of overlapping features
# -v: Idenfify regions that do not overlap


# By default you can see where the intersections occurred: 
bedtools intersect -a cpg.bed -b exons.bed | head -5

# To see which regions from A file are overlaped: 
bedtools intersect -a cpg.bed -b exons.bed -wa | head -n 5
# use the -wa options

# To see what overlaped in both files: 
bedtools intersect -a cpg.bed -b exons.bed -wa -wb | head -n 5

# To see how many base pairs overlap: 
bedtools intersect -a cpg.bed -b exons.bed -wo | head

# Find features that do not overlap: 
bedtools intersect -a cpg.bed -b exons.bed -v | head

# To specify what fraction of each feature in A should be overlapped in B: (50% of the region of A overlaps to B) 
bedtools intersect -a cpg.bed -b exons.bed -wo -f 0.50 | head

# intersect with many files: 
bedtools intersect -a exons.bed -b cpg.bed gwas.bed hesc.chromHmm.bed -sorted -wa -wb -names cpg gwas chromhmm | head -n 10000 | tail -n 10


#################################################################
################################################################

###  BEDTOOLS WINDOW

# Allows to specify a number of base pairs upstream and downstream for overlapping features

bedtools windows -a cpg.bed -b exons.bed -w 200
# -l: specify upstream 
# -r: specify number of basepairs downstream


#################################################################
#################################################################

###  BEDTOOLS CLOSEST 

# Search for overlapping features, but in the event that no feature B overlaps in A,search for the closest one.  
bedtools closest -a exons.bed -b gwas.bed -d | head
# -d: reports for the distance

##################################################################
##################################################################

###  BEDTOOLS MERGE

# We can have a file with many individual features that overlap one another, then we can combine them into a single and contiguous interval. 

bedtools merge -i exons.bed 

# Report the number of intervals that are integrated into the new, merged file: 
bedtools merge -i exons.bed -c 1 -o count | head
# -c: indicate the column that you want to summarize (in this case we use the chromosome to know if there is more than 1 interval)
# -o: operation to apply in the -c option. In this case count the number of chr that we have. 

# To keep track of the intervals that were merged: 
bedtools merge -i exons.bed -d 90 -c 1,4 -o count,collapse | head
# -d: merge intervals that do not overlap, but are close one another. 


###################################################################
###################################################################

###  BEDTOOLS COMPLEMENT 

# To know which interval of a file are NOT covered
bedtools complement -i exons.bed -g genome.txt | head


##################################################################
##################################################################

###  BEDTOOLS GENOMECOV

# Measure genome-wide coverage
bedtools genomecov -i exons.bed -g genome.txt



