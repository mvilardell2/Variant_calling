#####################################################################
###################       BCFTOOLS COMMANDS          ################
#####################################################################


# BCFtools is a command-line tool that provides a width range of commands for manipulating and analyzing VCF files. VCF and BCF files are used to store genetic variation data. 


##############################################################################
##############################################################################
###  BCFTOOLS BGZIP

# Before indexing the VCF file, it has to be compressed with bgzip:
bgzip input_file.vcf > input_file.vcf.gz

# Also with the bcftools view command: 
bcftools view input_file.vcf  -Oz -o input_file.vcf.gz
# view can be substituted by sort if the file is not sorted (sort + compress).


##############################################################################
##############################################################################
###  BCFTOOLS INDEX

bcftools index -t input_file.vcg.gz
# -t for TBI index file
# -c for CSI index file


##############################################################################
###############################################################################
###  BCFTOOLS SORT
bcftools sort input_file.vcf -o input_file_sorted.vcf


##############################################################################
##############################################################################
### BCFTOOLS VIEW

# Simple visualization:
bcftools view input_file.vcf.gz | less
bcftools view -r chr1 input_file.vcf.gz

# To count the number of snps and indels in a VCF file:
# -v: Select specific type of variants
# -V: Exclude types
bcftools view -v snps input_file.vcf.gz | grep -v '^#' | wc -l # Only SNPs
bcftools view -V snps input_file.vcf.gz | grep -v '^#' | wc -l # Exclude SNPs variants
bcftools view -v indels input_file.vcf.gz | grep -v '^#' | wc -l # Count only indels

######## To split or subset samples: # ##########

# View (and further select) information regarding specific samples: 
bcftools view -s SAMPLE_NAME input_file.vcf.gz | less
# can add more samples separated by comma
# can also input a sample file:
bcftools view -S subset_samples.txt input_file.vcg.gz -Oz -o samples.vcf.gz

# To exclude samples: 
bcftools view -s ^VM141 input_file.vcf.gz
bcftools view -s ^exclude.txt input_file.vcf.gz
# can also put in a sample file

# Split by INDIVIDUAL SAMPLES:
bcftools +split -S samples.txt input_file_all.vcf.gz -Oz -o split 
# -o: name of directory where to store all generated files.


###############################################################################
###############################################################################
###  BCFTOOLS QUERY

# This command extracts specific filed from VCF files. 

#### OPTIONS: 
# -f: extract information from columns
# -e and -i: exclude and include
# -r: to select specific genomic regions
# -o: output file
# -H: prints the header of the VCF file


# List all the samples/IDs for which VCF file contains data: 
bcftools query -l input_file.vcf.gz

# Extract chromosome information (and count):
bcftools query -f '%CHROM\n' input_file.vcf.gz | uniq -c

# Extract information from more than one field: 
bcftools query -f '%CHROM %POS %REF %ALT\n' input_file.vcf.gz | head

# View information from a specific chromosome, and chromosome region:  
bcftools view -r chr1 input_file.vcf.gz | less
bcftools view -r chr1:10000-50000 input_file.vcf.gz | less

# Extract information from a specific genomic region: 
bcftools query -f '%CHROM %POS\t' -r chr1:10000-50000 input_file.vcf.gz

# Extract/Filtering (including/excluding) rows: 
bcftools query -i 'INFO/QD>30' -f'%CHROM %POS %REF %ALT %QD\n' input_file.vcf.gz
# -e for exclusion

bcftools query -H -i 'QUAL>20 & INFO/DP>10' -f '%CHROM %POS %REF %ALT %QUAL %INFO/DP\n' input_file.vcf.gz | head 


############################################################################
############################################################################

###  BCFTOOLS FILTER

# Apply filters and moddify the FILTER column
bcftools filter -i 'QUAL>=60' -s 'HighQual' raw_indels.vcf.gz
# -i: including filter condition. Variants that meet this condition will be retained in the output file.
# -s: name of the filter

############################################################################
############################################################################

###  BCFTOOLS stats
bcftools stats input_file.vcg.gz -o input_file.stats


#############################################################################
#############################################################################

###  BCFTOOLS CONCAT

# Concatenate VCF files from the same samples (for example if you have chr in separate files or SNP and indels separately)



#############################################################################
#############################################################################

###  BCFTOOLS MERGE

# Merge multiple individual VCF files (one sample per file) into a single VCF

# Simple way:
bcftools merge input_file1.vcf.gz input_file2.vcf.gz -o merged_files.vcf.gz

# Using wild card
bcftools merge input_file*.vcf.gz -Oz -o merged_files.vcf.gz

# Using a file list
bcftools merge --file-list list_of_file.txt -Oz -o merged_files.vcf.gz

# Merge two files but only one specific region: 
bcftools merge file1.vcf.gz file2.vcf.gz -Oz -r 1 -o merged_files.vcf.gz


#############################################################################
#############################################################################

###  BCFTOOLS ANNOTATE

# Used to add or remove fields from the VCF file 

# Options: 
# -a: defines the annotatin file from which you want to add.
# -c: Columns that you want to add 
# -x: Columns that you want to remove from the input file.


# To remove a field: 
bcftools annotate -x FORMAT input_file.vcf.gz -Oz -o input_file_noformat.vcf.gz
bcftools annotate -x INFO/DP raw_indels.vcf.gz -Oz -o raw_indels_noDP.vcf.gz


# To add a column to the VCF file: 
bcftools annotate input_file_no_id.vcf.gz -a All.vcf.gz  -c ID -Oz -o input_file_added_id.vcf.gz


