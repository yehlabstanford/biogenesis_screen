
import sys
import csv
import subprocess
from subprocess import call


### PARAMETERS THAT CAN BE MODIFIED
threads = 20                       # computer threads to run for certain steps
trim5 = 10                         # bowtie - bases to trim from 5' end of reads
trim3 = 15                         # bowtie - bases to trim from 3' end of reads
max_length = 1000                  # bowtie - maximum length between paired end reads (outer-to-outer)
mismatches = 1                     # bowtie - maximum number of mismatches in a read
alignments = 1                     # bowtie - maximum number of possible alignments for a read
min_base_qual = 1                  # samtools mpileup - minimum base quality to consider
ploidy = 1                         # bcftools call - ploidy

min_depth_parent = 1               # bcftools filter - minimum depth to call parental snps
min_depth_sample = 5               # bcftools filter - minimum depth to call sample snps
min_allele_fraction_parent = 0.05  # bcftools filter - minimum non-reference allele fraction to call parental snps
min_allele_fraction_sample = 0.9   # bcftools filter - minimum non-reference allele fraction to call sample snps
###




parentfile = sys.argv[1]    # .csv file of parental samples to analyze (see example)
samplefile = sys.argv[2]    # .csv file of mutant samples to analyze (see example)
referencefile = sys.argv[3] # reference genome in .fasta format

parent_dicty = {}   # Create key for parent samples
sample_dicty = {}   # Create key for mutant samples


# Read through parent and sample ID files, creating a list [ID, read1 file, read2 file] for each sample
# Add that list to the parent or sample list, to create a dictionary from. Remove the header.
with open(parentfile, 'rt') as pfile:
    reader = csv.reader(pfile)
    parent_list = list(reader)
parent_list.pop(0)

with open(samplefile, 'rt') as sfile:
    reader = csv.reader(sfile)
    sample_list = list(reader)
sample_list.pop(0)


# From the parent and sample ID lists, create a dictionary with the ID as the key, and the paired-end files as the value
for parent in parent_list:
    ID = parent.pop(0)
    parent_dicty[ID] = parent

for sample in sample_list:
    ID = sample.pop(0)
    sample_dicty[ID] = sample


# Build the reference genome index, stored as 'REF'
subprocess.run(['bowtie-build', referencefile, 'REF'])


# Go through the parent dictionary to analyze each parent sequentially
# Create filenames for each step using the dictionary key (the parent IDs)
for key in parent_dicty:
    logfile = key + '_log.txt'
    samfile = key + '.sam'
    bamfile = key + '.bam'
    sortbamfile = key + 's.bam'
    dedupfile = key + 'sd.bam'
    covfile = key + 'cov.txt'
    bcffile = key + '.bcf'
    vcffile = key + '.vcf'
    filtvcffile = key + 'filt.vcf'

    # Assign paired-end read files from the dictionary values
    read1, read2 = parent_dicty[key]


    # Create a log file to store info for each sample from the stderr output of each command line step
    with open(logfile, 'wt') as log:

        # Align reads to the reference genome index (REF) using bowtie
        subprocess.run(['bowtie', '-S', '-p', str(threads), '-v', str(mismatches), '-m', str(alignments), '-5', str(trim5), '-3', str(trim3), '-X', str(max_length), 'REF', '-1', read1, '-2', read2, samfile], stderr = log)
        # -p <int> : number of threads to run
        # -v <int> : maximum number of mismatches to allow per read alignment (else discard)
        # -m <int> : maximum number of possible alignments per read (else discard)
        # -5 <int> : trim bases from 5' end of read
        # -3 <int> : trim bases from 3' end of read
        # -X <int> : maximum distance allowed between paired-end reads (outer-to-outer)
        # -1 <str> : input paired-end file 1
        # -2 <str> : input paired-end file 2


    # Append the stderr output log from each step (i.e. number of reads aligned, percent of PCR duplicates removed, etc) to a log file for each sample
    with open(logfile, 'a') as log:


        # Convert the SAM output to a BAM file
        subprocess.run(['samtools', 'view', '-@', str(threads), '-b', '-o', bamfile, samfile], stderr=log)
        # -@ <int> : number of threads to run
        # -b       : BAM output
        # -o <str> : output file name


        # Sort the BAM file
        subprocess.run(['samtools', 'sort', '-o', sortbamfile, bamfile], stderr=log)
        # -o <str> : output file name

        
        # Remove PCR duplicate reads from the sorted BAM file
        subprocess.run(['samtools', 'rmdup', sortbamfile, dedupfile], stderr=log)


        # Create a pileup file from the sorted and de-duplicated BAM file
        subprocess.run(['samtools', 'mpileup', '-I', '-B', '-f', referencefile, '-ug', '-Q', str(min_base_qual), '-t', 'DP,AD,INFO/AD,SP', '-o', bcffile, dedupfile], stderr=log)
        # -I       : dont call indels in BCF
        # -B       : disable read re-alignment based on quality scores (-B reduces false SNP calls) 
        # -f <str> : input reference FASTA file name
        # -ug      : output as uncompressed BCF
        # -Q <int> : minimum base quality to be considered 
        # -t <str> : output tags: depth, allelic depth, total allelic depths, strand bias
        # -o       : output file name


        # Call variants from pileup file
        subprocess.run(['bcftools', 'call', '-m', '--ploidy', str(ploidy), '-v', '-O', 'v', '-o', vcffile, bcffile], stderr=log)
        # -m             : multiallelic caller
        # --ploidy <str> : ploidy of genome
        # -v             : only output called variants
        # -O <char>      : output type, v is uncompressed VCF
        # -o <str>       : output file name


        # Create an expression to filter called variants if they have depth below cutoff, or allelic fraction below cutoff
        filt_expression = '%MIN(DP)<' + str(min_depth_parent) + ' || %MIN(AD[1]/(AD[0]+AD[1]))<' + str(min_allele_fraction_parent)

        # Filter variants from VCF using the filter expression above
        subprocess.run(['bcftools', 'filter', '-e', filt_expression, '-o', filtvcffile, vcffile], stderr=log)
        # -e <str> : input filter expression
        # -o <str> : output file name


    # Create and write a file with genome coverage statistics
    with open(covfile, 'wt') as cov:
        subprocess.run(['bedtools', 'genomecov', '-ibam', dedupfile], stdout = cov)
        # -ibam : input BAM file




# Create a list to append with each mutant's final VCF filename, to concatenate into a format that can be uploaded to the Plasmodium database
PlasmoDB_list = []

# Go through the mutant dictionary to analyze each sample sequentially
# Create filenames for each step using the dictionary key (the mutant IDs)
for key in sample_dicty:
    logfile = key + '_log.txt'
    samfile = key + '.sam'
    bamfile = key + '.bam'
    sortbamfile = key + 's.bam'
    dedupfile = key + 'sd.bam'
    covfile = key + 'cov.txt'
    bcffile = key + '.bcf'
    vcffile = key + '.vcf'
    filtvcffile = key + '.filt.vcf'
    read1, read2 = sample_dicty[key]
    
    # Create a log file to store info for each sample from the stderr output of each command line step
    with open(logfile, 'wt') as log:
        
        # Align reads to the reference genome index (REF) using bowtie
        subprocess.run(['bowtie', '-S', '-p', str(threads), '-v', str(mismatches), '-m', str(alignments), '-5', str(trim5), '-3', str(trim3), '-X', str(max_length), 'REF', '-1', read1, '-2', read2, samfile], stderr = log)
        # -p <int> : number of threads to run                                                    
        # -v <int> : maximum number of mismatches to allow per read alignment (else discard)     
        # -m <int> : maximum number of possible alignments per read (else discard)               
        # -5 <int> : trim bases from 5' end of read                                              
        # -3 <int> : trim bases from 3' end of read                                              
        # -X <int> : maximum distance allowed between paired-end reads (outer-to-outer)          
        # -1 <str> : input paired-end file 1                                                     
        # -2 <str> : input paired-end file 2


    # Append the stderr output log from each step (i.e. number of reads aligned, percent of PCR duplicates removed, etc) to a log file for each sample
    with open(logfile, 'a') as log:

        # Convert the SAM output to a BAM file
        subprocess.run(['samtools', 'view', '-@', str(threads), '-b', '-o', bamfile, samfile], stderr=log)
        # -@ <int> : number of threads to run    
        # -b       : BAM output                  
        # -o <str> : output file name 


        # Sort the BAM file
        subprocess.run(['samtools', 'sort', '-o', sortbamfile, bamfile], stderr=log)
        # -o <str> : output filename


        # Remove PCR duplicate reads from the sorted bam file
        subprocess.run(['samtools', 'rmdup', sortbamfile, dedupfile], stderr=log)


        # Create a pileup file from the sorted and de-duplicated BAM file
        subprocess.run(['samtools', 'mpileup', '-I', '-B', '-f', referencefile, '-ug', '-Q', str(min_base_qual), '-t', 'DP,AD,INFO/AD,SP', '-o', bcffile, dedupfile], stderr=log)
        # -I       : dont call indels in BCF                                                         
        # -B       : disable read re-alignment based on quality scores (-B reduces false SNP calls) 
        # -f <str> : input reference FASTA file name                                                 
        # -ug      : output as uncompressed BCF                                                      
        # -Q <int> : minimum base quality to be considered                                           
        # -t <str> : output tags: depth, allelic depth, total allelic depths, strand bias            
        # -o       : output file name


        # Call variants from pileup file
        subprocess.run(['bcftools', 'call', '-m', '--ploidy', str(ploidy), '-v', '-O', 'v', '-o', vcffile, bcffile], stderr=log)
        # -m             : multiallelic caller                   
        # --ploidy <str> : ploidy of genome                      
        # -v             : only output called variants           
        # -O <char>      : output type - v is uncompressed VCF   
        # -o <str>       : output file name 


        # Create an expression to filter called variants if they have depth below cutoff, or allelic fraction below cutoff
        filt_expression = '%MIN(DP)<' + str(min_depth_sample) + ' || %MIN(AD[1]/(AD[0]+AD[1]))<' + str(min_allele_fraction_sample)


        # Filter variants from the VCF using the filter expression above
        subprocess.run(['bcftools', 'filter', '-e', filt_expression, '-o', filtvcffile, vcffile], stderr=log)
        # -e <str> : input filter expression            
        # -o <str> : output file name


    # Create and write a file with genome coverage statistics
    with open(covfile, 'wt') as cov:                                                                                                   
        subprocess.run(['bedtools', 'genomecov', '-ibam', dedupfile], stdout = cov)        
        # -ibam : import BAM file


    # Name output files for each parental VCF subtraction step, and annotation and filtering steps for mutant VCFs
    sub1 = key + '-1.vcf'
    sub16 = key + '-1-16.vcf'
    sub17 = key + '-1-16-17.vcf'
    subdd2 = key + '-1-16-17-dd2'
    snpeff = key + 'snpeff.vcf'
    snpeffhtml = key + 'snpeff.html'
    snpeff_filt = key + 'snpeff.filt.vcf'
    snpeff_filt_var = key + 'snpeff.filt-var.vcf'


    # Sequentially subtract all parental VCFs from each mutant VCF and save to output file using stdout
    # -a <str> : mutant VCF
    # -b <str> : parent VCF to subtract
    with open(sub1, 'wt') as outfile:
        subprocess.run(['Bedtools', 'subtract', '-header', '-a', filtvcffile, '-b', 'ssra2_1filt.vcf'], stdout = outfile)

    with open(sub16, 'wt') as outfile:
        subprocess.run(['Bedtools', 'subtract', '-header', '-a', sub1, '-b', 'ssra2_16filt.vcf'], stdout = outfile)

    with open(sub17, 'wt') as outfile:
        subprocess.run(['Bedtools', 'subtract', '-header', '-a', sub16, '-b', 'ssra2_17filt.vcf'], stdout = outfile)

    with open(subdd2, 'wt') as outfile:
        subprocess.run(['Bedtools', 'subtract', '-header', '-a', sub17, '-b', 'ssra2_dd2filt.vcf'], stdout = outfile)
        # ssra2_dd2filt.vcf was generated identically to other parent strains except that trim5=20 and trim3=30


    # Annotate the subtracted mutant VCF by calling the P. falciparum genome annotation
    with open(snpeff, 'wt') as outfile:
        subprocess.run(['SnpEff', '-Xmx6g', '-v', '-s', snpeffhtml, 'Pf3D7v91', subdd2], stdout = outfile)
        # -Xmx6g   : use 6G of memory for this step
        # -v       : verbose mode
        # -s <str> : snpeff html file


    # Create an expression to filter annotated SNP's by their functional consequence - only non-synonymous coding variants
    snpsift_expression = "(ANN[*].IMPACT has 'HIGH') || (ANN[*].IMPACT has 'MODERATE')"

    # Filter the SnpEff VCF using the expression above
    with open(snpeff_filt, 'wt') as outfile:
        subprocess.run(['SnpSift', 'filter', snpsift_expression, snpeff], stdout = outfile)


    # Create an expression to filter out genes present in a separate file (known hypervariable antigen variation genes)
    vargene_expression = "! (ANN[0].GENE in SET[0])"

    # Filter out the known hypervariable antigen variation genes (from rif+stev+EMP.txt file)
    with open(snpeff_filt_var, 'wt') as outfile:
        subprocess.run(['SnpSift', 'filter', '-s', 'rif+stev+EMP.txt', vargene_expression, snpeff_filt], stdout = outfile)

    # Append this filtered annotated VCF filename to a list, to prepare for upload to the Plasmodium database
    PlasmoDB_list.append(snpeff_filt_var)


# Create file names for concatenating all mutant sample VCFs, and for formatting each variant call as a SNP ID   
# This SNP ID file can then be input to the Plasmodium database to filter out known SNPs                         
PlasmoDB_catfile = 'ssra2_PlasmoDBcat.txt'
PlasmoDB_infile = 'ssra2_PlasmoDBin.txt'

# Concatenate the variants from each mutant's final VCF into a single file
with open(PlasmoDB_catfile, 'wt') as writefile:
    for filename in PlasmoDB_list:
        with open(filename, 'rt') as readfile:
            for line in readfile:
                if not line.startswith('#'):    
                    writefile.write(line)

# For each SNP call in this concatenated VCF, write the SNP ID to a new file in the appropriate format
# This file can then be uploaded to the Plasmodium database PlasmoDB to filter out known SNPs
with open(PlasmoDB_infile, 'wt') as writefile, open(PlasmoDB_catfile, 'rt') as readfile:
    for line in readfile:
        chrom, pos, ID, ref, alt, qual, filt, comments, _, _ = line.strip('\n').split('\t')
        writefile.write('NGS_SNP.' + chrom + '.' + pos + '\n')


