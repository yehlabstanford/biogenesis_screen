# biogenesis_screen

Files included in repository:  
ssra_pipeline_final.py  - Custom script written for NGS analysis  
Parent_IDs.csv          - Example of ID file for parental samples and their corresponding fastq file names  
Sample_IDs.csv          - Example of ID file for mutant samples and their corresponding fastq file names  
rif+stev+EMP.txt        - Contains Pf gene ID’s for hyper-variable gene families  


This pipeline inputs parent and mutant sample ID files with their paired-end FASTQ file names, as well as a FASTA reference genome file. All of the files above, each of the fastq files listed in each ID.csv file, and the reference genome (Reference_file.fasta) are to be placed in a directory, and the script can be run via the command line from that directory using:  

> python3 ssra_pipeline_final.py Parent_IDs.csv Sample_IDs.csv Reference_file.fasta   

For each sample, the pipeline aligns the fastq reads to the reference genome, removes PCR duplicates, creates a pileup, calls variants, and filters these variants. It then subtracts the parental variants from the mutant sample variants to provide a final list of mutant SNV’s. These SNV’s are then annotated, and SNV’s in hyper-variable gene families are filtered out. Finally, the pipeline concatenates filtered SNV’s from each mutant sample and generates a file with each SNP ID in a format which can be uploaded to the PlasmoDB SNP ID search page for further analysis (http://plasmodb.org/plasmo/showQuestion.do?questionFullName=SnpQuestions.NgsSnpBySourceId).   

This pipeline creates intermediate files for each sample at each analysis step, in the initial directory. It also generates a log file containing all command line errors and statistics for each sample analyzed, and generates a genome coverage file. All parameters and filtering requirements are annotated within the script.   

Requirements (version used):  
Python3 (3.6.4)  
Bowtie (1.1.2)  
Samtools (1.3.1)  
Bcftools (1.3.1)  
Bedtools (2.26.0)  
SnpEff (4.2)  
SnpSift (4.2)  
