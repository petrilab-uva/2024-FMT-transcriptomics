# QUALITY CONTROL AND PSEUDOALIGNMENT OF RNA-SEQ FASTQs USING KALLISTO
# NAME: G. Brett Moreau
# DATE: July 7, 2023


# INTRODUCTION -----------------------------------------------------------------

# This script will perform QC and map reads from FASTQs from the Day 2 mouse FMT 
# bulk RNA sequencing experiment. QC will be evaluated using FastQC, while reads
# will be pseudoaligned to the mouse genome using Kallisto. All QC information
# will be summarized into a QC report using MultiQC.

# Before running this script, raw FASTQ files need to be downloaded and placed 
# into a "day2_RNAseq_fastqs" directory within the "data" directory. This script 
# should be run using the command line in a conda enviroment containing the 
# packages listed below, which will be used for this script.

# Permissions need to be granted to allow this script to be loaded by going into
# the src directory, then using the following command (no quotes):

#"chmod +x ./01-fastq_processing_and_readMapping_day2.sh" 


# PACKAGE VERSIONS USED --------------------------------------------------------
# The following package versions were use for this analysis.

### python: Version 3.8.16
### conda: Version 23.5.0
### fastqc: Version 0.12.1
### multiqc: Version 1.14
### bbmap: Version 38.18
### kallisto: Version 0.44.0


# SETTING UP THE ENVIRONMENT ---------------------------------------------------

# Change to the project directory (change as necessary)
cd ..

# Make output directory for fastqc and kallisto report files
mkdir ./results/fastqc
mkdir ./results/fastqc/raw_reads
mkdir ./results/fastqc/trimmed_reads
mkdir ./results/kallisto
mkdir ./results/kallisto/logs


# QUALITY CONTROL OF PROCESSED FASTQS USING FASTQC -----------------------------

# Run FastQC on processed FASTQ files.
for i in ./data/day2_RNAseq_fastqs/*fq.gz
do
	fastqc -o ./results/fastqc/raw_reads $i -t 16
done

# Summarize FastQC results using MultiQC
multiqc ./results \
-d -n report_day2_raw_reads.html \
-o ./results


# ADAPTER SEQUENCE TRIMMING USING BBDUK ----------------------------------------

# Trim adapter sequences using BBDuk
for i in $(basename ./data/day2_RNAseq_fastqs/*.fq.gz | awk '{print substr( $0, 0, 5)}' | uniq)
do 
	bbduk.sh \
	in1=./data/day2_RNAseq_fastqs/"$i"_1.fq.gz \
	in2=./data/day2_RNAseq_fastqs/"$i"_2.fq.gz \
	out1=./data/day2_RNAseq_fastqs/trimmed_reads/"$i"_1_trimmed.fq.gz \
	out2=./data/day2_RNAseq_fastqs/trimmed_reads/"$i"_2_trimmed.fq.gz \
	ref=./data/adapters.fa \
	ktrim=r \
	mink=11 \
	hdist=1 \
	tpe tbo \

done 2> ./results/bbduk_output.log


# Run FastQC on trimmed FASTQ files
for i in ./data/day2_RNAseq_fastqs/trimmed_reads/*fq.gz
do
	fastqc -o ./results/fastqc/trimmed_reads $i -t 16
done


# READ MAPPING WITH KALLISTO ---------------------------------------------------

# Build cDNA index from reference FASTA file for aligning reads. 
#kallisto index -i ./data/Mus_musculus.GRCm39.cdna.all.index \
#./data/Mus_musculus.GRCm39.cdna.all.fa

# NOTE: I'm using release 109 of the Mus musculus cDNA FASTA file, which can be
# accessed from the ENSEMBL depository: 
# (https://useast.ensembl.org/info/data/ftp/index.html).


for i in $(basename ./data/day2_RNAseq_fastqs/trimmed_reads/*.fq.gz | awk '{print substr( $0, 0, 5)}' | uniq)
do
	kallisto quant -i ./data/Mus_musculus.GRCm39.cdna.all.index \
	-o ./results/kallisto/"$i" -t 16 \
	./data/day2_RNAseq_fastqs/trimmed_reads/"$i"_1_trimmed.fq.gz \
	./data/day2_RNAseq_fastqs/trimmed_reads/"$i"_2_trimmed.fq.gz \
	&> ./results/kallisto/logs/"$i".log
	echo $i "finished"
done	


# SUMMARIZE REPORTS USING MULTIQC ----------------------------------------------

# Summarize using MultiQC
multiqc ./results \
--ignore ./results/fastqc/raw_reads \
-d -n report_day2_trimmed_reads_and_mapping.html \
-o ./results

echo "Finished"
