# Fecal Microbiota Transplantation Stimulates Type 2 and Tolerogenic Immune Responses in a Mouse Model
G. Brett Moreau, Farha Naz, William A. Petri, Jr.

## Project Summary
The following repository contains the data and code used for analysis of microbiome diversity and host gene expression during a mouse model of Fecal Microbiota Transplantation (FMT). Mice were treated with antibiotics and then treated with either FMT from healthy donor mice or PBS as a control. On Day 2 and Day 7, both 16S and bulk RNA sequencing were performed on cecal contents and cecal tissue, respectively. This work has been submitted to *Anaerobe* and is currently under review. It is accessable under DOI ____________.

## Highlights
* FMT restores microbiome alpha diversity within 7 days after antibiotic-induced dysbiosis.
* Pro-inflammatory and Type 2 immune genes are upregulated within 2 days of FMT treatment.
* Immune genes are downregulated at day 7, with increased expression of proliferation and signaling genes.
* FMT promotes increases in Type 2 and tolerogenic immune cell populations.

## Software Implementation.
This repository consists of four independent directories, one for each analysis. All source code required for each analysis is included within the `src` directory, while required metadata is included within the `data` directory. 

All data and source code included with this repository can be downloaded by cloning the git repository:
```
git clone https://github.com/petrilab-uva/2024-FMT-transcriptomics.git
```
Alternatively, you can download a [zip archive of this repository](https://github.com/petrilab-uva/2024-FMT-transcriptomics/archive/refs/heads/main.zip). 

## Dependencies
All source code for 16S analyses are written in R (Version 4.3.2). Installation of R and any dependent packages (which are outlined in each .R file) are required to perform these analyses.

Source code for RNA sequencing analyses include .R scripts as well as .sh scripts, which must be run on the command line. These .sh scripts require a Python environment to run the code. This environment was set up using [Anaconda](https://www.anaconda.com/download/), which provides the `conda` package manager. Conda virtual environments were used in these analyses to install required packages (FastQC, MultiQC, BBMap, and Kallisto) in isolation. To do this, create a new conda environment using the following code, then install the required packages within this environment.
```
conda create --name ENVIRONMENT_NAME_HERE
```

Before running the .sh script within the terminal, activate the conda environment using the following code.
```
conda activate ENVIRONMENT_NAME_HERE
```


## Setup
There are several files that are required for this analysis but not included. Instructions for downloading these files are outlined below:

1. Sequencing files have been uploaded to the NCBI Sequence Read Archive (SRA) and are available under Bioproject Accession Number PRJNA1078834. Instructions for each analysis are outlined below:
     * **16S_day2**: Download Day 2 16S FASTQ files and place them within a `raw_reads` folder in the `data` directory.
     * **16S_day7**: Download Day 7 16S FASTQ files and place them within a `raw_reads` folder in the `data` directory.
     * **RNAseq_day2**: Download Day 2 RNAseq FASTQ files and place them within a `day2_RNAseq_fastqs` folder in the `data` directory.
     * **RNAseq_day7**: Download Day 7 RNAseq FASTQ files and place them within a `day7_RNAseq_fastqs` folder in the `data` directory.
       
2. **For 16S Analysis**: Training sets are required for taxonomic assignment. I used Silva training set data that has been formatted for DADA2, which are available [here](https://zenodo.org/records/4587955). Training set and Species assignment .fa.gz files should be downloaded and placed within the `data` directory for 16S analysis.

3. **For RNAseq Analysis**: Two murine genomic data sets are required for these analyses. These files can be accessed from [Ensembl as an FTP download](https://useast.ensembl.org/info/data/ftp/index.html). I used release version 109, but more up-to-date releases may be currently available. These versions can be used, **but ensure that the same release version is used for both files!** Instructions for obtaining these files are outlined below.
     * A *Mus musculus* cDNA file containing all transcript sequences is required to generate an index file used for read mapping. I used release 109, which is available [here](https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/cdna/). This file should be placed within the `data` directory. Code for generating the .index file from this file is commented out within the .sh script.
     * A *Mus musculus* GTF gene set is required to match transcript and gene identifiers to mapped reads. I used release 109, which is available [here](https://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/). This file should be placed within the `data` directory.
  
## Contact Information
G. Brett Moreau - gbm5pn@virginia.edu

