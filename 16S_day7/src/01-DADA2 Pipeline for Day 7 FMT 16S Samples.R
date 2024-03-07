# PROCESSING OF DAY 7 FMT SAMPLES USING THE DADA2 PIPELINE ---------------------
# NAME: G. Brett Moreau
# DATE: December 8, 2023

# PACKAGES ---------------------------------------------------------------------
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
packageVersion("BiocManager") # I'm using version 1.30.22

#BiocManager::install("dada2")
library(dada2)
packageVersion("dada2") # I'm using version 1.30.0

#BiocManager::install("phyloseq")
library(phyloseq)
packageVersion("phyloseq") # I'm using version 1.46.0

#install.packages("tidyverse")
library(tidyverse)
packageVersion("tidyverse") # I'm using version 2.0.0


# INTRODUCTION -----------------------------------------------------------------

# DNA was isolated from cecal contents from mice treated with antibiotics and 
# administered either FMT or PBS as a control at Day 7 post-FMT. Donor mice from
# which the fecal slurry used as FMT was derived were also examined. The V4 
# region of the 16S gene was PCR amplified using the Kozich et al 2017 primer 
# strategy and sequenced using Illumina MiSeq.




# ORGANIZING READS -------------------------------------------------------------


# In order to process the FASTQ files using DADA2, I'll first set up the 
# directory paths to find the FASTQ files.

path <- "../data/raw_reads" 
list.files(path, recursive = TRUE) # All files are present

# Pull forward and reverse reads
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE, recursive = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE, recursive = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "-"), `[`, 1)




# READ QUALITY PROFILES --------------------------------------------------------
# The first step of the DADA2 pipeline is to check the quality read profiles for 
# forward and reverse reads. I'll look at these profiles in aggregate and use 
# the profiles to determine where reads need to be trimmed.


# Set seed for reproducibility
#sample(1:10000, 1)  # It selected 9988
set.seed(9988) 

### FORWARD READ PROFILES ###
plotQualityProfile(fnFs, aggregate = TRUE)
#ggsave("../results/figures/Quality_Profile_forward_reads.png", 
#       width = 5, height = 3)


### REVERSE READ PROFILES ### 
plotQualityProfile(fnRs, aggregate = TRUE) 
#ggsave("../results/figures/Quality_Profile_reverse_reads.png",
#       width = 5, height = 3) 


# The reverse reads have pretty good quality and a slow drop over the length of 
# the read. Reverse reads look significantly worse than expected, with dips in 
# quality within the first ~75bp. I'll trim the final 10bp from the forward 
# reads and trim the reverse reads at around 160bp. I'll also trim the first 
# 25bp from the forward reads to remove primer sequences and remove the first 
# 70bp from the reverse reads to remove the poor quality regions.




# FILTER AND TRIM READS --------------------------------------------------------
# Place filtered files in filtered subdirectory
filtFs <- file.path("../data/filtered_reads", 
                    paste0(sample.names, "_F_filt.fastq.gz"))

filtRs <- file.path("../data/filtered_reads", 
                    paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen = c(240,160),
                     trimLeft = c(20,70),
                     maxN = 0, 
                     maxEE = c(2,2), 
                     truncQ = 2, 
                     rm.phix = TRUE,
                     compress = TRUE, 
                     multithread = TRUE) 

head(out) # None of these samples look like they lost a lot of reads during the 
# filtering / trimming process.




# DEREPLICATE FILTERED SEQUENCES -----------------------------------------------

# Dereplicate filtered reads to reduce computational time.
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names




# LEARN ERROR RATES ------------------------------------------------------------

# Learn error rates, which will be used for sample inference.
errF <- learnErrors(derepFs, multithread = TRUE)
errR <- learnErrors(derepRs, multithread = TRUE)

# Plot error rates and save figures.
plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

ggsave("../results/figures/Error_Rates_reverse_reads.png", 
       width = 4, height = 4)
ggsave("../results/figures/Error_Rates_forward_reads.png", 
       width = 4, height = 4)

# Error rate for both reads have the same general trend as the line of best fit,
# but there is a substantial deviation of the lines from the line of best fit.




# INFER SAMPLE COMPOSITION -----------------------------------------------------
# The DADA2 sample inference algorithm will now be applied to the filtered, 
# trimmed, and dereplicated reads. This will count the number of unique sequence 
# variants in each sample.

dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)




# MERGE PAIRED READS -----------------------------------------------------------
# Now I'll merge forward and reverse reads together to generate full sequences. 
# By default, DADA2 requires that forward and reverse reads overlap by at least 
# 12 base pairs, and there must be no mismatches along these base pairs. I'll 
# keep these parameters.

mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# Inspect the merger data.frame from the first sample
head(mergers[[1]])






# CONSTRUCT SEQUENCE TABLE -----------------------------------------------------
# Finally, I'll construct a table of all Amplicon Sequence Variants (ASVs) 
# present in our samples.
seqtab <- makeSequenceTable(mergers)

dim(seqtab) # There are 1579 unique sequences in 20 samples.

table(nchar(getSequences(seqtab))) 

plot(table(nchar(getSequences(seqtab))), 
     xlab = "read length", 
     ylab="number of reads") 

# Almost all reads have the same read length (220bp). Only 1 read had a 
# different length (232). I will remove this read.

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% seq(220, 221)]

table(nchar(getSequences(seqtab2))) 

plot(table(nchar(getSequences(seqtab2))), 
     xlab = "read length", 
     ylab="number of reads") 




# REMOVE CHIMERAS --------------------------------------------------------------
# Now I will remove chimeras (sequencing artifacts that are combinations of 
# multiple unique sequences). This new table will include all sequences after 
# chimeras have been removed.

seqtab.nochim <- removeBimeraDenovo(seqtab2, 
                                    method = "consensus", 
                                    multithread = TRUE, 
                                    verbose = TRUE)

dim(seqtab.nochim) # The original data set had 1579 total ASVs in 20 samples. 
# After removing chimeras, the new data set includes only 456 ASVs, indicating 
# that 71.1%% of all ASVs from the original data set were chimeras.

sum(seqtab.nochim)/sum(seqtab) # Looking at the percentage of reads that were 
# chimeric, this is smaller (~12%), but 12% is still fairly high for this 
# pipeline. 




# TRACK READS ACROSS PIPELINE --------------------------------------------------

# Generate table of tracked reads.
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track) 

# Format table.
track <- as.data.frame(track)
track <- plyr::rename(track, replace = c("input" = "Input", 
                                         "filtered" = "Filtered",
                                         "denoisedF" = "DenoisedF",
                                         "nonchim" = "Non-Chimera"))

track$'Percent Input (%)' <- track$`Non-Chimera`/track$Input * 100

# No samples appear to have lost a significant number of reads at any one step,
# and samples maintained ~65-75% of reads throughout the pipeline.

# I'll print the table of tracked reads for future reference.
#write.csv(track, file = "../results/tables/tracked read counts_day7.csv")




# ASSIGN TAXONOMY --------------------------------------------------------------
# Now I'll assign taxonomy to the ASV table using an assignment file from the 
# Silva reference database. For this project I'll use version 138.1 for both 
# genus and species assignment.

# Make taxonomic assignments
taxa <- assignTaxonomy(seqtab.nochim, 
                       "../data/silva_nr99_v138.1_train_set.fa.gz", 
                       multithread = TRUE)

taxa <- addSpecies(taxa, 
                   "../data/silva_species_assignment_v138.1.fa.gz") 


taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL

head(taxa.print) 

# Characterize taxonomic assignment efficiency
taxa.merged <- as.data.frame(taxa.print) 
table(taxa.merged$Phylum, useNA = "ifany") # 1 / 456 (0.2%) are NAs at Phylum level.
table(taxa.merged$Family, useNA = "ifany") # 84 / 456 (18.4%) are NAs at Family level.
table(taxa.merged$Genus, useNA = "ifany") # 254 / 456 (55.7%) are NAs at Genus level

# There is a significant proportion of ASVs that are missing Family- and Genus-
# level designations. This is not unexpected for the Genus-level but is higher
# than expected for the Family level.




# ORGANIZE METADATA AND GENERATE PHYLOSEQ OBJECT -------------------------------

# The primary goal of this experiment is to test the sequencing protocol for
# it's consistency across replicates and accuracy of taxanomic assignments. 
# Because of that, there isn't much metadata for these samples, just the 
# sample identifiers for replicate samples.

metadata <- read.csv("../data/metadata_16S_d7.csv")

metadata$Sample.ID <- as.character(metadata$Sample.ID)
metadata <- dplyr::arrange(metadata, Sample.ID)
rownames(metadata) <- metadata$Sample.ID



# I'll now use this metadata along with the otu table and taxanomic assignments
# to generate a phyloseq object, which will be used for analysis.

# Generate phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               sample_data(metadata),
               tax_table(taxa))

# Assign arbitrary identifier numbers to ASVs for simplicity.
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))




# CHARACTERISTICS OF THE DATA SET ----------------------------------------------
ps.analysis <- ps # Generate new phyloseq object for downstream analysis.

# Keep only Bacterial reads. Remove Mitochondria/Chloroplast read and taxa with 
# 0 reads in these samples. Also remove ASVs with no Family designation.
ps.analysis <- subset_taxa(ps.analysis, 
                           Kingdom == "Bacteria" & 
                             Family != "Mitochondria" &
                             Family != "Chloroplast") 

ps.analysis <- prune_taxa(taxa_sums(ps.analysis) > 0, ps.analysis)

dim(ps.analysis@otu_table) # Overall, there are 20 samples in this analysis. 
# From these 20 samples there are 368 unique ASVs.

total.reads <- sample_sums(ps.analysis)

sum(total.reads) # There are a total of 611,757 reads across the 20 samples.
mean(total.reads) # The average read number is 30,588
median(total.reads) # The median read number is 28,366
range(total.reads) # The range of reads in samples is from 20,801-53,661 reads.


# Save image of workspace for ease of access.
save.image(file = "../results/DADA2_Pipeline_R_Workspace_Day7.RData")




