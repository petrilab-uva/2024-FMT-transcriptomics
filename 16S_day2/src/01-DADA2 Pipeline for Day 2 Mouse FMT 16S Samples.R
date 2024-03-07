# DADA2 PIPELINE FOR DAY 2 MOUSE 16S SAMPLES -----------------------------------
# NAME: G. Brett Moreau
# DATE: December 1, 2023

# PACKAGES ---------------------------------------------------------------------
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("dada2")
library(dada2)
packageVersion("dada2") # I'm using version 1.30.0

#install.packages("ggplot2")
library(ggplot2)
packageVersion("ggplot2") # I'm using version 3.4.4

#install.packages("phyloseq")
library(phyloseq)
packageVersion("phyloseq") # I'm using version 1.46.0

#install.packages("tidyverse")
library(tidyverse)
packageVersion("tidyverse") # I'm using version 2.0.0

# INTRODUCTION -----------------------------------------------------------------
# The goal of this experiment was to confirm differences in the intestinal 
# microbiome in mice after fecal mirobiota transplantation (FMT). Mice were 
# treated with antibiotics, then given either FMT or PBS as a control. DNA was
# isolated from cecal contents and the 16S V4 region was amplified by PCR before
# sequening using the Illumina MiSeq.




# ORGANIZING READS -------------------------------------------------------------

# Define the path to raw reads.
path <- "../data/raw_reads"
list.files(path) # We're in the correct directory

# Organize into forward and reverse reads.
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

set.seed(298) # Set seed for reproducibility. 



#### READ QUALITY PROFILES ################################################################################
# The first step of the process is the check the quality read profiles for forward and reverse reads. I'll
# look at these profiles in aggregate and use these profiles to determine where reads need to be trimmed.

# Plot forward reads.
plotQualityProfile(fnFs, aggregate = TRUE)
ggsave("../results/figures/Quality_Profile_forward_reads.png", 
       width = 6, height = 4)

# Plot reverse reads.
plotQualityProfile(fnRs, aggregate = TRUE)
ggsave("../results/figures/Quality_Profile_reverse_reads.png", 
       width = 6, height = 4)

# Overall, the forward reads look pretty good. Reverse reads were of worse 
# quality, but this is expected. I'll trim the last 10bp of the forward reads 
# and trim at 175bp on the reverse reads, where quality begins to drop. I'll 
# also trim the first 25bp from the front of both forward and reverse reads to 
# remove primer sequences.



# Place filtered files in filtered subdirectory
filtFs <- file.path("../data/filtered_reads", 
                    paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path("../data/filtered_reads", 
                    paste0(sample.names, "_R_filt.fastq.gz"))


out <- filterAndTrim(fnFs, 
                     filtFs, 
                     fnRs, 
                     filtRs, 
                     truncLen = c(240,175), 
                     trimLeft = c(25,25),
                     maxN = 0, 
                     maxEE = c(2,2), 
                     truncQ = 2, 
                     rm.phix = TRUE,
                     compress = TRUE, 
                     multithread = TRUE) 

head(out) # None of these samples look like they lost a ton of reads.




# DEREPLICATE FILTERED SEQUENCES -----------------------------------------------

# Dereplicate filtered reads.
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names




# LEARN ERROR RATES ------------------------------------------------------------
# DADA2 uses a parametric model to estimate error rates for the sequence reads. 
# These error rates are used downstream for sample inference.

# Learn error rates.
errF <- learnErrors(derepFs, multithread = TRUE)
errR <- learnErrors(derepRs, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)
ggsave("../results/figures/Error_Rates_forward_reads.png", 
       width = 4, height = 4)

plotErrors(errR, nominalQ = TRUE)
ggsave("../results/figures/Error_Rates_reverse_reads.png", 
       width = 4, height = 4)

# Error rate models largely follow the frequency observed in each sample and are
# consistent with the expected trend lines.




# INFER SAMPLE COMPOSITION -----------------------------------------------------
# The DADA2 sample inference algorithm will now be applied to the filtered,trimmed, and dereplicated reads.
# This will count the number of unique sequence variants in each sample.

# Perform sample inference.
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

# Merge paired end sequences.
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)




# CONSTRUCT SEQUENCE TABLE -----------------------------------------------------

# Construct sequence table of all Amplicon Sequence Variants (ASVs)
seqtab <- makeSequenceTable(mergers)

dim(seqtab) # There are 188 unique sequences in 22 samples.

# Check read length.
table(nchar(getSequences(seqtab))) 

plot(table(nchar(getSequences(seqtab))), 
     xlab = "read length", 
     ylab="number of reads") 

# All reads are 215bp long.




# REMOVE CHIMERAS --------------------------------------------------------------

# Remove chimeric sequences
seqtab.nochim <- removeBimeraDenovo(seqtab, 
                                    method = "consensus", 
                                    multithread = TRUE, 
                                    verbose = TRUE)

dim(seqtab.nochim) # The original data set had 188 total ASVs in 22 samples. 
# After removing chimeras, the new data set includes 144 ASVs, indicating that 
# 24% of all ASVs from the original data set were chimeras.

sum(seqtab.nochim)/sum(seqtab) # However, when you look at the proportion of 
# sequences, chimeras account for only 0.01% of all sequence reads.




# TRACK PIPELINE READS ---------------------------------------------------------
# For quality control I'm going to track the reads over the course of the pipeline. While I'd expect there
# to be a loss of sequences at the filtering step, there should not be a significantly lower number of 
# reads at any other step of the pipeline. If there is then that step should be examined in more detail.

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
track <- cbind("Sample" = rownames(track), track)
rownames(track) <- NULL
track <- as.data.frame(track)
track$nonchim <- as.numeric(track$nonchim)
track$input <- as.numeric(track$input)
track$percent.input <- track$nonchim / track$input *100

head(track) 

# No samples lost a large proportion of reads at any step, and the percent of 
# input reads for each sample was ~85-90%.

# Print tracked read table.
write.csv(track, file = "../results/tables/tracked reads_day2.csv")




# ASSIGN TAXONOMY --------------------------------------------------------------
# Now I'll assign taxonomy to the ASV table using an assignment file from the 
# Silva reference database. I'm using v138.1.

# Add taxonomic assignments.
taxa <- assignTaxonomy(seqtab.nochim, 
                       "../data/silva_nr99_v138.1_train_set.fa.gz", 
                       multithread = TRUE)

taxa <- addSpecies(taxa, 
                   "../data/silva_species_assignment_v138.1.fa.gz") 

taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL

# Characterize level of taxonomic assignment.
taxa.merged <- as.data.frame(taxa.print) 
table(taxa.merged$Phylum, useNA = "ifany") # 0 / 144 (0.0%) are NAs at Phylum level.
table(taxa.merged$Family, useNA = "ifany") # 3 / 144  (2.1%) are NAs at Family level.
table(taxa.merged$Genus, useNA = "ifany") # 49 / 144 (34.0%) are NAs at Genus level

# Most ASVs have a high level of taxonomic assignment, with only a few missing
# assignments at the Family level.




# HAND OFF THE PHYLOSEQ --------------------------------------------------------
# I'll now be transferring the analysis from the DADA2 package to the 
# Phyloseq package for more convenient analysis

# Import metadata.
metadata <- read.csv("../data/AF_Dec2018_Metadata_subset.csv") 


# Format metadata. 
row.names(metadata) <- metadata$sample_ID
metadata$sample_ID <- NULL

# Generate phyloseq object.
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               sample_data(metadata),
               tax_table(taxa))

# Give ASVs shortened identifiers for convenience.
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

# The Phyloseq object is now ready for use.





# CHARACTERIZE PHYLOSEQ DATA SET -----------------------------------------------

ps.analysis <- ps # Generate new phyloseq object for downstream analysis.

# Keep only Bacterial reads. Remove Mitochondria/Chloroplast reads and taxa with 
# 0 reads in these samples.
ps.analysis <- subset_taxa(ps.analysis, 
                           Kingdom == "Bacteria" & 
                             Family != "Mitochondria" &
                             Family != "Chloroplast") 

ps.analysis <- prune_taxa(taxa_sums(ps.analysis) > 0, ps.analysis)

### PLATE CHARACTERISTICS ###
dim(ps.analysis@otu_table) # Overall, there are 22 samples in this analysis. 
# From these 22 samples there are 139 unique ASVs.

total.reads <- sample_sums(ps.analysis)

sum(total.reads) # There are a total of 135,398 reads across the 22 samples.
mean(total.reads) # The average read number is 6,154
median(total.reads) # The median read number is 6,416
range(total.reads) # The range of reads in samples is from 10-7,695




# EVALUATE NEGATIVE CONTROL SAMPLE ---------------------------------------------

# Organize reads per sample.
NC_evaluate <- as.matrix(sample_data(ps.analysis))
NC_evaluate <- as.data.frame(NC_evaluate)
NC_evaluate$Read_number <- sample_sums(ps.analysis)
NC_evaluate <- NC_evaluate %>%
  dplyr::arrange(Read_number) %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::mutate(Negative = ifelse(Sample == "WaterNegD", "Yes", "No"))


ggplot(data = NC_evaluate, aes(x = Read_number, 
                               y = reorder(Sample, Read_number), 
                               fill = Negative)) + 
  geom_bar(stat = "identity") +
  labs(x = "Sample ID", y = "Number of Reads", fill = "Neg. Control") +
  theme_bw() 

ggsave("../results/figures/read counts per sample.png", width = 6, height = 5)

# The negative control had almost no reads, significantly fewer than those of 
# the samples or mock community control.




# EVALUATE MOCK COMMUNITY CONTROL ----------------------------------------------

# Convert phyloseq object to realtive abundances in mock community only.
ps.norm <- transform_sample_counts(ps, function(ASV) 100 * ASV/sum(ASV))
ps.norm <- prune_samples(sample_names(ps.norm) == "ZymoPosD", ps.norm)
ps.norm <- prune_taxa(taxa_sums(ps.norm) > 0, ps.norm) 


# Add ASV information
Mock.ASVs <- as.data.frame(ps.norm@otu_table)
Mock.ASVs <- t(Mock.ASVs)
Mock.ASVs <- as.data.frame(Mock.ASVs)
Mock.ASVs <- cbind(ASV = row.names(Mock.ASVs), Mock.ASVs)
Mock.ASVs <- plyr::rename(Mock.ASVs, replace = c("ZymoPosD" = "Observed Relative Abundance (%)"))

# Add taxonomic assignments.
Mock.taxa <- as.data.frame(ps.norm@tax_table)
Mock.taxa <- cbind(ASV = row.names(Mock.taxa), Mock.taxa)
Mock.taxa <- select(Mock.taxa, ASV, Genus)
Mock.taxa <- as.data.frame(Mock.taxa)

Mock.ASV.taxa <- full_join(Mock.taxa, Mock.ASVs, by = "ASV")


# Compare to predicted abundances.
Mock.predicted <- read.csv("../data/Mock Community Expected Abundance.csv")
Mock.predicted <- plyr::rename(Mock.predicted, 
                               replace = c("Expected.Relative.Abundance...." = 
                                             "Expected Relative Abundance (%)"))

Mock.predicted.observed <- full_join(Mock.predicted, Mock.ASV.taxa, by = "Genus")
Mock.predicted.observed$ASV <- NULL

View(Mock.predicted.observed)

# All 8 predicted mock community members were taxonomically assigned by DADA2,
# indicating that taxonomic assignment performed as expected. Mock community 
# members were detected at approximately the relative abundance levels expected.
# Two additional ASVs were observed, but both were at <0.2% relative abundance.

write.csv(Mock.predicted.observed, file = "../results/tables/Mock Community Composition.csv")


# Overall, our sequencing looks to be good based on our Negative and Mock 
# community controls. I will now proceed to downstream analysis.

# Save image of workspace for ease of access.
save.image(file = "../results/DADA2_Pipeline_R_Workspace_Day2.RData")
