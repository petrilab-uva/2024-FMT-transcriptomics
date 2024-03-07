# IMPORT OF MAPPED READS AND DESEQ2 PIPELINE -----------------------------------
# NAME: G. Brett Moreau
# DATE: July 26, 2023

# PACKAGES ---------------------------------------------------------------------
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
packageVersion("BiocManager") # I'm using version 1.30.22

#install.packages("readxl")
library(readxl)
packageVersion("readxl") # Using version 1.4.3

#install.packages("tidyverse")
library(tidyverse)
packageVersion("tidyverse") # Using version 2.0.0

#BiocManager::install("rhdf5")
library(rhdf5)
packageVersion("rhdf5") # Using version 2.46.1

#BiocManager::install("tximport")
library(tximport)
packageVersion("tximport") # Using version 1.30.0

#BiocManager::install("ensembldb")
library(ensembldb)
packageVersion("ensembldb") # Using version 2.26.0

#BiocManager::install("DESeq2")
library(DESeq2)
packageVersion("DESeq2") # Using version 1.42.0

#BiocManager::install("plyranges")
library(plyranges)
packageVersion("plyranges") # Using version 1.22.0

#install.packages("RColorBrewer")
library(RColorBrewer)
packageVersion("RColorBrewer") # Using version 1.1.3

#install.packages("pheatmap")
library(pheatmap)
packageVersion("pheatmap") # Using version 1.0.12


# INTRODUCTION -----------------------------------------------------------------

# Raw bulk RNAseq reads have gone through QC, trimming, and mapping to the 
# murine genome using Kallisto. I'll now import these data into R and use 
# DESeq2 to quantify and normalize counts per gene for each sample.




# IMPORT KALLISTO DATA INTO R  -------------------------------------------------

# Import sample data
metadata <- read.csv("../data/day7_RNAseq_metadata.csv")


# Generate file path for Kallisto files
path <- file.path("../results/kallisto", 
                  metadata$Novogene.Sample.ID, 
                  "abundance.tsv")

file.exists(path) # All files present
names(path) <- metadata$Novogene.Sample.ID


# Now I need to make a tx2gene file, which will allow us to match transcript and 
# gene identifiers during import of Kallisto abundance information.Because I 
# made my reference index from the Mus musculis cDNA FASTA on the Ensembl 
# website, I need to use a matching transcript database for matching gene 
# identifiers during import. I will therefore use the GTF from the Ensemble FTP 
# site (VERSION 109, MATCHING THE INDEX VERSION) for generating the transcript 
# database.

# Make Tx2gene file for Tximport
TxDB<- GenomicFeatures::makeTxDbFromGFF("../data/Mus_musculus.GRCm39.109.gtf.gz",
                                        format = "gtf")
k <- keys(TxDB, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(TxDB, k, "GENEID", "TXNAME")

# Import Kallisto data. Will also estimate counts using length-adjusted 
# abundance estimates.
txi.kallisto <- tximport(path, 
                         type = "kallisto", 
                         tx2gene = tx2gene, 
                         countsFromAbundance = "lengthScaledTPM",
                         ignoreTxVersion = TRUE)

View(txi.kallisto$counts) # Scaled counts with gene identifiers as row names.

# There are two methods of data normalization that are recommended by the makers
# of tximport before handing off to differential gene expression pipelines. The
# method used here scales counts using the average transcript length over 
# samples and the library size, then hands these adjusted counts (in the 
# 'counts' table) directly to the DGE pipeline of choice. I'll use DESeq2 for
# differential gene expression analysis.




# HAND OFF COUNTS TO DESEQ2 ----------------------------------------------------

# Format metadata for DESeq2
rownames(metadata) <- metadata$Novogene.Sample.ID
all(rownames(metadata) == colnames(txi.kallisto$counts)) # Identical names
metadata$Treatment[metadata$Treatment == "PBS Control"] <- "PBS_Control"
metadata$Treatment[metadata$Treatment == "FMT Recipient"] <- "FMT_Recipient"


metadata$Treatment <- factor(metadata$Treatment, 
                             levels = c("PBS_Control", "FMT_Recipient"))

# Make DESeq2 data set
dds <- DESeqDataSetFromTximport(txi = txi.kallisto, 
                                colData = metadata, 
                                design = ~ Treatment)

# Filter low abundance genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Generate DESeq object for differential gene expression
dds <- DESeq(dds)




# EVALUATE GROUP DIFFERENCES BY PRINCIPAL COMPONENT ANALYSIS -------------------

# Log transform and organize count data for PCA
rld <- rlog(dds, blind = TRUE)
plotpca <- plotPCA(rld, intgroup = c("Treatment"), returnData = TRUE)
percentVar <- round(100 * attr(plotpca, "percentVar"))

plotpca$Group[plotpca$Treatment == "PBS_Control"] <- "PBS Control"
plotpca$Group[plotpca$Treatment == "FMT_Recipient"] <- "FMT Recipient"
plotpca$Group <- factor(plotpca$Group, 
                        levels = c("PBS Control", "FMT Recipient"))


# Plot PCA
ggplot(data = plotpca, aes(x = PC1, y = PC2
                           #color = metadata$Treatment,
                           #label = metadata$name
)) + 
  geom_point(aes(color = Group), size = 6) +
  #geom_text(hjust=0, vjust=0) +
  labs(color = "Group") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values = c("#808080", "#ff8080")) +
  theme_bw() +
  #scale_color_discrete(name = "Treatment", labels = c("Vehicle Control", "FMT Recipient")) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 30), 
        axis.title.y = element_text(size = 30),
        legend.title = element_blank(),
        legend.text = element_text(size = 30),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        aspect.ratio = 1)

# Save PCA plot.
ggsave("../results/figures/PCA plot_day7.png", 
       dpi = 300)




### HEATMAP OF WITHIN-SAMPLE VARIABILITY ###

sampleDists <- dist(t(assay(rld)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$Novogene.Sample.ID, " (", rld$Treatment, ")")
colnames(sampleDistMatrix) <- paste(rld$Novogene.Sample.ID, " (", rld$Treatment, ")")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

# Visualize heatmap (remove comments to save plot.)
#png("../results/figures/heatmap_sample distances_day7.png", 
#    width = 750, height = 750)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)
#dev.off()


# Generate results files for each comparison.
res_FMTvsVEHICLE <- results(dds, contrast = c("Treatment",
                                              "FMT_Recipient", 
                                              "PBS_Control"))




# QUALITY CONTROL FOR RNA SEQUENCING DATA --------------------------------------

### DISPERSION ANALYSIS ###

# Plot Dispersion Estimates (remove comments to save plots)
#png("../results/figures/dispersion estimates_day7.png", 
#    width = 750, height = 750)
plotDispEsts(dds)
#dev.off()


# The final dispersion estimates are shrunk compared to the gene-wise estimates.
# This is the expected result from this analysis.


### MA-PLOT ###

# Generate and save MA-plots (remove comments to save plots)
#png("../results/figures/MA_FMT_vs_control_day7.png", 
#    width = 750, height = 750)
DESeq2::plotMA(res_FMTvsVEHICLE)
#dev.off()


# Save only results files for downstream analysis
rm(list = setdiff(ls(), c("res_FMTvsVEHICLE", "rld", "dds", "metadata")))

save.image(file = "../results/DESeq2 results_day7.RData")


