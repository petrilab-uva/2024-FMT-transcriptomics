# DIFFERENTIALLY EXPRESSED GENE AND PATHWAY ANALYSIS ---------------------------
# NAME: G. Brett Moreau
# DATE: July 25, 2023

# PACKAGES ---------------------------------------------------------------------
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
packageVersion("BiocManager") # I'm using version 1.30.22

#install.packages("tidyverse")
library(tidyverse)
packageVersion("tidyverse") # Using version 2.0.0

#install.packages("msigdbr")
library(msigdbr)
packageVersion("msigdbr") # Using version 7.5.1

#BiocManager::install("fgsea")
library(fgsea)
packageVersion("fgsea") # Using version 1.28.0

#BiocManager::install("DESeq2")
library(DESeq2)
packageVersion("DESeq2") # Using version 1.42.0

#BiocManager::install("plyranges")
library(plyranges)
packageVersion("plyranges") # Using version 1.22.0

#install.packages("pheatmap")
library(pheatmap)
packageVersion("pheatmap") # Using version 1.0.12

#install.packages("ggrepel")
library(ggrepel)
packageVersion("ggrepel") # Using version 0.9.5

# INTRODUCTION -----------------------------------------------------------------

# Bulk RNAseq counts have been normalized and quantified for each gene. I will
# now examine differentially expressed genes (DEGs) for each comparison along
# with enrichment analysis to identify biological pathways associated with 
# these changes.




# GENERATING TABLES OF DIFFERENTIALLY EXPRESSED GENES --------------------------

# Import DESeq2 results files
load("../results/DESeq2 results_day7.RData")

# Import gene name information for appending results files.
gene_names <- read_gff("../data/Mus_musculus.GRCm39.109.gtf.gz") 
gene_names <- plyranges::select(gene_names, gene_id, gene_name)
gene_names <- as.data.frame(gene_names)
gene_names <- gene_names %>%
  dplyr::select(gene_name, gene_id) %>%
  dplyr::distinct()


# Add gene names to results files
res_FMTvsVEHICLE_table <- as.data.frame(res_FMTvsVEHICLE)
res_FMTvsVEHICLE_table$gene_id <- rownames(res_FMTvsVEHICLE)
res_FMTvsVEHICLE_table <- res_FMTvsVEHICLE_table %>%
  dplyr::left_join(., gene_names, by = "gene_id") %>%
  dplyr::relocate(gene_id, gene_name, .before = baseMean)


# Save results file.
write.csv(res_FMTvsVEHICLE_table, 
          file = "..//results/tables/FMT_vs_control_DEG_day7.csv", 
          row.names = FALSE)



# Overview of differentially expressed genes

# Count number of DEGs overall
sum(res_FMTvsVEHICLE_table$padj < 0.05, na.rm = TRUE) # 233 DEGs

# Count Number of DEGs per group
res_FMTvsVEHICLE_table %>%
  dplyr::filter(padj < 0.05) %>% 
  dplyr::filter(log2FoldChange > 0) %>%
  dplyr::select(gene_name) # 103 DEGs upregulated in FMT


res_FMTvsVEHICLE_table %>%
  dplyr::filter(padj < 0.05) %>% 
  dplyr::filter(log2FoldChange < 0) %>%
  dplyr::select(gene_name) # 130 DEGs upregulated in control




# VOLCANO PLOT OF DIFFERENTIALLY EXPRESSED GENES -------------------------------

# Generate column for colors
res_FMTvsVEHICLE_volcano <- res_FMTvsVEHICLE_table %>%
  mutate(Significance = ifelse(padj > 0.05, "Not Significant", 
                               ifelse(padj < 0.05 & log2FoldChange >= 1, 
                                      "Upregulated", 
                                      ifelse(padj < 0.05 & log2FoldChange <= -1,
                                             "Downregulated", "No FC")))) %>%
  mutate(Significance = factor(Significance, levels = c("Not Significant",
                                                        "No FC",
                                                        "Upregulated",
                                                        "Downregulated")))


# Plot volcano plot
ggplot(data = res_FMTvsVEHICLE_volcano, aes(x = log2FoldChange,
                                          y = -log10(padj), 
                                          color = Significance)) +
  geom_point(size = 5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = -1, linetype = "dashed") +
  labs(x = "Log2(Fold Change)", y = "Log10(p-value)") +
  scale_color_manual(values = c("gray", "gray", "#00BFC4", "#F8766D")) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_text(size = 40), 
        axis.text.y = element_text(size = 40), 
        axis.title.x = element_text(size = 40), 
        axis.title.y = element_text(size = 40),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 38))

# Save volcano plot.
ggsave("../results/figures/volcano_FMT_vs_control_day7.png", 
       height = 9,
       dpi = 300)




# HEATMAPS OF TOP DIFFERENTIALLY EXPRESSED GENES -------------------------------

### FMT vs CONTROL ###
# Select 50 most differentially expressed genes.
heatmap_FMTvsVEHICLE <- 
  assay(rld)[ head(order(res_FMTvsVEHICLE$padj), 50), ]

# Organize data set to remove unnecessary samples and add gene names.
heatmap_FMTvsVEHICLE <- as.data.frame(heatmap_FMTvsVEHICLE)
heatmap_FMTvsVEHICLE$gene_id <- rownames(heatmap_FMTvsVEHICLE)
heatmap_FMTvsVEHICLE <- heatmap_FMTvsVEHICLE %>%
  dplyr::left_join(., gene_names, by = "gene_id")

# Convert gene names to row names
rownames(heatmap_FMTvsVEHICLE) <- heatmap_FMTvsVEHICLE$gene_name
heatmap_FMTvsVEHICLE <- dplyr::select(heatmap_FMTvsVEHICLE, -gene_id, -gene_name)

# Convert count data to change in expression across the average of all samples.
heatmap_FMTvsVEHICLE <- heatmap_FMTvsVEHICLE - rowMeans(heatmap_FMTvsVEHICLE)
heatmap_FMTvsVEHICLE <- as.matrix(heatmap_FMTvsVEHICLE)

# Format metadata to add to plot
key <- metadata %>% 
  dplyr::select(Treatment)

# Plot heatmap (remove comments to save plot.)
#png("../results/figures/heatmap_FMT_vs_control_day72.png", 
#    width = 7, height = 10, units = "in", res = 300)
pheatmap(heatmap_FMTvsVEHICLE, 
         annotation_col = key)
#dev.off()




# ORGANIZATION FOR GENE SET ENRICHMENT ANALYSIS --------------------------------

# GSEA will be performed with the fgsea package. I will look for enrichment 
# using the Hallmark gene set from the MSigDB Collection. In order to perform
# GSEA, I need to have a ranked list of genes and the gene set to compare it to.

# I'll start by importing the Hallmark gene set, then organize individual ranked
# gene lists.

# Set seed for reproducibility
#sample(1:10000, 1) # 7625
set.seed(7625)

# Load Hallmark pathways for mouse
h_gene_set <- msigdbr(species = "Mus musculus", category = "H")

# Organize Hallmark pathway set for GSEA
h_gene_set_list <- split(x = h_gene_set$gene_symbol, 
                         f = h_gene_set$gs_name)

View(h_gene_set_list) # List with genes for each gene set name.


# Organize ranked gene list.

# Check for missing Wald statistic values
table(is.na(res_FMTvsVEHICLE_table$stat)) # no missing values

# Check for duplicate gene names
table(duplicated(res_FMTvsVEHICLE_table$stat)) # 21 duplicates

# I'll remove all genes with NA values in the Wald statistic, then average 
# together all ENSEMBL IDs that map to the same gene name. 

# Generate ranked gene list.
ranked_FMTvsVEHICLE <- res_FMTvsVEHICLE_table %>%
  dplyr::select(gene_name, stat) %>%
  na.omit() %>%
  group_by(gene_name) %>%
  summarize(stat = mean(stat)) %>%
  deframe()

head(ranked_FMTvsVEHICLE, 20)




# PERFORM GENE SET ENRICHMENT ANALYSIS -----------------------------------------

# Run GSEA
fgsea_FMTvsVEHICLE <- fgsea(pathways = h_gene_set_list,
                            stats = ranked_FMTvsVEHICLE,
                            eps = 0)


# Collapse leadingEdge column from list for table printing
fgsea_FMTvsVEHICLE_table <- fgsea_FMTvsVEHICLE %>%  
  rowwise() %>%
  mutate(leadingEdge = paste(leadingEdge, collapse=',')) %>%
  ungroup()

# Count number of commas in leading edge column as proxy for number of genes.
fgsea_FMTvsVEHICLE_table$enriched_genes <- 
  lengths(gregexpr(",", fgsea_FMTvsVEHICLE_table$leadingEdge)) + 1

# Calculate the percent of pathway genes enriched for each pathway
fgsea_FMTvsVEHICLE_table$percent_coverage <- 
  fgsea_FMTvsVEHICLE_table$enriched_genes / 
  fgsea_FMTvsVEHICLE_table$size * 100

# Save GSEA output table.
write.csv(fgsea_FMTvsVEHICLE_table, 
          file = "../results/tables/GSEA_FMT_vs_control_day7.csv", 
          row.names = FALSE)


# Plot all significantly enriched pathways
fgsea_plot <- fgsea_FMTvsVEHICLE_table %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::mutate(Enrichment = ifelse(NES > 0, "FMT Recipient", "PBS Control")) %>%
  dplyr::mutate(Enrichment = factor(Enrichment, 
                                    levels = c("PBS Control", "FMT Recipient"))) 

# Rename gene sets for plot
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_ADIPOGENESIS"] <- "Adipogenesis"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_ALLOGRAFT_REJECTION"] <- "Allograft Rejection"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_ANGIOGENESIS"] <- "Angiogenesis"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_APICAL_JUNCTION"] <- "Apical Junction"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_APICAL_SURFACE"] <- "Apical Surface"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_APOPTOSIS"] <- "Apoptosis"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_BILE_ACID_METABOLISM"] <- "Bile Acid Metabolism"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_CHOLESTEROL_HOMEOSTASIS"] <- "Cholesterol Homeostais"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_COAGULATION"] <- "Coagulation"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_COMPLEMENT"] <- "Complement"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_E2F_TARGETS"] <- "E2F Targets"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"] <- "Epithelial-Mesenchymal Transition"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_ESTROGEN_RESPONSE_EARLY"] <- "Estrogen Response: Early"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_ESTROGEN_RESPONSE_LATE"] <- "Estrogen Response: Late"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_FATTY_ACID_METABOLISM"] <- "Fatty Acid Metabolism"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_GLYCOLYSIS"] <- "Glycolysis"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_HEDGEHOG_SIGNALING"] <- "Hedgehog Signaling"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_HEME_METABOLISM"] <- "Heme Metabolism"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_IL2_STAT5_SIGNALING"] <- "IL-2/STAT5 Signaling"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_IL6_JAK_STAT3_SIGNALING"] <- "IL-6/JAK-STAT3 Signaling"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_INFLAMMATORY_RESPONSE"] <- "Inflammatory Response"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_INTERFERON_ALPHA_RESPONSE"] <- "IFNa Response"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE"] <- "IFNg Response"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_KRAS_SIGNALING_UP"] <- "KRAS Signaling: Up"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_MTORC1_SIGNALING"] <- "mTORC1 Signaling"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_MYOGENESIS"] <- "Myogenesis"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION"] <- "Oxidative Phosphorylation"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_P53_PATHWAY"] <- "p53 Pathway"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_PEROXISOME"] <- "Peroxisome"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_TGF_BETA_SIGNALING"] <- "TGFb Signaling"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"] <- "TGFa Signaling via NFkB"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"] <- "Unfolded Protein Response"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_UV_RESPONSE_DN"] <- "UV Response: Down"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_XENOBIOTIC_METABOLISM"] <- "Xenobiotic Metabolism"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_ANDROGEN_RESPONSE"] <- "Androgen Response"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_HYPOXIA"] <- "Hypoxia"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_KRAS_SIGNALING_DN"] <- "KRAS Signaling: Down"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_PI3K_AKT_MTOR_SIGNALING"] <- "PI3K/AKT/mTOR Signaling"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_UV_RESPONSE_UP"] <- "UV Response: Up"


# Visualize bubble plot of enriched pathways.
ggplot(fgsea_plot, aes(x = NES, y = reorder(pathway, NES), size = enriched_genes, color = Enrichment)) +
  geom_point() +
  scale_size(range = c(1, 10)) +
  labs(x = "NES", y = NULL,
       size = "# Genes") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.direction = "vertical", 
        legend.justification = "left", 
        axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 13), 
        axis.title.x = element_text(size = 13), 
        axis.title.y = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10))

# Save bubble plot.
ggsave("../results/figures/GSEA_FMT_vs_control_Hallmark_allsig_bubble.png", 
       dpi = 300)


# GSEA Plot as Horizontal Bar Plot
ggplot(fgsea_plot, aes (x = reorder(pathway, NES), y = NES, 
                        fill = Enrichment)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, y = "NES", color = "Treatment") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 24.5, linetype = "dashed") +
  scale_fill_manual(values = c("#808080", "#ff8080")) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Save plot.
ggsave("../results/figures/GSEA_FMTvsVEHICLE_Top10_Enrichment_Day7_horizontal.png", 
       dpi = 300, height = 8, width = 11)




# GENERATE BOXPLOTS OF SELECT LEADING EDGE GENES -------------------------------

# Organize normalized count data
count_table <- counts(dds, normalized = TRUE)

# Organize data and add groupings
count_table <- count_table %>%
  as.data.frame() %>%
  dplyr::mutate(gene_id = rownames(.)) %>%
  tidyr::pivot_longer(cols = !gene_id, 
                      names_to = "Novogene.Sample.ID", 
                      values_to = "Normalized_Counts") %>%
  dplyr::left_join(., gene_names, by = "gene_id") %>%
  dplyr::left_join(., metadata, by = "Novogene.Sample.ID") %>%
  dplyr::mutate(Treatment = factor(Treatment, 
                                   levels = c("PBS_Control", "FMT_Recipient")))


## CSF3R ##

# Visualize counts.
count_table %>%
  dplyr::filter(gene_name == "Csf3r") %>%
  ggplot(aes(x = Treatment, 
             y = Normalized_Counts, 
             color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Csf3r", y = "Normalized Counts") +
  scale_color_manual(values = c("#808080", "#ff8080")) +
  theme_classic() +
  theme(legend.position = "none", 
        #plot.title = element_text(hjust = 0.5, size = 20), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# Save boxplot.
ggsave("../results/figures/final_plots/Csf3r.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics.
count_table %>%
  dplyr::filter(gene_name == "Csf3r") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.024



## PTGS2 ##

# Visualize counts.
count_table %>%
  dplyr::filter(gene_name == "Ptgs2") %>%
  ggplot(aes(x = Treatment, 
             y = Normalized_Counts, 
             color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Ptgs2", y = NULL) +
  scale_color_manual(values = c("#808080", "#ff8080")) +
  theme_classic() +
  theme(legend.position = "none", 
        #plot.title = element_text(hjust = 0.5, size = 20), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

# Save boxplot.
ggsave("../results/figures/final_plots/Ptgs2.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics.
count_table %>%
  dplyr::filter(gene_name == "Ptgs2") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.028



## PDE4B ##

# Visualize counts.
count_table %>%
  dplyr::filter(gene_name == "Pde4b") %>%
  ggplot(aes(x = Treatment, 
             y = Normalized_Counts, 
             color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Pde4b", y = NULL) +
  scale_color_manual(values = c("#808080", "#ff8080")) +
  theme_classic() +
  theme(legend.position = "none", 
        #plot.title = element_text(hjust = 0.5, size = 20), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

# Save boxpot.
ggsave("../results/figures/final_plots/Pde4b.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics.
count_table %>%
  dplyr::filter(gene_name == "Pde4b") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.0061



## CXCL9 ##

# Visualize counts.
count_table %>%
  dplyr::filter(gene_name == "Cxcl9") %>%
  ggplot(aes(x = Treatment, 
             y = Normalized_Counts, 
             color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Cxcl9", y = NULL) +
  scale_color_manual(values = c("#808080", "#ff8080")) +
  theme_classic() +
  theme(legend.position = "none", 
        #plot.title = element_text(hjust = 0.5, size = 20), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# Save boxplot.
ggsave("../results/figures/final_plots/Cxcl9.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics.
count_table %>%
  dplyr::filter(gene_name == "Cxcl9") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.0565



## EREG ##

# Visualize counts
count_table %>%
  dplyr::filter(gene_name == "Ereg") %>%
  ggplot(aes(x = Treatment, 
             y = Normalized_Counts, 
             color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Ereg", y = NULL) +
  scale_color_manual(values = c("#808080", "#ff8080")) +
  theme_classic() +
  theme(legend.position = "none", 
        #plot.title = element_text(hjust = 0.5, size = 20), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# Save boxplot.
ggsave("../results/figures/final_plots/Ereg.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics.
count_table %>%
  dplyr::filter(gene_name == "Ereg") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.024



## AREG ##

# Visualize counts.
count_table %>%
  dplyr::filter(gene_name == "Areg") %>%
  ggplot(aes(x = Treatment, 
             y = Normalized_Counts, 
             color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Areg", y = NULL) +
  scale_color_manual(values = c("#808080", "#ff8080")) +
  theme_classic() +
  theme(legend.position = "none", 
        #plot.title = element_text(hjust = 0.5, size = 20), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

# Save boxplots.
ggsave("../results/figures/final_plots/Areg.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics.
count_table %>%
  dplyr::filter(gene_name == "Areg") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.047


