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

# I have performed differential gene expression analysis and identified Hallmark
# gene sets associated with each group. There were limited gene sets upregulated
# by FMT, so I want to see if more are identified using a list of gene sets that
# is more granular.




# LOAD AND ORGANIZE DEG DATA  --------------------------------------------------

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




# GENERATE GENE SET FOR GO GSEA ------------------------------------------------

# I'll repeat GSEA using the GO Biological Processes gene set,  using Go terms, 
# which should give a broader picture of biological pathways of interest.


# Make gene set with GO terms for mouse
### BIOLOGICAL PROCESSES ###
GO_BP_set <- msigdbr(species = "Mus musculus", 
                     category = "C5", 
                     subcategory = "GO:BP")

GO_BP_set_list <- split(x = GO_BP_set$gene_symbol, 
                        f = GO_BP_set$gs_name)


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


# Run GSEA
fgsea_FMTvsVEHICLE_GO_BP <- fgsea(pathways = GO_BP_set_list, 
                                  stats = ranked_FMTvsVEHICLE, 
                                  eps = 0)

# Collapse leadingEdge column from list for table printing
fgsea_FMTvsVEHICLE_GO_BP_table <- fgsea_FMTvsVEHICLE_GO_BP %>%  
  rowwise() %>%
  mutate(leadingEdge = paste(leadingEdge, collapse=',')) %>%
  ungroup()

# Save GSEA output
write.csv(fgsea_FMTvsVEHICLE_GO_BP_table, 
          file = "../results/tables/GSEA_GO_BP_FMT_vs_control_day7.csv", 
          row.names = FALSE)


# Count number of commas in leading edge column as proxy for number of genes.
fgsea_FMTvsVEHICLE_GO_BP_table$enriched_genes <- 
  lengths(gregexpr(",", fgsea_FMTvsVEHICLE_GO_BP_table$leadingEdge)) + 1

# Calculate the percent of pathway genes enriched for each pathway
fgsea_FMTvsVEHICLE_GO_BP_table$percent_coverage <- 
  fgsea_FMTvsVEHICLE_GO_BP_table$enriched_genes / 
  fgsea_FMTvsVEHICLE_GO_BP_table$size * 100



# Plot top 10 upregulated gene sets
fgsea_GO_plot <- fgsea_FMTvsVEHICLE_GO_BP_table %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(desc(NES)) %>%
  dplyr::slice_head(n = 10)

# Rename gene sets for plot
fgsea_GO_plot$pathway[fgsea_GO_plot$pathway == "GOBP_CYTOPLASMIC_TRANSLATION"] <- "Cytoplasmic Translation"
fgsea_GO_plot$pathway[fgsea_GO_plot$pathway == "GOBP_CEREBELLAR_CORTEX_MORPHOGENESIS"] <- "Cerebellar Cortex Morphogenesis"
fgsea_GO_plot$pathway[fgsea_GO_plot$pathway == "GOBP_HINDBRAIN_MORPHOGENESIS"] <- "Hindbrain Morphogenesis"
fgsea_GO_plot$pathway[fgsea_GO_plot$pathway == "GOBP_CEREBELLAR_CORTEX_DEVELOPMENT"] <- "Cerebellar Cortex Development"
fgsea_GO_plot$pathway[fgsea_GO_plot$pathway == "GOBP_FORELIMB_MORPHOGENESIS"] <- "Forelimb Morphogenesis"
fgsea_GO_plot$pathway[fgsea_GO_plot$pathway == "GOBP_CELL_DIFFERENTIATION_IN_HINDBRAIN"] <- "Cell Differentiation in Hindbrain"
fgsea_GO_plot$pathway[fgsea_GO_plot$pathway == "GOBP_POSITIVE_REGULATION_OF_SYNAPSE_ASSEMBLY"] <- "Positive Regulation of Synapse Assembly"
fgsea_GO_plot$pathway[fgsea_GO_plot$pathway == "GOBP_REGULATION_OF_SYNAPSE_ASSEMBLY"] <- "Regulation of Synapse Assembly"
fgsea_GO_plot$pathway[fgsea_GO_plot$pathway == "GOBP_TRNA_PROCESSING"] <- "tRNA Processing"
fgsea_GO_plot$pathway[fgsea_GO_plot$pathway == "GOBP_NEUROPEPTIDE_SIGNALING_PATHWAY"] <- "Neuropeptide Signaling Pathway"


# Plot top upregulated gene sets as bar plot.
ggplot(fgsea_GO_plot, aes (x = NES, y= reorder(pathway, NES), 
                          fill = Enrichment)) +
  geom_bar(stat = "identity", fill = "#ec8783") +
  labs(x = "NES", y = NULL, color = "Treatment") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Save bar graph.
ggsave("../results/figures/GSEA_FMTvsVEHICLE_Top10_Enrichment_Day7_GO-BP.png", 
       dpi = 300, height = 10, width = 8)




# UNIVARIATE ANALYSIS OF SELECT LEADING EDGE GENES -----------------------------

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


### INTESTINAL HOMEOSTASIS ###

### ZBTB16 ###

# Visualize counts.
count_table %>%
  dplyr::filter(gene_name == "Zbtb16") %>%
  ggplot(aes(x = Treatment, 
             y = Normalized_Counts, 
             color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Zbtb16", y = "Normalized Counts") +
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
ggsave("../results/figures/final_plots/Zbtb16.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics.
count_table %>%
  dplyr::filter(gene_name == "Zbtb16") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.00105



### ALDH1A2 ###

# Visualize counts.
count_table %>%
  dplyr::filter(gene_name == "Aldh1a2") %>%
  ggplot(aes(x = Treatment, 
             y = Normalized_Counts, 
             color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  ylim(0, 250) +
  labs(x = "Aldh1a2", y = NULL) +
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
ggsave("../results/figures/final_plots/Aldh1a2.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics
count_table %>%
  dplyr::filter(gene_name == "Aldh1a2") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.001944



### SERPINE2 ###

# Visualize counts.
count_table %>%
  dplyr::filter(gene_name == "Serpine2") %>%
  ggplot(aes(x = Treatment, 
             y = Normalized_Counts, 
             color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  ylim(80, 200) +
  labs(x = "Serpine2", y = NULL) +
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
ggsave("../results/figures/final_plots/Serpine2.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics
count_table %>%
  dplyr::filter(gene_name == "Serpine2") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.0016



### SHH ###

# Visualize counts.
count_table %>%
  dplyr::filter(gene_name == "Shh") %>%
  ggplot(aes(x = Treatment, 
             y = Normalized_Counts, 
             color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Shh", y = NULL) +
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
ggsave("../results/figures/final_plots/Shh.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics
count_table %>%
  dplyr::filter(gene_name == "Shh") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.0079



### SYNAPSE ASSEMBLY ###

### LRRTM1 ###

# Visualize counts.
count_table %>%
  dplyr::filter(gene_name == "Lrrtm1") %>%
  ggplot(aes(x = Treatment, 
             y = Normalized_Counts, 
             color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Lrrtm1", y = "Normalized Counts") +
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
ggsave("../results/figures/final_plots/Lrrtm1.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics.
count_table %>%
  dplyr::filter(gene_name == "Lrrtm1") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.016



### FLRT2 ###

# Visualize counts.
count_table %>%
  dplyr::filter(gene_name == "Flrt2") %>%
  ggplot(aes(x = Treatment, 
             y = Normalized_Counts, 
             color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Flrt2", y = NULL) +
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
ggsave("../results/figures/final_plots/Flrt2.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics
count_table %>%
  dplyr::filter(gene_name == "Flrt2") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.0055



### ADGRB1 ###

# Visualize counts.
count_table %>%
  dplyr::filter(gene_name == "Adgrb1") %>%
  ggplot(aes(x = Treatment, 
             y = Normalized_Counts, 
             color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Adgrb1", y = NULL) +
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
ggsave("../results/figures/final_plots/Adgrb1.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics
count_table %>%
  dplyr::filter(gene_name == "Adgrb1") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.022



### SYNDIG1 ###

# Visualize counts.
count_table %>%
  dplyr::filter(gene_name == "Syndig1") %>%
  ggplot(aes(x = Treatment, 
             y = Normalized_Counts, 
             color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Syndig1", y = NULL) +
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
ggsave("../results/figures/final_plots/Syndig1.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics
count_table %>%
  dplyr::filter(gene_name == "Syndig1") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.0197



### NEUROPEPTIDE SIGNALING ###

### CPE ###

# Visualize counts.
count_table %>%
  dplyr::filter(gene_name == "Cpe") %>%
  ggplot(aes(x = Treatment, 
             y = Normalized_Counts, 
             color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Cpe", y = "Normalized Counts") +
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
ggsave("../results/figures/final_plots/Cpe.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics.
count_table %>%
  dplyr::filter(gene_name == "Cpe") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.0103



### GRP ###

# Visualize counts.
count_table %>%
  dplyr::filter(gene_name == "Grp") %>%
  ggplot(aes(x = Treatment, 
             y = Normalized_Counts, 
             color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Grp", y = NULL) +
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
ggsave("../results/figures/final_plots/Grp.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics
count_table %>%
  dplyr::filter(gene_name == "Grp") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.022



### NXPH3 ###

# Visualize counts.
count_table %>%
  dplyr::filter(gene_name == "Nxph3") %>%
  ggplot(aes(x = Treatment, 
             y = Normalized_Counts, 
             color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Nxph3", y = NULL) +
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
ggsave("../results/figures/final_plots/Nxph3.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics
count_table %>%
  dplyr::filter(gene_name == "Nxph3") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.01698



### NPY4R ###

# Visualize counts.
count_table %>%
  dplyr::filter(gene_name == "Npy4r") %>%
  ggplot(aes(x = Treatment, 
             y = Normalized_Counts, 
             color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Npy4r", y = NULL) +
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
ggsave("../results/figures/final_plots/Npy4r.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics
count_table %>%
  dplyr::filter(gene_name == "Npy4r") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.0169



