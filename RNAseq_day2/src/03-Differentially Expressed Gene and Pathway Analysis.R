# DIFFERENTIALLY EXPRESSED GENE AND PATHWAY ANALYSIS ---------------------------
# NAME: G. Brett Moreau
# DATE: July 25, 2023

# PACKAGES ---------------------------------------------------------------------
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
load("../results/DESeq2 results.RData")

# Import gene name information for appending results files.
gene_names <- read_gff("../data/Mus_musculus.GRCm39.109.gtf.gz") 
gene_names <- plyranges::select(gene_names, gene_id, gene_name)
gene_names <- as.data.frame(gene_names)
gene_names <- gene_names %>%
  dplyr::select(gene_name, gene_id) %>%
  dplyr::distinct()

# Add gene names to results files
res_FMTvsVEHICLE_table <- as.data.frame(res_FMTvsVEHICLE)
res_FMTvsVEHICLE_table <- res_FMTvsVEHICLE_table %>%
  tibble::rownames_to_column(., var = "gene_id") %>%
  dplyr::left_join(., gene_names, by = "gene_id") %>%
  dplyr::relocate(gene_id, gene_name, .before = baseMean)

# Print DEG table
write.csv(res_FMTvsVEHICLE_table, 
          file = "../results/tables/FMT_vs_control_DEG.csv",
          row.names = FALSE)



# Count Number of DEGs
res_FMTvsVEHICLE_table %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::select(gene_name) # 720 DEGs

# Count Number of DEGs per group
res_FMTvsVEHICLE_table %>%
  dplyr::filter(padj < 0.05) %>% 
  dplyr::filter(log2FoldChange > 0) %>%
  dplyr::select(gene_name) # 351 DEGs upregulated in FMT


res_FMTvsVEHICLE_table %>%
  dplyr::filter(padj < 0.05) %>% 
  dplyr::filter(log2FoldChange < 0) %>%
  dplyr::select(gene_name) # 369 DEGs upregulated in control




# VOLCANO PLOT OF DIFFERENTIALLY EXPRESSED GENES -------------------------------

# Generate column for colors
res_FMTvsVEHICLE_volcano <- res_FMTvsVEHICLE_table %>%
  drop_na(padj) %>%
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

# Save volcano plot
ggsave("../results/figures/volcano_FMT_vs_control_day2.png", 
       dpi = 300, height = 9)




# HEATMAPS OF TOP DIFFERENTIALLY EXPRESSED GENES -------------------------------

# Select 50 most differentially expressed genes.
heatmap_FMTvsVEHICLE <- 
  assay(rld)[ head(order(res_FMTvsVEHICLE$padj), 50), ]

# Organize data set to remove unnecessary samples and add gene names.
heatmap_FMTvsVEHICLE <- as.data.frame(heatmap_FMTvsVEHICLE)
heatmap_FMTvsVEHICLE$gene_id <- rownames(heatmap_FMTvsVEHICLE)
heatmap_FMTvsVEHICLE <- heatmap_FMTvsVEHICLE %>%
  dplyr::select(gene_id, "G1_07", "G1_08", "G1_09", "G1_10", "G2_17", "G2_18", 
                "G2_19", "G2_20") %>%
  dplyr::left_join(., gene_names, by = "gene_id")

# Convert gene names to row names
rownames(heatmap_FMTvsVEHICLE) <- heatmap_FMTvsVEHICLE$gene_name
heatmap_FMTvsVEHICLE <- dplyr::select(heatmap_FMTvsVEHICLE, -gene_id, -gene_name)

# Convert count data to change in expression across the average of all samples.
heatmap_FMTvsVEHICLE <- heatmap_FMTvsVEHICLE - rowMeans(heatmap_FMTvsVEHICLE)
heatmap_FMTvsVEHICLE <- as.matrix(heatmap_FMTvsVEHICLE)

# Format metadata to add to plot
key <- metadata %>% 
  dplyr::select(Treatment) %>%
  plyr::rename(replace = c("Treatment" = "Group"))

key$Treatment[key$Group == "FMT_Recipient"] <- "ABX + FMT"
key$Treatment[key$Group == "Vehicle_Control"] <- "ABX + PBS"

key$Treatment <- factor(key$Treatment, levels = c("ABX + FMT", "ABX + PBS"))
key$Group <- NULL

ann_colors = list(Treatment = c(`ABX + FMT` ="#4097AA", `ABX + PBS` = "#FA9483"))

# Plot heatmap (remove comments to save plot)
#png("../results/figures/heatmap_FMT_vs_control_day2.png", 
#    width = 3000, height = 5000, res = 600, pointsize = 40)
pheatmap(heatmap_FMTvsVEHICLE, 
         annotation_col = key,
         annotation_colors = ann_colors)
#dev.off()




# ORGANIZATION FOR GENE SET ENRICHMENT ANALYSIS --------------------------------

# GSEA will be performed with the fgsea package. I will look for enrichment 
# using the Hallmark gene set from the MSigDB Collection. In order to perform
# GSEA, I need to have a ranked list of genes and the gene set to compare it to.

# I'll start by importing the Hallmark gene set, then organize individual ranked
# gene lists.

# Load Hallmark pathways for mouse
h_gene_set <- msigdbr(species = "Mus musculus", category = "H")

# Organize Hallmark pathway set for GSEA
h_gene_set_list <- split(x = h_gene_set$gene_symbol, 
                         f = h_gene_set$gs_name)

View(h_gene_set_list) # List with genes for each gene set name.


# Organize ranked gene list for each comparison

# Check for missing Wald statistic values
table(is.na(res_FMTvsVEHICLE_table$stat)) # no missing values

# Check for duplicate gene names
table(duplicated(res_FMTvsVEHICLE_table$stat)) # 26 duplicates

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
                            nPermSimple = 10000,
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

# Print GSEA output
write.csv(fgsea_FMTvsVEHICLE_table, 
          file = "../results/tables/GSEA_FMT_vs_control.csv", 
          row.names = FALSE)


# Plot all significantly enriched pathways
fgsea_plot <- fgsea_FMTvsVEHICLE_table %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::mutate(Enrichment = ifelse(NES > 0, "FMT Recipient", "Vehicle Control")) %>%
  dplyr::mutate(Enrichment = factor(Enrichment, 
                                    levels = c("Vehicle Control", "FMT Recipient")))

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
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_MYC_TARGETS_V2"] <- "Myc Targets v2"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_MYOGENESIS"] <- "Myogenesis"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_OXIDATIVE_PHOSPHORYLATION"] <- "Oxidative Phosphorylation"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_P53_PATHWAY"] <- "p53 Pathway"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_PEROXISOME"] <- "Peroxisome"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_TGF_BETA_SIGNALING"] <- "TGFb Signaling"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_TNFA_SIGNALING_VIA_NFKB"] <- "TGFa Signaling via NFkB"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"] <- "Unfolded Protein Response"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_UV_RESPONSE_DN"] <- "UV Response: Down"
fgsea_plot$pathway[fgsea_plot$pathway == "HALLMARK_XENOBIOTIC_METABOLISM"] <- "Xenobiotic Metabolism"

# Filter to include only top 5 pathways in each group by NES
fgsea_FMTvsVEHICLE_Top5up <- fgsea_plot %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(desc(NES)) %>%
  dplyr::slice_head(n = 5)

fgsea_FMTvsVEHICLE_Top5down <- fgsea_plot %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(NES) %>%
  dplyr::slice_head(n = 5)

fgsea_FMTvsVEHICLE_Top10 <- rbind(fgsea_FMTvsVEHICLE_Top5up, fgsea_FMTvsVEHICLE_Top5down)


# Plot enrichment ranking
ggplot(fgsea_FMTvsVEHICLE_Top10, aes(x = NES, y = reorder(pathway, NES), 
                                     size = enriched_genes, color = Enrichment)) +
  geom_point() +
  scale_size(range = c(5, 15)) +
  labs(x = "Normalized Enrichment Score", y = NULL, color = "Treatment", size = "# Genes") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Save plot.
ggsave("../results/figures/GSEA_FMTvsVEHICLE_Top10_Enrichment_Day2.png", 
       dpi = 300, height = 6, width = 10)


# Plot enrichment as bar plot only
ggplot(fgsea_plot, aes (x = reorder(pathway, NES), y = NES, 
                                      fill = Enrichment)) +
  geom_bar(stat = "identity") +
  labs(x = NULL, y = "NES", color = "Treatment") +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 10.5, linetype = "dashed") +
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

# Save bar plot.
ggsave("../results/figures/GSEA_FMTvsVEHICLE_Top10_Enrichment_Day2_horizontal.png", 
       dpi = 300, height = 8, width = 11)




# ENRICHMENT OF IL-4/IL-13 GENE SET BY GSEA ------------------------------------

# I'd like to see if Type 2 immunity genes are broadly upregulated at Day 2 
# along with inflammatory pathways. I'll use the "Interleukin-4 and 
# Interleukin-13 Signaling" gene set from Reactome for this enrichment.

# Load Reactome pathways for humans.
reactome_all <- msigdbr(species = "Mus musculus", 
                        category = "C2", 
                        subcategory = "CP:REACTOME")

# Filter for IL-4/IL-13 Reactome pathway.
reactome_IL4 <- reactome_all %>%
  dplyr::filter(gs_exact_source == "R-HSA-6785807") %>%
  dplyr::distinct(ensembl_gene, .keep_all = TRUE) # 113 unique genes

# Check for overlap with RNAseq gene list.
overlap_FMTvsVEHICLE <- res_FMTvsVEHICLE_table %>%
  dplyr::select(gene_id, gene_name, stat) %>%
  na.omit() %>%
  plyr::rename(replace = c("gene_id" = "ensembl_gene", 
                           "gene_name" = "gene_symbol")) %>%
  dplyr::inner_join(., reactome_IL4, by = "gene_symbol") # 110/113 overlaps.


# Generate ranked gene list for GSEA.
ranked_FMTvsVEHICLE <- res_FMTvsVEHICLE_table %>%
  dplyr::select(gene_name, stat) %>%
  plyr::rename(replace = c("gene_name" = "gene_symbol")) %>%
  na.omit() %>%
  group_by(gene_symbol) %>%
  summarize(stat = mean(stat)) %>%
  deframe()

# Organize Reactome pathway set for GSEA
reactome_IL4_list <- split(x = reactome_IL4$gene_symbol, 
                           f = reactome_IL4$gs_name)


# Run GSEA
fgsea_IL4 <- fgsea(pathways = reactome_IL4_list, 
                   stats = ranked_FMTvsVEHICLE, 
                   eps = 0)

# GSEA says the the IL-4/IL-13 reactome pathway is significantly enriched.

plotEnrichment(reactome_IL4_list[["REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING"]],
               ranked_FMTvsVEHICLE) + 
  labs(x = "Gene Rank", y = "Enrichment Score", title = "IL4 / IL-13 Signaling") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 25, hjust = 0.5),
        axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15), 
        axis.title.x = element_text(size = 25), 
        axis.title.y = element_text(size = 25))

# Save GSEA plot for Il-4/IL-13 pathway.
ggsave("../results/figures/IL4 enrichment plot.png", 
       dpi = 300, width = 6)




# MAKE BOXPLOTS FOR SELECT LEADING EDGE GENES ----------------------------------

# Organize normalized count data
count_table <- counts(dds_object, normalized = TRUE)

count_table <- count_table %>%
  as.data.frame() %>%
  dplyr::mutate(gene_id = rownames(.)) %>%
  tidyr::pivot_longer(cols = !gene_id, 
                      names_to = "Name", 
                      values_to = "Normalized_Counts") %>%
  dplyr::left_join(., gene_names, by = "gene_id") %>%
  dplyr::left_join(., metadata, by = "Name") %>%
  dplyr::filter(Treatment != "FMT_Donor")

count_table$Treatment <- factor(count_table$Treatment, 
                                levels = c("Vehicle_Control", "FMT_Recipient"))



### IMMUNE ACTIVATION GENES ###

### BCL3 ###

# Visualize distribution.
count_table %>%
  dplyr::filter(gene_name == "Bcl3") %>%
  ggplot(aes(x = Treatment, y = Normalized_Counts, color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Bcl3", y = "Normalized Counts") +
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
ggsave("../results/figures/boxplots/Bcl3.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics.
count_table %>%
  dplyr::filter(gene_name == "Bcl3") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.0012



### MYD88 ###

# Visualize distribution.
count_table %>%
  dplyr::filter(gene_name == "Myd88") %>%
  ggplot(aes(x = Treatment, y = Normalized_Counts, color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Myd88", y = NULL) +
  scale_color_manual(values = c("#808080", "#ff8080")) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# Save boxplot.
ggsave("../results/figures/boxplots/MyD88.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics.
count_table %>%
  dplyr::filter(gene_name == "Myd88") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.0028



### STAT3 ###

# Visualize distribution.
count_table %>%
  dplyr::filter(gene_name == "Stat3") %>%
  ggplot(aes(x = Treatment, y = Normalized_Counts, color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Stat3", y = NULL) +
  scale_color_manual(values = c("#808080", "#ff8080")) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

# Save boxplot.
ggsave("../results/figures/boxplots/Stat3.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics.
count_table %>%
  dplyr::filter(gene_name == "Stat3") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.0091



### TISSUE REMODELING ###

### MMP3 ###

# Visualize distribution.
count_table %>%
  dplyr::filter(gene_name == "Mmp3") %>%
  ggplot(aes(x = Treatment, y = Normalized_Counts, color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Mmp3", y = "Normalized Counts") +
  scale_color_manual(values = c("#808080", "#ff8080")) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

# Save boxplot.
ggsave("../results/figures/boxplots/Mmp3.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics.
count_table %>%
  dplyr::filter(gene_name == "Mmp3") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.027




### TIMP1 ###

# Visualize distribution.
count_table %>%
  dplyr::filter(gene_name == "Timp1") %>%
  ggplot(aes(x = Treatment, y = Normalized_Counts, color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Timp1", y = NULL) +
  scale_color_manual(values = c("#808080", "#ff8080")) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) # p = 0.0269

# Save boxplot.
ggsave("../results/figures/boxplots/Timp1.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics.
count_table %>%
  dplyr::filter(gene_name == "Timp1") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.0269



### SERPINE1 ###

# Visualize distribution.
count_table %>%
  dplyr::filter(gene_name == "Serpine1") %>%
  ggplot(aes(x = Treatment, y = Normalized_Counts, color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Serpine1", y = NULL) +
  scale_color_manual(values = c("#808080", "#ff8080")) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

# Save boxplot.
ggsave("../results/figures/boxplots/Serpine1.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics.
count_table %>%
  dplyr::filter(gene_name == "Serpine1") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.10



### TYPE 2 IMMUNITY ###

### IL13RA1 ###

# Visualize distribution.
count_table %>%
  dplyr::filter(gene_name == "Il13ra1") %>%
  ggplot(aes(x = Treatment, y = Normalized_Counts, color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Il13ra1", y = "Normalized Counts") +
  scale_color_manual(values = c("#808080", "#ff8080")) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# Save boxplot.
ggsave("../results/figures/boxplots/Il13ra1.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics.
count_table %>%
  dplyr::filter(gene_name == "Il13ra1") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.014



### IL4RA ###

# Visualize distribution.
count_table %>%
  dplyr::filter(gene_name == "Il4ra") %>%
  ggplot(aes(x = Treatment, y = Normalized_Counts, color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Il4ra", y = NULL) +
  scale_color_manual(values = c("#808080", "#ff8080")) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

# Save boxplot.
ggsave("../results/figures/boxplots/Il4ra.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics.
count_table %>%
  dplyr::filter(gene_name == "Il4ra") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.019



### IL33 ###

# Visualize distribution.
count_table %>%
  dplyr::filter(gene_name == "Il33") %>%
  ggplot(aes(x = Treatment, y = Normalized_Counts, color = Treatment)) +
  geom_boxplot() +
  geom_point(size = 3) +
  labs(x = "Il33", y = NULL) +
  scale_color_manual(values = c("#808080", "#ff8080")) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# Save boxplot.
ggsave("../results/figures/boxplots/Il33.png", 
       dpi = 300, height = 3, width = 3)

# Calculate statistics.
count_table %>%
  dplyr::filter(gene_name == "Il33") %>%
  t.test(Normalized_Counts ~ Treatment, data = .) # p = 0.0026


