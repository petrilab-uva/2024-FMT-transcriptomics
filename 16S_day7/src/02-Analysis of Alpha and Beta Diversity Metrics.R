# ANALYSIS OF ALPHA AND BETA DIVERSITY METRICS ---------------------------------
# NAME: G. Brett Moreau
# DATE: December 8, 2023

# PACKAGES ---------------------------------------------------------------------
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install()
packageVersion("BiocManager") # I'm using version 1.30.22

#BiocManager::install("phyloseq")
library(phyloseq)
packageVersion("phyloseq") # I'm using version 1.46.0

#install.packages("tidyverse")
library(tidyverse)
packageVersion("tidyverse") # I'm using version 2.0.0

#BiocManager::install("microbiome")
library(microbiome) 
packageVersion("microbiome") # I'm using version 1.24.0

#install.packages("vegan")
library(vegan)
packageVersion("vegan") # I'm using version 2.6.4


# INTRODUCTION -----------------------------------------------------------------

# 16S sequencing data from the Day 7 FMT experiment has been run through the 
# DADA2 pipeline. I will now evaluate alpha and beta diversity between groups.




# ORGANIZE THE DATA SET --------------------------------------------------------


# Load environment from DADA2 analysis
load(file = "../results/DADA2_Pipeline_R_Workspace_Day7.RData")
rm(list = ls()[!ls() %in% c("ps.analysis")])

# Organize phyloseq object
ps.analysis@sam_data$Group <- factor(ps.analysis@sam_data$Group,
                                     levels = c("Donor", 
                                                "FMT",
                                                "PBS"))




# ALPHA DIVERSITY: RICHNESS ----------------------------------------------------

# I'll start by looking at alpha diversity. I'll evaluate richness using 
# Observed OTUs and evenness using Pielou's Evenness Index.


### RICHNESS ###
plot_richness(ps.analysis, x = "Group", measures = c("Observed"), 
              color = "Group") + 
  geom_boxplot() +
  geom_point(size = 3) +
  theme_bw() +
  ylim(0, 150) +
  labs(x = NULL, y = "Observed ASVs", title = "Observed: Day 7") +
  scale_color_manual(values = c("#8080ff", "#ff8080", "#808080")) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),  
        legend.position = "none", 
        #aspect.ratio = 1, 
        plot.title = element_text(hjust = 0.5, size = 20), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18))

ggsave("../results/figures/richness-observed_d7.png", width = 4, height = 5, dpi = 300)


# Calculate statistics
richness <- estimate_richness(ps.analysis)
pairwise.wilcox.test(richness$Observed, p.adjust.method = "bonferroni", sample_data(ps.analysis)$Group)

# The control (PBS) group was significantly lower than both Donor and FMT groups
# after Bonferroni correction for multiple comparisons. There was no significant
# difference between Donor and FMT recipient groups.




# ALPHA DIVERSITY: EVENNESS ----------------------------------------------------
# Evenness measures are more robust to noise from rare sequence variants and can be more easily 
# interpreted using this pipeline. I'll use Pielou's Evenness as my evenness measure.

# First I'll organize the data.
sample.data <- ps.analysis@sam_data
sample.data <- data.frame(sample.data)

evenness <- evenness(ps.analysis, index = c("pielou"))
evenness <- cbind("Sample.ID" = rownames(evenness), evenness)

sample.data <- full_join(sample.data, evenness, by = "Sample.ID")

### PIELOU ###
ggplot(sample.data, aes(x = Group, y = pielou, color = Group)) +
  geom_boxplot() +
  geom_point(size = 3) +
  theme_bw() +
  labs(x = NULL, y = "Pielou Evenness Index", title = "Evenness: Day 7") +
  scale_color_manual(values = c("#8080ff", "#ff8080", "#808080")) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),  
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 20), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20), 
        legend.text = element_text(size = 18))

#ggsave("../results/figures/evenness-pielou_Timepoint.png", 
#       width = 5, height = 4, dpi = 300)


pairwise.wilcox.test(evenness$pielou, p.adjust.method = "bonferroni", sample_data(ps.analysis)$Group)

# Evenness was significantly lower in the PBS control group compared to Donor 
# mice. It trended lower compared to FMT recipient mice (p=0.08), but there was
# a lot of variability between FMT recipient samples. There was no significant 
# difference between Donor and FMT recipient groups.




# BETA DIVERSITY ---------------------------------------------------------------
# I'll compare beta diversity between samples using Bray-Curtis Dissimilarity.

# Set seed for reproducibility.
#sample(1:10000, 1) # It selected 3365
set.seed(3365) 

# Ordinate using NMDS
ord.nmds.bray <- ordinate(ps.analysis, method = "NMDS", distance = "bray")

# NMDS Plot
ordplot.bray <- plot_ordination(ps.analysis, ord.nmds.bray, color = "Group") +
  geom_point(size = 3) +
  #geom_text(aes(label = Sample.ID, hjust = 2)) +
  labs(title = "Bray-Curtis: Day 7", color = "Treatment") +
  scale_color_manual(values = c("#8080ff", "#ff8080", "#808080" )) +
  theme_bw() +
  theme(aspect.ratio = 1, 
        legend.position = "none", 
        plot.title = element_text(hjust = 0.5, size = 20), 
        axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18)) 

ordplot.bray

#ggsave("../results/figures/beta-diversity-bc.png", 
#       width = 5, height = 4, dpi = 300)



# Add ordination ellipses
ordplot.bray +
  stat_ellipse(type = "norm") +
  theme_bw() +
  theme(axis.text = element_text(size = 10), 
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5))  +
  coord_fixed()


# Run PERMANOVA
adonis2(distance(ps.analysis, method = "bray") ~sample_data(ps.analysis)$Group)


# There was a significant difference; will now look at pairwise comparisons.

permanova.combinations <- combn(x = levels(ps.analysis@sam_data$Group), m = 2)
permanova.p.bray <- vector(length = 3)
for(i in 1:ncol(permanova.combinations)){
  ps.subs.bray <- subset_samples(ps.analysis, Group %in% permanova.combinations[,i])
  metadata.subs.bray <- as(sample_data(ps.subs.bray), "data.frame")
  pairwise.bray <- adonis(distance(ps.subs.bray, method = "bray") ~ Group, data = metadata.subs.bray)
  permanova.p.bray[i] <- pairwise.bray$aov.tab[1,ncol(pairwise.bray$aov.tab)]
}
combinations <- c("Donor + PBS", "FMT + PBS", "Donor + FMT")

permanova.FDR.bray <- p.adjust(permanova.p.bray, method = "bonferroni")
permanova.bray.table <- cbind(combinations, permanova.p.bray)
permanova.bray.table <- cbind(permanova.bray.table, permanova.FDR.bray)

View(permanova.bray.table) # All significant from each other.


