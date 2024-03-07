# COMPARISON OF MICROBIAL DIVERSITY FOR DAY 2 FMT SAMPLES ----------------------
# NAME: G. Brett Moreau
# DATE: December 1, 2023

# PACKAGES ---------------------------------------------------------------------

#BiocManager::install("phyloseq")
library(phyloseq)
packageVersion("phyloseq") # Version 1.46.0

#install.packages("tidyverse")
library(tidyverse)
packageVersion("tidyverse") # Version 2.0.0

#install.packages("vegan")
library(vegan)
packageVersion("vegan") # Version 2.6.4


# INTRODUCTION -----------------------------------------------------------------
# 16S sequencing reads from the Day 2 FMT experiment have been processed using
# DADA2. I will now analyze differences in alpha and beta diversity across 
# treatment groups.




# SETTING UP THE ENVIRONMENT ---------------------------------------------------

# Load phyloseq data and keep essential objects for analysis.
load("../results/DADA2_Pipeline_R_Workspace_Day2.RData")
rm(list = ls()[!ls() %in% c("ps.analysis")])

# Remove control samples
ps.analysis <- subset_samples(ps.analysis, experiment_date != "EMPTY")
ps.analysis <- prune_taxa(taxa_sums(ps.analysis) > 0, ps.analysis)

View(ps.analysis@sam_data)
# ALPHA DIVERSITY --------------------------------------------------------------

# I'll compare microbial richness using Observed ASVs and microbial evenness 
# using Pielou's Evenness Index. I hypothesize that antibiotic treatment will
# significantly reduce richness, and FMT treatment will increase richness.

# Rename groups for graphing.
ps.analysis@sam_data$experimental_group <- gsub(fixed("FMT_DONOR"), "Donor", ps.analysis@sam_data$experimental_group)
ps.analysis@sam_data$experimental_group <- gsub(fixed("FMT_Recipient"), "FMT", ps.analysis@sam_data$experimental_group)
ps.analysis@sam_data$experimental_group <- gsub(fixed("Vehicle_recipient"), "PBS", ps.analysis@sam_data$experimental_group)

ps.analysis@sam_data$experimental_group <- factor(ps.analysis@sam_data$experimental_group,
                                                  levels = c("Donor", "FMT", "PBS"))


### RICHNESS ###
plot_richness(ps.analysis, x = "experimental_group", measures = c("Observed"), 
              color = "experimental_group") + 
  geom_boxplot() +
  geom_point(size = 3) +
  theme_bw() +
  ylim(0, 125) +
  labs(x = NULL, y = "Observed ASVs", title = "Observed: Day 2") +
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

#ggsave("../results/figures/richness-observed_day2.png", 
#       width = 4, height = 5, dpi = 300)

# As expected, richness was reduced in the PBS control group and closer to the
# levels seen in donor mice in the FMT-treated group.


# Perform statistics
richness <- estimate_richness(ps.analysis)
pairwise.wilcox.test(richness$Observed, p.adjust.method = "bonferroni", 
                     sample_data(ps.analysis)$experimental_group)

#All groups were significantly different from each other using a Wilcoxon rank
# sum test with Bonferroni correction for multiple comparisons.




# BETA DIVERSITY ---------------------------------------------------------------

# Bray-Curtis dissimilarity will be used to examine differences in beta 
# diversity between groups.


set.seed(405) # Set seed for reproducibility. 

# Organize data
ps.analysis.prop <- transform_sample_counts(ps.analysis, function(ASV) ASV/sum(ASV))
ord.nmds.bray <- ordinate(ps.analysis.prop, method = "NMDS", distance = "bray")

# Plot NMDS Plot
ordplot.bray <- plot_ordination(ps.analysis.prop, ord.nmds.bray, color = "experimental_group") +
  geom_point(size = 3) +
  labs(title = "Bray-Curtis: Day 2", color = "Treatment") +
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

#ggsave("../results/figures/bray-curtis_day2.png", 
#       width = 5, height = 4, dpi = 300)


# Run PERMANOVA
adonis2(distance(ps.analysis.prop, method = "bray") ~sample_data(ps.analysis)$experimental_group)

# Significant difference between groups.

# Perform pairwise comparisons
permanova.combinations <- combn(x = levels(ps.analysis.prop@sam_data$experimental_group), m = 2)
permanova.p.bray <- vector(length = 3)
for(i in 1:ncol(permanova.combinations)){
  ps.subs.bray <- subset_samples(ps.analysis.prop, experimental_group %in% permanova.combinations[,i])
  metadata.subs.bray <- as(sample_data(ps.subs.bray), "data.frame")
  pairwise.bray <- adonis(distance(ps.subs.bray, method = "bray") ~ experimental_group, data = metadata.subs.bray)
  permanova.p.bray[i] <- pairwise.bray$aov.tab[1,ncol(pairwise.bray$aov.tab)]
}
combinations <- c("Donor/FMT", "Donor/PBS", "FMT/PBS")

permanova.FDR.bray <- p.adjust(permanova.p.bray, method = "bonferroni")
permanova.bray.table <- cbind(combinations, permanova.p.bray)
permanova.bray.table <- cbind(permanova.bray.table, permanova.FDR.bray)

View(permanova.bray.table)

# All comparisons are significant after correcting for multiple comparisons, 
# indicating that Bray-Curtis dissimilarity is significantly different for each
# group compared to every other group.
