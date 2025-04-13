#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(indicspecies)
library(patchwork) 
library(magick)
library(ggplot2)
library(grid)
library(stringr)

#### Load data ####
load("longmeta_phyloseq_unrare.RData")

#### DESeq ####
Mono_Infected <- subset_samples(mpt_final, hiv_status_clean == "HIV+"& hcv == "NO")
mpt_plus1 <- transform_sample_counts(Mono_Infected, function(x) x+1)
mpt_deseq <- phyloseq_to_deseq2(mpt_plus1, ~`INSTI_drug_current`)
DESEQ_mpt <- DESeq(mpt_deseq)
res <- results(DESEQ_mpt, tidy=TRUE, 
               #this will ensure that No is your reference group
               contrast = c("INSTI_drug_current","YES","NO"))
#View(res)

# Look at results 

## Volcano plot: effect size VS significance
ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
vol_plot <- res %>%
  mutate(significance = case_when(
      padj < 0.01 & log2FoldChange > 2 ~ "Increased",
      padj < 0.01 & log2FoldChange < -2 ~ "Decreased",
      TRUE ~ "Not significant"),
    significance = factor(significance, levels = c("Increased", "Decreased", "Not significant"))) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Increased" = "darkred", "Decreased" = "navyblue", "Not significant" = "gray70")) +
  labs(title = "Differential Abundance Analysis by \n INSTI Status on HIV+ Individuals",
    x = "Log2 Fold Change",
    y = "Adjusted P-value",
    color = "Significance") +
  theme_minimal() +
  theme(legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray80"))

vol_plot


# Bar plot with significant ASVs
# To get table of results
sigASVs <- res %>% 
  filter(padj < 0.01 & abs(log2FoldChange) > 1.7) %>%
  dplyr::rename(ASV = row)

# Get only ASV names 
sigASVs_vec <- sigASVs %>% pull(ASV)

# Prune phyloseq file
mpt_DESeq <- prune_taxa(sigASVs_vec, Mono_Infected)

# Now join with taxonomy info
sigASVs_with_tax <- tax_table(mpt_DESeq) %>% 
  as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs, by="ASV") %>%
  arrange(log2FoldChange) %>%
  filter(!is.na(Genus)) %>%  
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = gsub("g__", "", Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

# Create bar plot
sigASVs_bar <- ggplot(sigASVs_with_tax, aes(x = Genus, y = log2FoldChange, fill = log2FoldChange > 0)) +
  geom_bar(stat = "identity", width = 0.8) +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  scale_fill_manual(values = c("FALSE" = "orange", "TRUE" = "darkgreen"), 
                    labels = c("FALSE" = "Decreased", "TRUE" = "Increased"),
                    name = "Abundance") +
  labs(title = NULL,
    x = NULL,
    y = "Log2 Fold Change") +
  theme_minimal() +
  theme(legend.position = "none",text = element_text(size = 14)) +
  coord_flip()  

sigASVs_bar

library(cowplot)
library(ggplot2)

# First, try loading the image again with a direct path
img_path <- "../03_indicator_species/Indicator_species2.png"


# Create a ggdraw object with the image
img_plot <- ggdraw() + draw_image(img_path)

# Create the combined figure
combined_plot <- img_plot + sigASVs_bar + 
  plot_annotation(tag_levels = 'A',title = "") & 
  theme(plot.tag = element_text(size = 18))

combined_plot

#### Saving files #####
ggsave("Volcano_plot.png",vol_plot)
ggsave("sigASVs_bar.png", sigASVs_bar)
ggsave("FINAL_combined_DESeq.png", 
       combined_plot,
       width = 12,     
       height = 6,    
       dpi = 300,     
       units = "in",  
       limitsize = FALSE,  
       bg = "white")
  

