#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(DESeq2)
library(indicspecies)


#### Load data ####
load("mpt_final.RData")

# Uncomment the following line to load Monoinfected mpt final with updated file path
# load("mpt_final_monoinfected.RData")

#### DESeq ####
Mono_Infected <- subset_samples(mpt_final, hiv_status_clean == "HIV+"& hcv == "NO")
mpt_plus1 <- transform_sample_counts(Mono_Infected, function(x) x+1)
mpt_deseq <- phyloseq_to_deseq2(mpt_plus1, ~`INSTI_drug_current`)
DESEQ_mpt <- DESeq(mpt_deseq)
res <- results(DESEQ_mpt, tidy=TRUE, 
               #this will ensure that No is your reference group
               contrast = c("INSTI_drug_current","YES","NO"))
View(res)

# Look at results 

## Volcano plot: effect size VS significance
ggplot(res) +
  geom_point(aes(x=log2FoldChange, y=-log10(padj)))

## Make variable to color by whether it is significant + large change
vol_plot <- res %>%
  mutate(significant = padj<0.01 & abs(log2FoldChange)>2) %>%
  ggplot() +
  geom_point(aes(x=log2FoldChange, y=-log10(padj), col=significant))

vol_plot

# Bar plot with significant ASVs
# To get table of results
sigASVs <- res %>% 
  filter(padj<0.01 & abs(log2FoldChange)>1.7) %>%
  dplyr::rename(ASV=row)
View(sigASVs)
# Get only asv names
sigASVs_vec <- sigASVs %>%
  pull(ASV)

# Prune phyloseq file
mpt_DESeq <- prune_taxa(sigASVs_vec,mpt_final)
sigASVs <- tax_table(mpt_DESeq) %>% as.data.frame() %>%
  rownames_to_column(var="ASV") %>%
  right_join(sigASVs) %>%
  arrange(log2FoldChange) %>%
  mutate(Genus = make.unique(Genus)) %>%
  mutate(Genus = factor(Genus, levels=unique(Genus)))

ggplot(sigASVs) +
  geom_bar(aes(x=Genus, y=log2FoldChange), stat="identity")+
  geom_errorbar(aes(x=Genus, ymin=log2FoldChange-lfcSE, ymax=log2FoldChange+lfcSE)) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

#### Saving files #####
ggsave(filename="vol_plot.png",vol_plot)
