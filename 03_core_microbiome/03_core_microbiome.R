### Setup ###
# load packages
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(sf)

# load data
load("mpt_final.RData")

# Keep only HIV monoinfected
mpt_mono <- subset_samples(mpt_final, hcv == "NO" & hiv_status_clean == "HIV+")

### Core Microbiome ###
# Convert to relative abundance
mpt_RA <- transform_sample_counts(mpt_mono, fun=function(x) x/sum(x))

# Subset INSTI and no INSTI
INSTI_pos <- subset_samples(mpt_RA, INSTI_drug_current == "YES")
INSTI_neg <- subset_samples(mpt_RA, INSTI_drug_current == "NO")

# Check ASVs
pos_ASVs <- core_members(INSTI_pos, detection = 0.005, prevalence = 0.35)
neg_ASVs <- core_members(INSTI_neg, detection = 0.005, prevalence = 0.35)

# detection of 0.01 and prevalence of 0.35 is too high

pos_pruned <- prune_taxa(pos_ASVs, mpt_final) |>
  tax_table()
neg_pruned <- prune_taxa(neg_ASVs, mpt_final) |>
  tax_table()

# Plot relative abundance
prune_taxa(pos_ASVs, mpt_final) |> 
  plot_bar(fill="Family") + 
  facet_wrap(.~`INSTI_drug_current`, scales ="free")

# Venn diagram
insti_venn <- ggVennDiagram(x=list(INSTI = pos_pruned, NoINSTI = neg_pruned))
insti_venn
