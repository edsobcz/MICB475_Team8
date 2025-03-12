### Setup ###
# load packages
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)

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
pos_ASVs <- core_members(INSTI_pos, detection = 0.01, prevalence = 0.7)
neg_ASVs <- core_members(INSTI_neg, detection = 0, prevalence = 0.7)

### this is where errors show up: when detection is changed to above 0, the prune_taxa lines that
### call on that object throw up an error. somebody please figure out what I'm doing wrong and tell
### me I'm being an idiot. please and thank you.

prune_taxa(pos_ASVs, mpt_final) |>
  tax_table()
prune_taxa(neg_ASVs, mpt_final) |>
  tax_table()

# Plot relative abundance
prune_taxa(pos_ASVs, mpt_final) |> 
  plot_bar(fill="Family") + 
  facet_wrap(.~`INSTI_drug_current`, scales ="free")
