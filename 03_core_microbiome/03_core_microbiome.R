### Setup ###
# load packages
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(sf)
library(metagMisc)

# load data
load("mpt_final.RData")

# Keep only HIV monoinfected
mpt_mono <- subset_samples(mpt_rare, hcv == "NO" & hiv_status_clean == "HIV+")

# Glom to Genus and converty to relative abundance
# mpt_genus <- tax_glom(mpt_mono, "Genus", NArm = FALSE)
mpt_RA <- transform_sample_counts(mpt_mono, fun=function(x) x/sum(x))

### Core Microbiome ###
# Subset INSTI and no INSTI
INSTI_pos <- subset_samples(mpt_RA, INSTI_drug_current == "YES")
INSTI_neg <- subset_samples(mpt_RA, INSTI_drug_current == "NO")

# Check ASVs
pos_ASVs <- core_members(INSTI_pos, detection = 0.005, prevalence = 0.35)
neg_ASVs <- core_members(INSTI_neg, detection = 0.005, prevalence = 0.35)

# detection of 0.005 and prevalence of 0.35 for prev analysis

pos_pruned <- prune_taxa(pos_ASVs, mpt_final) |>
  tax_table()
neg_pruned <- prune_taxa(neg_ASVs, mpt_final) |>
  tax_table()

# Coerce OTU table with taxa to df
core_pos_taxa <- prune_taxa(pos_ASVs, mpt_RA)
pruned_pos_df <- phyloseq_to_df(core_pos_taxa, addtax = T, addtot = T, addmaxrank = F) |> 
  mutate(INSTI_group = "INSTI")

core_neg_taxa <- prune_taxa(pos_ASVs, mpt_RA)
pruned_neg_df <- phyloseq_to_df(core_neg_taxa, addtax = T, addtot = T, addmaxrank = F) |> 
  mutate(INSTI_group = "no INSTI")

pruned_df <- rbind(pruned_pos_df, pruned_neg_df)

# Coerce taxa table to df
# OTU_pos <- as(otu_table(core_pos_taxa), "matrix")
# if(taxa_are_rows(core_pos_taxa)){OTU_pos <- t(OTU_pos)}
# pruned_OTU_pos_df <- as.data.frame(OTU_pos)
  
# Plot core taxa 
core_microbiome_bar <- pruned_df |>
  ggplot(aes(x = Class, fill = Family)) +
  geom_bar() +
  facet_wrap(.~INSTI_group, scales = "free") +
  labs(y = "Count", title = "INSTI Core Microbiome Members")
core_microbiome_bar

ggsave("core_microbiome_bar.png", core_microbiome_bar)

# Venn diagram
insti_venn <- ggVennDiagram(x=list(INSTI = pos_pruned, NoINSTI = neg_pruned)) +
  
insti_venn

ggsave("core_microbiome_venn.png", insti_venn)
