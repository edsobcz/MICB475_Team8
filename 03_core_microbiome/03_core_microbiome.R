### Setup ###
# load packages
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(sf)
library(metagMisc)
set.seed(1234)

# load data
load("longmeta_phyloseq_unrare.RData")

# Keep only HIV monoinfected
mono <- subset_samples(mpt_final, hcv == "NO" & hiv_status_clean == "HIV+")

# Glom to Genus and converty to relative abundance
phyloseq_RA <- mono %>% tax_glom("Genus", NArm = FALSE) %>% transform_sample_counts(fun=function(x) x/sum(x))

### Core Microbiome ###
# Subset INSTI and no INSTI
INSTI_pos <- subset_samples(phyloseq_RA, INSTI_drug_current == "YES")
INSTI_neg <- subset_samples(phyloseq_RA, INSTI_drug_current == "NO")

# Check ASVs
pos_ASVs <- core_members(INSTI_pos, detection = 0.005, prevalence = 0.35)
neg_ASVs <- core_members(INSTI_neg, detection = 0.005, prevalence = 0.35)


diffs <- setdiff(pos_ASVs, neg_ASVs)
core_pos_taxa <- prune_taxa(diffs, phyloseq_RA)

# detection of 0.005 and prevalence of 0.35 for prev analysis

pos_pruned <- prune_taxa(pos_ASVs, mono) |>
  tax_table()

neg_pruned <- prune_taxa(neg_ASVs, mono) |>
  tax_table()


# Coerce OTU table with taxa to df
core_pos_taxa <- prune_taxa(pos_ASVs, phyloseq_RA)
pruned_pos_df <- phyloseq_to_df(core_pos_taxa, addtax = T, addtot = T, addmaxrank = F) |> 
  mutate(INSTI_group = "INSTI")

core_neg_taxa <- prune_taxa(neg_ASVs, phyloseq_RA)
pruned_neg_df <- phyloseq_to_df(core_neg_taxa, addtax = T, addtot = T, addmaxrank = F) |> 
  mutate(INSTI_group = "no INSTI")

pruned_df <- rbind(pruned_pos_df, pruned_neg_df)

neg_vec <- pruned_neg_df |> pull(Genus)
pos_vec <- pruned_pos_df |> pull(Genus)

# find the names of the unique insti taxa (by genus)
unique_insti_genus <- subset_taxa(core_pos_taxa, !(Genus %in% neg_vec)) |> 
  phyloseq_to_df(addtax = T, addtot = T, addmaxrank = F)

unique_insti_genus_ps <- subset_taxa(core_pos_taxa, !(Genus %in% neg_vec)) |> 
  transform_sample_counts(fun=function(x) x/sum(x)) |> 
  subset_samples(INSTI_drug_current == "YES")

# Plot core taxa 
core_microbiome_bar <- unique_insti_genus |> 
  ggplot(aes(x = Class, fill = Genus)) +
  geom_bar(position = "stack") +
 # facet_wrap(.~INSTI_group, scales = "free") +
  labs(y = "Count", title = "INSTI Core Microbiome Members")
core_microbiome_bar


taxabar <- plot_bar(unique_insti_genus_ps, fill = "Genus")
taxabar

# Venn diagram
insti_venn <- ggVennDiagram(x=list(INSTI = pos_pruned, NoINSTI = neg_pruned),
                            category.names = c("INSTI",
                                               "no INSTI")) +
  scale_fill_gradient(low = "lightgray", high = "magenta4") +
  theme(legend.position = "none", text = element_text(size = 2)) +
  coord_flip()
insti_venn

# save objects
ggsave("core_microbiome_venn_genus.png", insti_venn)
ggsave("core_microbiome_bar.png", core_microbiome_bar)
write_csv(unique_insti_genus, "insti_core_microbiome.csv")
