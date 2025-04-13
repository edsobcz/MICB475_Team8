### Setup ###
# load packages
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(sf)
library(metagMisc)
library(patchwork)
library(RColorBrewer)
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
  labs(y = "Count", title = NULL, fill = NULL)+
  scale_x_discrete(labels = c("Bacilli", "Bacteroidia", "Clostridia")) +
  scale_fill_brewer(labels = c("Eubacteriaceae",
                               "Alloprevotella",
                               "Anaerostipes",
                               "Fusicatenibacter",
                               "Holdemanella"),
                    palette = "Set3") +
  theme(
    legend.position = "top", 
    text = element_text(size = 10),
    legend.key.width = unit(0.5, "cm"),
    legend.text = element_text(size = 10),
    legend.spacing.x = unit(0.2, 'cm'),
    legend.box.spacing = unit(0.1, "cm")) +
  guides(fill = guide_legend(ncol = 2))

core_microbiome_bar


taxabar <- plot_bar(unique_insti_genus_ps, fill = "Genus")
taxabar

# Venn diagram
# Get unique genera (not ASVs) for each group
pos_genera <- unique(tax_table(subset_taxa(INSTI_pos, taxa_names(INSTI_pos) %in% pos_ASVs))[,"Genus"])
neg_genera <- unique(tax_table(subset_taxa(INSTI_neg, taxa_names(INSTI_neg) %in% neg_ASVs))[,"Genus"])

# Create the Venn diagram at genus level
genus_venn <- ggVennDiagram(x=list(INSTI = pos_genera, NoINSTI = neg_genera),
                            category.names = c("INSTI+", "INSTI-")) +
  scale_fill_gradient(low = "lightgray", high = "magenta4") +
  theme(legend.position = "none", text = element_text(size = 20)) +
  coord_flip()

genus_venn

# Create the combined figure
combined_plot <- genus_venn + core_microbiome_bar + 
  plot_annotation(tag_levels = 'A',title = "") & 
  theme(plot.tag = element_text(size = 18))

combined_plot


# save objects
ggsave("core_microbiome_venn_genus.png", genus_venn)
ggsave("core_microbiome_bar.png", core_microbiome_bar)
ggsave("FINAL_combined_core_microbiome.png", 
       combined_plot,
       width = 12,     
       height = 6,    
       dpi = 300,     
       units = "in",  
       limitsize = FALSE,  
       bg = "white")  

write_csv(unique_insti_genus, "insti_core_microbiome.csv")
