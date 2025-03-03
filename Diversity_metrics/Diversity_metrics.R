#!/usr/bin/env Rscript
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)

# You guys can just load the RData objects I created and start the analysis
# or you can recreate the objects from scratch, you just need to download
# the files from the server and put them in the same folder as the Rproject
# or update file paths accordingly.

#### Load in RData ####
load("mpt_rare.RData")
load("mpt_final.RData")

#### Create RObject from the beginning #### 

## Load data ##
# Change file paths as necessary
meta <- read_delim("modified_metadata.tsv", delim="\t")
otu <- read_delim(file = "feature-table.txt", delim="\t", skip=1)
tax <- read_delim("taxonomy.tsv", delim="\t")
phylotree <- read.tree("tree.nwk")

## Format Files ##

# OTU table
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

#Metadata 
samp_df <- as.data.frame(meta[,-1])
rownames(samp_df)<- meta$'sample-id'
SAMP <- sample_data(samp_df)

#Taxonomy 
# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature ID`
TAX <- tax_table(tax_mat)

# Create phyloseq object
mpt <- phyloseq(OTU, SAMP, TAX, phylotree)

######### ANALYZE ##########
# Remove non-bacterial sequences, if any
mpt_filt <- subset_taxa(mpt,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
mpt_filt_nolow <- filter_taxa(mpt_filt, function(x) sum(x)>5, prune = TRUE)
# Remove Blanks from the dataset 
mpt_no_blanks <- subset_samples(mpt_filt_nolow, sample_type != "control blank") 
# Remove samples with less than 100 reads
mpt_final <- prune_samples(sample_sums(mpt_no_blanks)>100, mpt_no_blanks) 
mpt_rare <- rarefy_even_depth(mpt_final, rngseed = 124, sample.size = 25228)


#### Alpha diversity exploratory attempts ####

# ----- Alpha diversity HIV+ vs HIV-&HCV- --------

mpt_rare_HIVpos <- subset_samples(mpt_rare, hiv_status_clean == "HIV+")
mpt_rare_healthy <- subset_samples(mpt_rare, hiv_status_clean == "HIV-" & hcv == "NO") 

# Calculate alpha diversity for each phyloseq object
rich_HIVpos <- estimate_richness(mpt_rare_HIVpos, measures = c("Shannon"))
rich_healthy <- estimate_richness(mpt_rare_healthy, measures = c("Shannon"))

# Add group labels
rich_HIVpos$group <- "HIV+"
rich_healthy$group <- "HIV- & HCV-"

# Add sample names as a column
rich_HIVpos$sample <- rownames(rich_HIVpos)
rich_healthy$sample <- rownames(rich_healthy)

# Combine the data
combined_richness <- rbind(rich_HIVpos, rich_healthy)

# Plot
ggplot(combined_richness, aes(x = group, y = Shannon, color = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  theme_bw() +
  labs(title = "Shannon Diversity by HIV Status",
       x = "Group",
       y = "Shannon Diversity Index") +
  theme(legend.position = "none")

# ----- Alpha diversity HIV+ dep+ vs HIV+ dep+ --------

mpt_HIVpos_depPos <- subset_samples(mpt_rare_HIVpos, bdi_ii > 20)
mpt_HIVpos_depNeg <- subset_samples(mpt_rare_HIVpos, bdi_ii <= 20)

# Calculate alpha diversity for each phyloseq object
rich_hp_dp <- estimate_richness(mpt_HIVpos_depPos, measures = c("Shannon"))
rich_hp_dn <- estimate_richness(mpt_HIVpos_depNeg, measures = c("Shannon"))

# Add group labels
rich_hp_dp$group <- "HIV+/dep+"
rich_hp_dn$group <- "HIV+/dep-"

# Combine the data
combined_richness_dep <- rbind(rich_hp_dp, rich_hp_dn)

# Plot
ggplot(combined_richness_dep, aes(x = group, y = Shannon, color = group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  theme_bw() +
  labs(title = "Shannon Diversity HIV and Depression",
       x = "Group",
       y = "Shannon Diversity Index") +
  theme(legend.position = "none")

# ----- Alpha diversity HIV+ dep+ vs HIV+ dep+ --------

mpt_HIVpos_instiPos <- subset_samples(mpt_rare_HIVpos, INSTI_drug_current == "YES")
mpt_HIVpos_instiNeg <- subset_samples(mpt_rare_HIVpos, INSTI_drug_current == "NO")

# Calculate alpha diversity for each phyloseq object
rich_hp_ip <- estimate_richness(mpt_HIVpos_instiPos, measures = c("Shannon"))
rich_hp_in <- estimate_richness(mpt_HIVpos_instiNeg, measures = c("Shannon"))

# Add group labels
rich_hp_ip$group <- "HIV+/INSTI+"
rich_hp_in$group <- "HIV+/INSTI-"

# Combine the data
combined_richness_insti <- rbind(rich_hp_ip, rich_hp_in)

# Plot
ggplot(combined_richness_insti, aes(x = group, y = Shannon, color = group)) +
  geom_boxplot() +
  #geom_jitter(width = 0.2, alpha = 0.6) +
  theme_bw() +
  labs(title = "Shannon Diversity HIV and INSTI",
       x = "Group",
       y = "Shannon Diversity Index") +
  theme(legend.position = "none")


##### Saving #####
save(mpt_final, file="mpt_final.RData")
save(mpt_rare, file="mpt_rare.RData")
