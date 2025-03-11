#!/usr/bin/env Rscript
library(phyloseq)
library(ape) 
library(tidyverse)
library(vegan)
library(stringr)
library(nlme)
library(picante)

# If you ran the script and already saved the RData files, you can just reload them.

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

# Subset samples for subsequent analysis
Infected <- subset_samples(mpt_rare, hiv_status_clean == "HIV+" & hcv == "NO")
Healthy <- subset_samples(mpt_rare, hiv_status_clean == "HIV-" & hcv == "NO")
INSTI_pos <- subset_samples(mpt_rare, hiv_status_clean == "HIV+"& INSTI_drug_current == "YES")
INSTI_neg <- subset_samples(mpt_rare, hiv_status_clean == "HIV+"& INSTI_drug_current == "NO")


#### Pie Chart ####
# ----- HIV+ vs Healthy -----
# Count the number of samples in each group
count_HIVpos <- nsamples(Infected)
count_healthy <- nsamples(Healthy)
count_remaining <- nsamples(mpt_rare) - count_healthy - count_HIVpos

# Create a data frame for plotting
pie_data <- data.frame(
  group = c("Infected", "Healthy", "Remaining"),
  count = c(count_HIVpos, count_healthy, count_remaining)
)

# Calculate percentages
pie_data$percentage <- round(pie_data$count / sum(pie_data$count) * 100, 1)
pie_data$label <- paste0(pie_data$group, "\n", pie_data$percentage, "%")

# Create pie chart with ggplot2
ggplot(pie_data, aes(x = "", y = count, fill = group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("Infected" = "firebrick", "Healthy" = "steelblue", "Remaining" = "yellow")) +
  labs(title = "Percentage of HIV+ vs Healthy Samples") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.3)) +
  theme(legend.position = "none")

# ----- HIV+ INSTI + vs HIV+ INSTI- -----

# Count the number of samples in each group
count_INSTIpos <- nsamples(INSTI_pos)
count_INSTIneg <- nsamples(INSTI_neg)

# Create a data frame for plotting
pie_data <- data.frame(
  group = c("INSTI+", "INSTI-"),
  count = c(count_INSTIpos,count_INSTIneg))

# Calculate percentages
pie_data$percentage <- round(pie_data$count / sum(pie_data$count) * 100, 1)
pie_data$label <- paste0(pie_data$group, "\n", pie_data$percentage, "%")

# Create pie chart with ggplot2
ggplot(pie_data, aes(x = "", y = count, fill = group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("INSTI+" = "darkgreen", "INSTI-" = "darkorange")) +
  labs(title = "Percentage of INSTI+ vs INSTI- \n        in the HIV+ Samples") +
  geom_text(aes(label = label), position = position_stack(vjust = 0.3)) +
  theme(legend.position = "none")

# ---- Chi-squared test -----

# Further divide into depressed and non-depressed groups
INSTI_pos_dep <- subset_samples(INSTI_pos, bdi_ii > 20)
INSTI_pos_nodep <- subset_samples(INSTI_pos, bdi_ii <= 20)

# Calculate richness for each group
rich_pos_dep <- estimate_richness(INSTI_pos_dep, measures = c("Shannon"))
rich_pos_nodep <- estimate_richness(INSTI_pos_nodep, measures = c("Shannon"))

# Categorize richness into "High" and "Low" based on median of all INSTI+ samples
all_richness <- c(rich_pos_dep$Shannon, rich_pos_nodep$Shannon)
median_richness <- median(all_richness, na.rm = TRUE)

# Create category vectors
rich_category_dep <- ifelse(rich_pos_dep$Shannon > median_richness, "High", "Low")
rich_category_nodep <- ifelse(rich_pos_nodep$Shannon > median_richness, "High", "Low")

# Create contingency table
contingency_table <- table(
  Depression = c(rep("Depressed", length(rich_category_dep)), 
                 rep("Not Depressed", length(rich_category_nodep))),
  Richness = c(rich_category_dep, rich_category_nodep)
)

# Perform chi-squared test
chi_result <- chisq.test(contingency_table)
print(chi_result)


#----- Wilcoxon test for INSTI+/INSTI- ------
# Calculate Shannon diversity
rich_pos <- estimate_richness(INSTI_pos, measures = c("Shannon"))
rich_neg <- estimate_richness(INSTI_neg, measures = c("Shannon"))

wilcox.test(rich_pos$Shannon, rich_neg$Shannon)


#### Alpha Diversity ####

# ----- Shannon Diversity HIV+/Healthy  --------

# Calculate alpha diversity for each phyloseq object
rich_HIVpos <- estimate_richness(Infected, measures = c("Shannon"))
rich_healthy <- estimate_richness(Healthy, measures = c("Shannon"))

# Add group labels
rich_HIVpos$group <- "HIV+"
rich_healthy$group <- "Healthy"

# Add sample names as a column
rich_HIVpos$sample <- rownames(rich_HIVpos)
rich_healthy$sample <- rownames(rich_healthy)

# Combine the data
combined_richness <- rbind(rich_HIVpos, rich_healthy)

# Plot
ggplot(combined_richness, aes(x = group, y = Shannon, fill = group)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Shannon Diversity by HIV Status",
       y = "Shannon DiversityD",
       x = "") +
  scale_fill_manual(values = c("HIV+" = "firebrick", "Healthy" = "steelblue"))

#----- Faiths Phylogenetic Diversity HIV+/Healthy ------

# Extract OTU tables and tree
otu_table_infected <- otu_table(Infected)
otu_table_healthy <- otu_table(Healthy)
tree <- phy_tree(mpt_rare)

# If taxa are rows in your OTU table, transpose them
if(taxa_are_rows(otu_table_infected)) {
  otu_table_infected <- t(otu_table_infected)
}
if(taxa_are_rows(otu_table_healthy)) {
  otu_table_healthy <- t(otu_table_healthy)
}

# Calculate Faith's PD
pd_infected <- pd(otu_table_infected, tree)
pd_healthy <- pd(otu_table_healthy, tree)

# Wilcoxon test
wilcox.test(pd_infected$PD, pd_healthy$PD)

# Create a boxplot to visualize
pd_combined <- data.frame(
  Group = c(rep("HIV+", nrow(pd_infected)), rep("Healthy", nrow(pd_healthy))),
  Faith_PD = c(pd_infected$PD, pd_healthy$PD)
)

ggplot(pd_combined, aes(x = Group, y = Faith_PD, fill = Group)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Faith's Phylogenetic Diversity by HIV Status",
       y = "Faith's PD",
       x = "") +
  scale_fill_manual(values = c("HIV+" = "firebrick", "Healthy" = "steelblue"))

#### Beta Diversity ####

# ---- Unweighted Unifrac HIV+/Healthy -----
# Calculate unweighted UniFrac distances
unifrac_dist <- phyloseq::distance(mpt_rare, method = "unifrac", weighted = FALSE)

# Filter to include only HIV+ and healthy (HIV-/HCV-) samples
hiv_healthy_samples <- sample_names(subset_samples(mpt_rare, 
                                                   hiv_status_clean == "HIV+" | 
                                                     (hiv_status_clean == "HIV-" & hcv == "NO")))

# Subset the distance matrix to include only these samples
unifrac_dist_subset_hh <- as.dist(as.matrix(unifrac_dist)[hiv_healthy_samples, hiv_healthy_samples])

# Get the metadata for these samples
meta_subset_hh <- data.frame(sample_data(mpt_rare))[hiv_healthy_samples, ]

# Create a factor for HIV status (making sure "HIV-" shows as "Healthy")
meta_subset_hh$group_factor <- ifelse(
  meta_subset_hh$hiv_status_clean == "HIV+", 
  "HIV+", 
  "Healthy"
)

# PERMANOVA test
permanova_result_hh <- adonis2(unifrac_dist_subset_hh ~ group_factor, data = meta_subset_hh)
print(permanova_result_hh)

# Create a phyloseq subset for visualization
hiv_healthy_phyloseq <- subset_samples(mpt_rare, 
                                       hiv_status_clean == "HIV+" | 
                                         (hiv_status_clean == "HIV-" & hcv == "NO"))

# Add a clear group label for plotting
sample_data(hiv_healthy_phyloseq)$Group <- ifelse(
  sample_data(hiv_healthy_phyloseq)$hiv_status_clean == "HIV+", 
  "HIV+", 
  "Healthy"
)

# Run the ordination
pcoa_hh <- ordinate(hiv_healthy_phyloseq, 
                    method = "PCoA", 
                    distance = "unifrac", 
                    weighted = FALSE)

# Plot the ordination
hh_plot <- plot_ordination(hiv_healthy_phyloseq, 
                           pcoa_hh, 
                           color = "Group") + 
  theme_bw() +
  stat_ellipse(type = "norm") +
  labs(title = "Unweighted UniFrac PCoA: HIV+ vs Healthy") +
  scale_color_manual(values = c("HIV+" = "firebrick", "Healthy" = "turquoise")) +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(hh_plot)


# ---- Unweighted Unifrac INSTI+/INSTI- -----
# Extract sample names for INSTI+ and INSTI- groups
insti_samples <- sample_names(subset_samples(mpt_rare, 
                                             (hiv_status_clean == "HIV+" & INSTI_drug_current == "YES")| 
                                               (hiv_status_clean == "HIV+" & INSTI_drug_current == "NO")))

# Get metadata for these samples
meta_subset_insti <- data.frame(sample_data(mpt_rare))[insti_samples, ]

# Create a factor for INSTI status directly from the metadata
# This ensures the factor matches exactly with the samples in your distance matrix
meta_subset_insti$insti_factor <- ifelse(meta_subset_insti$INSTI_drug_current == "YES", "INSTI+", "INSTI-")

# Subset the distance matrix
unifrac_dist_subset_insti <- as.dist(as.matrix(unifrac_dist)[insti_samples, insti_samples])

# Run PERMANOVA
permanova_result <- adonis2(unifrac_dist_subset_insti ~ insti_factor, data = meta_subset_insti)
print(permanova_result)

# Create a subset of your phyloseq object with only HIV+ samples that have INSTI data
insti_phyloseq <- subset_samples(mpt_rare, 
                                 (hiv_status_clean == "HIV+" & INSTI_drug_current == "YES") | 
                                   (hiv_status_clean == "HIV+" & INSTI_drug_current == "NO"))

# Add a clear INSTI label for plotting
sample_data(insti_phyloseq)$INSTI_status <- ifelse(
  sample_data(insti_phyloseq)$INSTI_drug_current == "YES", 
  "INSTI+", 
  "INSTI-"
)

# Run the ordination using your pre-calculated UniFrac distances
# We need to extract the relevant subset of the distance matrix
insti_samples <- sample_names(insti_phyloseq)
unifrac_insti <- as.dist(as.matrix(unifrac_dist)[insti_samples, insti_samples])

# Create the ordination
pcoa_insti <- ordinate(insti_phyloseq, 
                       method = "PCoA", 
                       distance = unifrac_insti)

# Plot the ordination
insti_plot <- plot_ordination(insti_phyloseq, 
                              pcoa_insti, 
                              color = "INSTI_status") + 
  theme_bw() +
  stat_ellipse(type = "norm") +
  labs(title = "Unweighted UniFrac PCoA by INSTI Status in \n   HIV+ Individuals") +
  scale_color_manual(values = c("INSTI+" = "darkgreen", "INSTI-" = "darkorange"))

# Display the plot
print(insti_plot)


# ---- Unweighted Unifrac INSTI+/Healthy -----
# Extract sample names for INSTI+ and healthy groups
insti_healthy_samples <- sample_names(subset_samples(mpt_rare, 
                                             (hiv_status_clean == "HIV+" & INSTI_drug_current == "YES")| 
                                               (hiv_status_clean == "HIV-" & hcv == "NO")))

# Subset the distance matrix to include only these samples
unifrac_dist_subset_ih <- as.dist(as.matrix(unifrac_dist)[insti_healthy_samples, insti_healthy_samples])

# Get the metadata for these samples
meta_subset_ih <- data.frame(sample_data(mpt_rare))[insti_healthy_samples, ]

# Create a factor for grouping
meta_subset_ih$group_factor <- ifelse(
  meta_subset_ih$hiv_status_clean == "HIV+" & meta_subset_ih$INSTI_drug_current == "YES", 
  "INSTI+", 
  "Healthy"
)

# PERMANOVA test
permanova_result_ih <- adonis2(unifrac_dist_subset_ih ~ group_factor, data = meta_subset_ih)
print(permanova_result_ih)

# Create a subset of your phyloseq object for visualization
insti_healthy_phyloseq <- subset_samples(mpt_rare, 
                                         (hiv_status_clean == "HIV+" & INSTI_drug_current == "YES") | 
                                           (hiv_status_clean == "HIV-" & hcv == "NO"))

# Add a clear group label for plotting
sample_data(insti_healthy_phyloseq)$Group <- ifelse(
  sample_data(insti_healthy_phyloseq)$hiv_status_clean == "HIV+" & 
    sample_data(insti_healthy_phyloseq)$INSTI_drug_current == "YES", 
  "INSTI+", 
  "Healthy"
)

# Run the ordination using your pre-calculated UniFrac distances
# We need to extract the relevant subset of the distance matrix
unifrac_ih <- as.dist(as.matrix(unifrac_dist)[insti_healthy_samples, insti_healthy_samples])

# Create the ordination
pcoa_ih <- ordinate(insti_healthy_phyloseq, 
                    method = "PCoA", 
                    distance = unifrac_ih)

# Plot the ordination
ih_plot <- plot_ordination(insti_healthy_phyloseq, 
                           pcoa_ih, 
                           color = "Group") + 
  theme_bw() +
  stat_ellipse(type = "norm") +
  labs(title = "Unweighted UniFrac PCoA: INSTI+ vs Healthy") +
  scale_color_manual(values = c("INSTI+" = "darkgreen", "Healthy" = "turquoise")) +
  theme(plot.title = element_text(hjust = 0.5))

# Display the plot
print(ih_plot)



##### Saving #####
save(mpt_final, file="mpt_final.RData")
save(mpt_rare, file="mpt_rare.RData")
