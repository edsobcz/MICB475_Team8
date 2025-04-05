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
# meta <- read_delim("modified_metadata.tsv", delim="\t")
otu <- read_delim(file = "feature-table.txt", delim="\t", skip=1)
tax <- read_delim("taxonomy.tsv", delim="\t")
phylotree <- read.tree("tree.nwk")

# Metadata file with 7 more columns
meta_long <- read_delim("modified_metadata_long.tsv", delim="\t")
## Format Files ##

# OTU table
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 

#Metadata 
#samp_df <- as.data.frame(meta[,-1])
#rownames(samp_df)<- meta$'sample-id'
#SAMP <- sample_data(samp_df)

#Metadata Long
samp_df <- as.data.frame(meta_long[,-1])
rownames(samp_df)<- meta_long$'sample-id'
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
set.seed(1234)

# Remove non-bacterial sequences, if any
mpt_filt <- subset_taxa(mpt,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
mpt_filt_nolow <- filter_taxa(mpt_filt, function(x) sum(x)>5, prune = TRUE)
# Remove Blanks from the dataset 
mpt_no_blanks <- subset_samples(mpt_filt_nolow, sample_type != "control blank") 
# Remove samples with less than 100 reads
mpt_final <- prune_samples(sample_sums(mpt_no_blanks)>100, mpt_no_blanks) 
mpt_rare <- rarefy_even_depth(mpt_final, rngseed = 124, sample.size = 25228)
mpt_final_monoinfected <- subset_samples(mpt_final, hiv_status_clean == "HIV+"& hcv == "NO")

# Subset samples for subsequent analysis
Infected <- subset_samples(mpt_rare, hiv_status_clean == "HIV+"& hcv == "NO")
Healthy <- subset_samples(mpt_rare, hiv_status_clean == "HIV-"& hcv == "NO")
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
  group = c("HIV mono-infected", "HIV-/HCV-", "HCV infected"),
  count = c(count_HIVpos, count_healthy, count_remaining)
)

# Calculate percentages
pie_data$percentage <- round(pie_data$count / sum(pie_data$count) * 100, 1)
pie_data$label <- paste0(pie_data$group, "\n", pie_data$percentage, "%")


# Create pie chart 
pie_chart_types <- ggplot(pie_data, aes(x = "", y = count, fill = group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("HIV mono-infected" = "firebrick", 
                               "HIV-/HCV-" = "steelblue", 
                               "HCV infected" = "gold3")) +
  #labs(title = "Distribution of Sample Types") +
  geom_text(aes(label = label, 
                hjust = ifelse(group == "HIV mono-infected", 0.5, 0.6)), 
            position = position_stack(vjust = 0.5), 
            size = 5.2) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

pie_chart_types

# ----- HIV+ INSTI + vs HIV+ INSTI- -----
# Count the number of samples in each group
count_INSTIpos <- nsamples(INSTI_pos)
count_INSTIneg <- nsamples(INSTI_neg)

# Create a data frame for plotting
pie_data_INSTI <- data.frame(
                    group = c("INSTI+", "INSTI-"),
                    count = c(count_INSTIpos,count_INSTIneg))

# Calculate percentages
pie_data_INSTI$percentage <- round(pie_data_INSTI$count / sum(pie_data_INSTI$count) * 100, 1)
pie_data_INSTI$label <- paste0(pie_data_INSTI$group, "\n", pie_data_INSTI$percentage, "%")


# Create pie chart 
pie_chart_INSTI <- ggplot(pie_data_INSTI, aes(x = "", y = count, fill = group)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = c("INSTI+" = "darkgreen", 
                               "INSTI-" = "darkorange")) +
  geom_text(aes(label = label, 
                vjust = ifelse(group == "INSTI-", 1.1, 
                               ifelse(group == "INSTI+", -0.2, 0.5))),
            position = position_stack(vjust = 0.5), 
            size = 5.2) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))
pie_chart_INSTI

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
rich_healthy$group <- "HIV-/HCV-"

# Add sample names as a column
rich_HIVpos$sample <- rownames(rich_HIVpos)
rich_healthy$sample <- rownames(rich_healthy)

# Combine the data
combined_richness <- rbind(rich_HIVpos, rich_healthy)

# Plot
shannon_hh<- ggplot(combined_richness, aes(x = group, y = Shannon, fill = group)) +
              geom_boxplot() +
              theme_bw() +
              labs(y = "Shannon Diversity",
                   x = "") +
              scale_fill_manual(values = c("HIV+" = "firebrick", "HIV-/HCV-" = "steelblue"))+
              theme(legend.position = "none",
                    axis.text.x = element_text(size = 15),
                    axis.text.y = element_text(size = 10),
                    axis.title.y = element_text(size = 15)) 

shannon_hh

#Wilcox Test
wilcox.test(rich_HIVpos$Shannon, rich_healthy$Shannon, exact = FALSE)

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
# wilcox.test(pd_infected$PD, pd_healthy$PD)

# Create a boxplot to visualize
pd_combined <- data.frame(
  Group = c(rep("HIV+", nrow(pd_infected)), rep("HIV-/HCV-", nrow(pd_healthy))),
  Faith_PD = c(pd_infected$PD, pd_healthy$PD)
)

faiths_hh <- ggplot(pd_combined, aes(x = Group, y = Faith_PD, fill = Group)) +
              geom_boxplot() +
              theme_bw() +
              labs(title = "Faith's Phylogenetic Diversity by HIV Status",
                   y = "Faith's PD",
                   x = "") +
              scale_fill_manual(values = c("HIV+" = "firebrick", "HIV-/HCV-" = "steelblue"))+
              theme(legend.position = "none",
                    axis.text.x = element_text(size = 15),
                    axis.text.y = element_text(size = 10),
                    axis.title.y = element_text(size = 15))
faiths_hh

# ----- Shannon Diversity HIV+ dep+/HIV+ dep-  --------
# Calculate alpha diversity with depression into account
depression <- subset_samples(Infected, bdi_ii > 20)
no_depression <- subset_samples(Infected, bdi_ii <= 20)

rich_dep <- estimate_richness(depression, measures = c("Shannon"))
rich_dep_no <- estimate_richness(no_depression, measures = c("Shannon"))

# Add group labels
rich_dep$group <- "Depressed"
rich_dep_no$group <- "Non-depressed"

# Add sample names as a column
rich_dep$sample <- rownames(rich_dep)
rich_dep_no$sample <- rownames(rich_dep_no)

# Combine the data
combined_richness <- rbind(rich_dep, rich_dep_no)

# Plot
shannon_hdp <- ggplot(combined_richness, aes(x = group, y = Shannon, fill = group)) +
                geom_boxplot() +
                theme_bw() +
                labs(y = "Shannon Diversity",
                     x = "") +
                scale_fill_manual(values = c("Depressed" = "tan", "Non-depressed" = "mediumorchid4"))+
                theme(legend.position = "none",
                      axis.text.x = element_text(size = 15),
                      axis.text.y = element_text(size = 10),
                      axis.title.y = element_text(size = 15))
shannon_hdp

#Wilcox Test
wilcox.test(rich_dep$Shannon, rich_dep_no$Shannon)

#----- Faiths Phylogenetic HIV+ dep+/HIV+ dep- ------
# Extract OTU tables and tree
otu_table_depressed <- otu_table(depression)
otu_table_no_depressed <- otu_table(no_depression)
tree <- phy_tree(Infected)

# If taxa are rows in your OTU table, transpose them
if(taxa_are_rows(otu_table_depressed)) {
  otu_table_depressed <- t(otu_table_depressed)
}
if(taxa_are_rows(otu_table_no_depressed)) {
  otu_table_no_depressed <- t(otu_table_no_depressed)
}

# Calculate Faith's PD
pd_depressed <- pd(otu_table_depressed, tree)
pd_no_depressed <- pd(otu_table_no_depressed, tree)

# Wilcoxon test
#wilcox.test(pd_depressed$PD, pd_no_depressed$PD)

# Create a boxplot to visualize Faith's PD for depression status
pd_combined <- data.frame(
  Group = c(rep("Depressed", nrow(pd_depressed)), rep("Non-depressed", nrow(pd_no_depressed))),
  Faith_PD = c(pd_depressed$PD, pd_no_depressed$PD)
)

# Plot the boxplot with appropriate title and colors
faiths_hdp <- ggplot(pd_combined, aes(x = Group, y = Faith_PD, fill = Group)) +
                geom_boxplot() +
                theme_bw() +
                labs(title = "Faith's Phylogenetic Diversity by Depression \n Status in HIV+ Individuals",
                     y = "Faith's PD",
                     x = "") +
                scale_fill_manual(values = c("Depressed" = "tan", "Non-depressed" = "chocolate"))+
                theme(legend.position = "none",
                      plot.title = element_text(hjust = 0.5))
faiths_hdp

# ----- Shannon Diversity INSTI+/INSTI-  --------
# Calculate alpha diversity with depression into account
rich_insti <- estimate_richness(INSTI_pos, measures = c("Shannon"))
rich_no_insti <- estimate_richness(INSTI_neg, measures = c("Shannon"))

# Add group labels
rich_insti$group <- "INSTI+"
rich_no_insti$group <- "INSTI-"

# Add sample names as a column
rich_insti$sample <- rownames(rich_insti)
rich_no_insti$sample <- rownames(rich_no_insti)

# Combine the data
combined_richness <- rbind(rich_insti, rich_no_insti)

# Plot
shannon_insti <- ggplot(combined_richness, aes(x = group, y = Shannon, fill = group)) +
  geom_boxplot() +
  theme_bw() +
  labs(y = "Shannon Diversity",
       x = "") +
  scale_fill_manual(values = c("INSTI+" = "darkgreen", "INSTI-" = "darkorange"))+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 15))

shannon_insti

#Wilcox Test
wilcox.test(rich_insti$Shannon, rich_no_insti$Shannon)

#----- Faiths Phylogenetic INSTI+/INSTI- ------
# Extract OTU tables and tree
otu_table_insti <- otu_table(INSTI_pos)
otu_table_no_insti <- otu_table(INSTI_neg)
tree <- phy_tree(Infected)

# If taxa are rows in your OTU table, transpose them
if(taxa_are_rows(otu_table_insti)) {
  otu_table_insti <- t(otu_table_insti)
}
if(taxa_are_rows(otu_table_no_insti)) {
  otu_table_no_insti <- t(otu_table_no_insti)
}

# Calculate Faith's PD
pd_insti <- pd(otu_table_insti, tree)
pd_no_insti <- pd(otu_table_no_insti, tree)

# Wilcoxon test
#wilcox.test(pd_insti$PD, pd_no_insti$PD)

# Create a boxplot to visualize Faith's PD for depression status
pd_combined <- data.frame(
  Group = c(rep("INSTI+", nrow(pd_insti)), rep("INSTI-", nrow(pd_no_insti))),
  Faith_PD = c(pd_insti$PD, pd_no_insti$PD)
)

# Plot the boxplot with appropriate title and colors
faiths_insti <- ggplot(pd_combined, aes(x = Group, y = Faith_PD, fill = Group)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Faith's Phylogenetic Diversity by INSTI Status",
       y = "Faith's PD",
       x = "") +
  scale_fill_manual(values = c("INSTI+" = "darkblue", "INSTI-" = "darkgrey"))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

faiths_insti

#### Beta Diversity ####
# Calculate distances
unifrac_dist <- phyloseq::distance(mpt_rare, method = "unifrac", weighted = FALSE)
weighted_unifrac_dist <- phyloseq::distance(mpt_rare, method = "unifrac", weighted = TRUE)

# This are for testing different mattices for pcoa (So far best variance are Weighted unifrac and Bray)
bray_dist <- phyloseq::distance(mpt_rare, method = "bray")
jaccard_dist <- phyloseq::distance(mpt_rare, method = "jaccard")

bray_dist_hh <- as.dist(as.matrix(bray_dist)[hiv_healthy_samples, hiv_healthy_samples])
jaccard_dist_hh <- as.dist(as.matrix(jaccard_dist)[hiv_healthy_samples, hiv_healthy_samples])

# ---- Unweighted Unifrac HIV+/Healthy -----
# Create phyloseq object first
hiv_healthy_phyloseq <- subset_samples(mpt_rare, 
                                       hiv_status_clean == "HIV+" | 
                                      (hiv_status_clean == "HIV-" & hcv == "NO"))

# Extract sample names
hiv_healthy_samples <- sample_names(hiv_healthy_phyloseq)

# Subset the distance matrix to include only these samples
unifrac_dist_subset_hh <- as.dist(as.matrix(unifrac_dist)[hiv_healthy_samples, hiv_healthy_samples])

# Add a clear group label for plotting
sample_data(hiv_healthy_phyloseq)$Group <- ifelse(
  sample_data(hiv_healthy_phyloseq)$hiv_status_clean == "HIV+", 
  "HIV+", 
  "HIV-/HCV-"
)

# Run the ordination
pcoa_hh <- ordinate(hiv_healthy_phyloseq, 
                    method = "PCoA", 
                    distance = unifrac_dist_subset_hh)

# Plot the ordination
pcoa_hh_plot <- plot_ordination(hiv_healthy_phyloseq, 
                                pcoa_hh, 
                                color = "Group") + 
                  theme_bw() +
                  stat_ellipse(type = "norm") +
                  labs(title = "Unweighted UniFrac PCoA by HIV Status") +
                  scale_color_manual(values = c("HIV+" = "firebrick", "HIV-/HCV-" = "turquoise")) +
                  theme(plot.title = element_text(hjust = 0.5))

pcoa_hh_plot

# Run PERMANOVA
# Get the metadata for these samples
meta_subset_hh <- data.frame(sample_data(mpt_rare))[hiv_healthy_samples, ]

# Create a factor for grouping
meta_subset_hh$group_factor <- ifelse(
  meta_subset_hh$hiv_status_clean == "HIV-" & meta_subset_hh$hcv == "NO", 
  "HIV-/HCV-", 
  "HIV+"
)

# PERMANOVA test
permanova_result_hh <- adonis2(unifrac_dist_subset_hh ~ group_factor, data = meta_subset_hh)
print(permanova_result_hh)


# ---- Unweighted Unifrac HIV+ dep+/HIV+ dep- -----
# Create phyloseq object first
hiv_phyloseq <- subset_samples(mpt_rare,hiv_status_clean == "HIV+" & hcv == "NO" & !is.na(bdi_ii))

# Extract sample names
hiv_samples <- sample_names(hiv_phyloseq)

# Subset the distance matrix to include only these samples
unifrac_dist_subset_hdp <- as.dist(as.matrix(unifrac_dist)[hiv_samples, hiv_samples])

# Add a clear group label for plotting
sample_data(hiv_phyloseq)$Group <- ifelse(
  sample_data(hiv_phyloseq)$hiv_status_clean == "HIV+" & 
    sample_data(hiv_phyloseq)$bdi_ii > 20, 
  "Depressed", 
  "Non-depressed"
)

# Run the ordination
pcoa_hdp <- ordinate(hiv_phyloseq, 
                    method = "PCoA", 
                    distance = unifrac_dist_subset_hdp)

# Plot the ordination
pcoa_hdp_plot <- plot_ordination(hiv_phyloseq, 
                                pcoa_hdp, 
                                color = "Group") + 
                  theme_bw() +
                  stat_ellipse(type = "norm") +
                  labs(title = "Unweighted UniFrac PCoA by \n Depression Status in HIV+ Individuals") +
                  scale_color_manual(values = c("Depressed" = "tan", "Non-depressed" = "chocolate")) +
                  theme(plot.title = element_text(hjust = 0.5))


pcoa_hdp_plot

# Run PERMANOVA
# Get the metadata for these samples
meta_subset_hdp <- data.frame(sample_data(mpt_rare))[hiv_samples, ]

# Create a factor for grouping
meta_subset_hdp$group_factor <- ifelse(
  meta_subset_hdp$bdi_ii > 20, 
  "Depressed", 
  "Non-depressed"
)

# PERMANOVA test
permanova_result_hdp <- adonis2(unifrac_dist_subset_hdp ~ group_factor, data = meta_subset_hdp)
print(permanova_result_hdp)


# ---- Unweighted Unifrac INSTI+/INSTI- -----
# Create a subset object first
insti_phyloseq <- subset_samples(mpt_rare, 
                                 (hiv_status_clean == "HIV+" & INSTI_drug_current == "YES")| 
                                   (hiv_status_clean == "HIV+" & INSTI_drug_current == "NO"))

# Extract sample names for INSTI+ and INSTI- groups
insti_samples <- sample_names(insti_phyloseq)

# Subset the distance matrix
unifrac_dist_subset_insti <- as.dist(as.matrix(unifrac_dist)[insti_samples, insti_samples])

# Add a clear INSTI label for plotting
sample_data(insti_phyloseq)$INSTI_status <- ifelse(
  sample_data(insti_phyloseq)$INSTI_drug_current == "YES", 
  "INSTI+", 
  "INSTI-"
)

# Create the ordination
pcoa_insti <- ordinate(insti_phyloseq, 
                       method = "PCoA", 
                       distance = unifrac_dist_subset_insti)

# Plot the ordination
pcoa_insti_plot <- plot_ordination(insti_phyloseq, 
                              pcoa_insti, 
                              color = "INSTI_status") + 
                    theme_bw() +
                    stat_ellipse(type = "norm") +
                    labs(title = "Unweighted UniFrac PCoA by INSTI Status in \nHIV+ Individuals") +
                    scale_color_manual(values = c("INSTI+" = "darkgreen", "INSTI-" = "darkorange"))+
                    theme(plot.title = element_text(hjust = 0.5))


pcoa_insti_plot

# Run PERMANOVA
# Get metadata for these samples
meta_subset_insti <- data.frame(sample_data(mpt_rare))[insti_samples, ]

# Create a factor for INSTI status directly from the metadata
meta_subset_insti$insti_factor <- ifelse(meta_subset_insti$INSTI_drug_current == "YES", "INSTI+", "INSTI-")

# PERMANOVA test

permanova_result_insti <- adonis2(unifrac_dist_subset_insti ~ insti_factor, data = meta_subset_insti)
print(permanova_result_insti)


# ---- Unweighted Unifrac INSTI+/Healthy -----
# Create a phyloseq subset first
insti_healthy_phyloseq <- subset_samples(mpt_rare, 
                                         (hiv_status_clean == "HIV+" & INSTI_drug_current == "YES")| 
                                           (hiv_status_clean == "HIV-" & hcv == "NO"))

# Extract sample names
insti_healthy_samples <- sample_names(insti_healthy_phyloseq)

# Subset the distance matrix to include only these samples
unifrac_dist_subset_ih <- as.dist(as.matrix(unifrac_dist)[insti_healthy_samples, insti_healthy_samples])

# Add a clear group label for plotting
sample_data(insti_healthy_phyloseq)$Group <- ifelse(
  sample_data(insti_healthy_phyloseq)$hiv_status_clean == "HIV+" & 
    sample_data(insti_healthy_phyloseq)$INSTI_drug_current == "YES", 
  "INSTI+", 
  "HIV-/HCV-"
)

# Create the ordination
pcoa_ih <- ordinate(insti_healthy_phyloseq, 
                    method = "PCoA", 
                    distance = unifrac_dist_subset_ih)

# Plot the ordination
pcoa_ih_plot <- plot_ordination(insti_healthy_phyloseq, 
                           pcoa_ih, 
                           color = "Group") + 
            theme_bw() +
            stat_ellipse(type = "norm") +
            labs(title = "Unweighted UniFrac PCoA: INSTI+ vs HIV-/HCV-") +
            scale_color_manual(values = c("INSTI+" = "darkgreen", "HIV-/HCV-" = "turquoise")) +
            theme(plot.title = element_text(hjust = 0.5))

pcoa_ih_plot

# Run PERMANOVA
# Get the metadata for these samples
meta_subset_ih <- data.frame(sample_data(mpt_rare))[insti_healthy_samples, ]

# Create a factor for grouping
meta_subset_ih$group_factor <- ifelse(
  meta_subset_ih$hiv_status_clean == "HIV+" & meta_subset_ih$INSTI_drug_current == "YES", 
  "INSTI+", 
  "HIV-/HCV-"
)

# PERMANOVA test
permanova_result_ih <- adonis2(unifrac_dist_subset_ih ~ group_factor, data = meta_subset_ih)
print(permanova_result_ih)


# ---- Weighted Unifrac HIV+/Healthy -----
# Subset the distance matrix to include only these samples
weighted_unifrac_dist_hh <- as.dist(as.matrix(weighted_unifrac_dist)[hiv_healthy_samples, hiv_healthy_samples])

# Run the ordination
pcoa_hh_w <- ordinate(hiv_healthy_phyloseq, 
                    method = "PCoA", 
                    distance = weighted_unifrac_dist_hh)

# Plot the ordination
pcoa_hh_w_plot <- plot_ordination(hiv_healthy_phyloseq, 
                                pcoa_hh_w, 
                                color = "Group") + 
                  theme_bw() +
                  stat_ellipse(type = "norm") +
                  scale_color_manual(values = c("HIV+" = "firebrick", "HIV-/HCV-" = "steelblue")) +
                  theme( axis.text.x = element_text(size = 10),
                         axis.text.y = element_text(size = 10),
                         axis.title.y = element_text(size = 15),
                         axis.title.x = element_text(size = 15),
                         legend.text = element_text(size = 13),
                         legend.title = element_text(size = 14))

pcoa_hh_w_plot

# Run PERMANOVA
set.seed(123)
permanova_result_hh_2 <- adonis2(weighted_unifrac_dist_hh ~ group_factor, data = meta_subset_hh)
print(permanova_result_hh_2)


# ---- Weighted Unifrac HIV+ dep+/HIV+ dep- -----
# Subset the distance matrix to include only these samples
weighted_unifrac_dist_subset_hdp <- as.dist(as.matrix(weighted_unifrac_dist)[hiv_samples, hiv_samples])

# Run the ordination
pcoa_hdp_w <- ordinate(hiv_phyloseq, 
                     method = "PCoA", 
                     distance = weighted_unifrac_dist_subset_hdp)

# Plot the ordination
pcoa_hdp_w_plot <- plot_ordination(hiv_phyloseq, 
                                  pcoa_hdp_w, 
                                 color = "Group") + 
                    theme_bw() +
                    stat_ellipse(type = "norm") +
                    scale_color_manual(values = c("Depressed" = "tan", "Non-depressed" = "mediumorchid4")) +
                    theme(axis.text.x = element_text(size = 10),
                          axis.text.y = element_text(size = 10),
                          axis.title.y = element_text(size = 15),
                          axis.title.x = element_text(size = 15),
                          legend.text = element_text(size = 13),
                          legend.title = element_text(size = 14))

pcoa_hdp_w_plot

# Run PERMANOVA
set.seed(123)
permanova_result_hdp_2 <- adonis2(weighted_unifrac_dist_subset_hdp ~ group_factor, data = meta_subset_hdp)
print(permanova_result_hdp_2)

# ---- Weighted Unifrac INSTI+/INSTI- -----
# Subset the distance matrix
weighted_unifrac_dist_subset_insti <- as.dist(as.matrix(weighted_unifrac_dist)[insti_samples, insti_samples])

# Create the ordination
pcoa_insti_w <- ordinate(insti_phyloseq, 
                       method = "PCoA", 
                       distance = weighted_unifrac_dist_subset_insti)

# Plot the ordination
pcoa_insti_w_plot <- plot_ordination(insti_phyloseq, 
                                   pcoa_insti_w, 
                                   color = "INSTI_status") + 
                theme_bw() +
                stat_ellipse(type = "norm") +
                scale_color_manual(values = c("INSTI+" = "darkgreen", "INSTI-" = "darkorange"))+
                theme(axis.text.x = element_text(size = 10),
                      axis.text.y = element_text(size = 10),
                      axis.title.y = element_text(size = 15),
                      axis.title.x = element_text(size = 15),
                      legend.text = element_text(size = 13),
                      legend.title = element_text(size = 14))


pcoa_insti_w_plot

# Run PERMANOVA
set.seed(123)
permanova_result_insti_2 <- adonis2(weighted_unifrac_dist_subset_insti ~ insti_factor, data = meta_subset_insti)
print(permanova_result_insti_2)


##### Saving #####
save(mpt_final, file="../01_phyloseq_objects/longmeta_phyloseq_unrare.RData")
save(mpt_rare, file="../01_phyloseq_objects/longmeta_phyloseq_rare.RData")
save(mpt_final_monoinfected, file="mpt_final_monoinfected.RData")
ggsave("Pie_chart_types.png", plot = pie_chart_types)
ggsave("Pie_chart_INSTI.png", plot = pie_chart_INSTI)
ggsave("HH_Shannon.png", plot = shannon_hh)
ggsave("HH_Faiths.png", plot = faiths_hh)
ggsave("HDP_Shannon.png", plot = shannon_hdp)
ggsave("HDP_Faiths.png", plot = faiths_hdp)
ggsave("INSTI_Shannon.png", plot = shannon_insti)
ggsave("INSTI_Faiths.png", plot = faiths_insti)
ggsave("HH_unweighted_pcoa.png", plot = pcoa_hh_plot)
ggsave("HDP_unweighted_pcoa.png", plot = pcoa_hdp_plot)
ggsave("INSTI_unweighted_pcoa.png", plot = pcoa_insti_plot)
ggsave("IH_unweighted_pcoa.png", plot = pcoa_ih_plot)
ggsave("HH_weighted_pcoa.png", plot = pcoa_hh_w_plot)
ggsave("HDP_weighted_pcoa.png", plot = pcoa_hdp_w_plot)
ggsave("INSTI_weighted_pcoa.png", plot = pcoa_insti_w_plot)
