##### Install packages #####
# Start by installing all necessary packages when asked if you want to install
# from source, please just type Yes in the terminal below

# If you don't have BiocManager, here is the code to install it
# A lot of you probably already have this so you can skip
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Create a list of all the packages you need to install
pkgs <- c("ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "metagenomeSeq", 
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")

# Use the above list to install all the packages using a for loop
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}
# when asked if you want to update all, some or none please type "n" for none

# After installing all of its above dependencies, install ggpicrust2
install.packages("GGally")
install.packages("MicrobiomeStat")
install.packages("ggh4x")
install.packages("ggprism")
install.packages("ggpicrust2")



#### Load packages ####
# Load all necessary libraries
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(DESeq2)
library(ggh4x)
library(dplyr)

# First, remove the existing package
#remove.packages("ggpicrust2")

#BiocManager::install("ggpicrust2", force = TRUE)

install.packages("ggpicrust2_1.7.3.tar.gz", repos = NULL, type = "source")

#### Import files and preparing tables ####
# Importing the pathway abundance data from PICrust2
abundance_data <- read_delim("pathway_abundance.tsv", delim = "\t", col_names = TRUE, trim_ws = TRUE)
abundance_data <- as.data.frame(abundance_data)

# Import metadata file
metadata <- read_delim("modified_metadata_long.tsv")
metadata <- metadata[!is.na(metadata$INSTI_drug_current), ]

# Filtering to keep only HIV+HCV-, and combine INSTI columns
metadata <- metadata %>%
  mutate(INSTI = if_else(INSTI_drug_current == "YES" | INSTI_drug_prior == "YES", "YES", "NO")) %>%
  filter(hiv_status_clean == "HIV+" & hcv == "NO")

# Filtering the abundance table to only include samples in filtered metadata
sample_names <- metadata$'sample-id'
sample_names <- append(sample_names, "pathway")
abundance_data_filtered <- abundance_data[, colnames(abundance_data) %in% sample_names]

# Removing samples with all zeros
abundance_data_filtered <- abundance_data_filtered[, colSums(abundance_data_filtered != 0) > 0]

# Ensuring rownames are properly reset
rownames(abundance_data_filtered) <- NULL

# Verify samples in metadata match samples in abundance_data
abun_samples <- rownames(t(abundance_data_filtered[, -1]))
metadata <- metadata[metadata$`sample-id` %in% abun_samples, ]
metadata <- metadata %>%
  select(-current_regimen, -bdi_ii, -host_height, -host_weight)

# Prepare data for DESeq2
abundance_data_filtered <- abundance_data_filtered %>% column_to_rownames("pathway")
abundance_data_filtered_matrix <- round(as.matrix(abundance_data_filtered))

#### Pre-filter data for DESeq2 analysis ####
# Remove pathways with all zeros
abundance_filtered <- abundance_data_filtered_matrix[rowSums(abundance_data_filtered_matrix) > 0, ]

# Remove samples with all zeros
abundance_filtered <- abundance_filtered[, colSums(abundance_filtered) > 0]

# Update metadata to match filtered samples
metadata_filtered <- metadata[metadata$`sample-id` %in% colnames(abundance_filtered), ]

#### DESEq ####
#Perform pathway DAA using DESEQ2 method
abundance_daa_results_df <- pathway_daa(abundance = abundance_filtered, 
                                        metadata = metadata_filtered, group = "INSTI", daa_method = "DESeq2")


# Annotate MetaCyc pathway so they are more descriptive
metacyc_daa_annotated_results_df <- pathway_annotation(pathway = "MetaCyc", 
                                                       daa_results_df = abundance_daa_results_df, ko_to_kegg = FALSE)

# Filter p-values to only significant ones
feature_with_p_0.05 <- abundance_daa_results_df %>% filter(p_adjust < 0.001)


#Changing the pathway column to description for the results 
feature_desc = inner_join(feature_with_p_0.05,metacyc_daa_annotated_results_df, by = "feature")
feature_desc$feature = feature_desc$description
feature_desc = feature_desc[,c(1:7)]
colnames(feature_desc) = colnames(feature_with_p_0.05)

#Adding pathway column
abundance_data_filtered <- abundance_data_filtered %>% tibble::rownames_to_column("pathway")  ##JD added, double check it does not mess up for later DeSEQ

#Changing the pathway column to description for the abundance table
abundance = abundance_data_filtered %>% filter(pathway %in% feature_with_p_0.05$feature)
colnames(abundance)[1] = "feature"
abundance_desc = inner_join(abundance,metacyc_daa_annotated_results_df, by = "feature")
abundance_desc$feature = abundance_desc$description
#508 represents the number of samples in the filtered abundance table
abundance_desc = abundance_desc[,-c(508:ncol(abundance_desc))] 

# Should it be the filtered one (metadata)?
# Generate a heatmap
pathway_heatmap(abundance = abundance_desc %>% column_to_rownames("feature"), metadata = metadata, group = "INSTI")

# Generate pathway PCA plot
pathway_pca(abundance = abundance_data_filtered %>% column_to_rownames("pathway"), metadata = metadata, group = "INSTI")

# Generating a bar plot representing log2FC from the custom deseq2 function

# Go to the Deseq2 function script and update the metadata category of interest

# Lead the function in
source("DESeq2_function.R")

### Runs up to this point
# Run the function on your own data
res =  DEseq2_function(abundance_data_filtered, metadata, "INSTI")
res$feature =rownames(res)
res_desc = inner_join(res,metacyc_daa_annotated_results_df, by = "feature")
res_desc = res_desc[, -c(8:13)]
View(res_desc)

# Filter to only include significant pathways
sig_res = res_desc %>%
  filter(pvalue < 0.05) %>%
  filter(log2FoldChange > 1 | log2FoldChange < -1 )
    
# You can also filter by Log2fold change

sig_res <- sig_res[order(sig_res$log2FoldChange),]

diff_path <- ggplot(data = sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")+
  scale_y_discrete(labels = c(
    "coenzyme B biosynthesis" = "coenzyme B\nbiosynthesis",
    "methanogenesis from H2 and CO2" = "methanogenesis from \n H2 and CO2",
    "archaetidylserine and archaetidylethanolamine biosynthesis" = 
      "archaetidylserine and\narchaetidylethanolamine\nbiosynthesis")) +
  theme_minimal(base_size = 20) +
  theme(axis.text = element_text(size = 20))

diff_path


ggsave("FINAL_functional_analysis.png", 
       diff_path,
       width = 14,     
       height = 8,    
       dpi = 300,     
       units = "in",  
       limitsize = FALSE,  
       bg = "white")