library(tidyverse)
library(phyloseq)

load("longmeta_phyloseq_unrare.RData")
Mono_Infected <- subset_samples(mpt_final, hiv_status_clean == "HIV+"& hcv == "NO")

# Extract components from the phyloseq object
tax_table <- as(tax_table(Mono_Infected), "matrix")
tax_table_df <- as.data.frame(tax_table)

# Find most resolved taxonomic rows
tax_table_most_resolved <- tax_table_df %>%
  mutate(level = rowSums(!is.na(across(everything())))) %>%
  filter(level == max(.$level)) %>%
  select(-level)

# Extract OTU table
otu_table <- as(otu_table(Mono_Infected), "matrix")

# Optionally, if you want to save the processed phyloseq object
ps_processed <- phyloseq(
  sample_data(Mono_Infected),
  tax_table(as.matrix(tax_table_most_resolved)),
  otu_table(otu_table, taxa_are_rows = TRUE)
)

# Save the processed phyloseq object
saveRDS(ps_processed, "processed_phyloseq.rds")
