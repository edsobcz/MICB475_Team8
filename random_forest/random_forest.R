
library(phyloseq)
library(ape) 
library(tidyverse)
#### Create RObject from the beginning #### 

## Load data ##
# Using long metadata
meta_long <- read_delim("modified_metadata_long.tsv", delim="\t")
otu <- read_delim(file = "feature-table.txt", delim="\t", skip=1)
tax <- read_delim("taxonomy.tsv", delim="\t")
phylotree <- read.tree("tree.nwk")

## Format Files ##

# OTU table
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`#OTU ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 


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