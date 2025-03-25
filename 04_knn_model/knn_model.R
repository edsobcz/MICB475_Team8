library(phyloseq)
library(tidyverse)
library(tidymodels)
library(microbiome)
library(metagMisc)

# loading in phyloseq objects
load("longmeta_phyloseq_rare.RData")

rarefied <- mpt_rare

indicspecies <- read_csv("../03_indicator_species/insti_indicspecies.csv") |> 
  mutate(insti = as.factor(case_when(index == 1 ~ "No",
                           index == 2 ~ "Yes")))

### CREATING A DATA FRAME WITH METADATA AND INDIC SPECIES TAXA DATA ONLY ###

# pulling a list of indicspecies with all information 
indicspecies_all <- indicspecies |> 
  unite(col = "Name", Phylum:Species, remove = FALSE) |> 
  pull(Name)
indicspecies_all

# pulling the Genus of the indicator species
indicspecies_genus <- indicspecies |> 
  filter(!is.na(Genus)) |> 
  pull(Genus)

# Normalize microbiome data (only the indicator species)
df_clr <- rarefied |>  microbiome::transform('clr') |> # CLR transformation
  subset_taxa(Genus %in% indicspecies_genus) |>  
  psmelt() |>  
  # Add Z transformation to each species individually
  group_by(Species) |> 
  mutate(Abundance = scale(Abundance)) |> 
  ungroup() |> 
  filter(hiv_status_clean == "HIV+" & hcv == "NO") |> 
  unite(col = "Name", Phylum:Species, remove = FALSE) |>
  filter(Name %in% indicspecies_all)


### RESHAPING THE DATAFRAME FOR USE IN KNN MODELS ###

# number of unique patients in the dataset after filtering for only HIV+ HCV- that have at least one indicator species
length(unique(df_clr$Sample))

# select relevant columns, mutate non applicable values to NAs
unshaped <- df_clr |> 
  select(OTU, Sample, INSTI_drug_current, diabetes, host_age, sex, host_height, host_weight, Name) |> 
  mutate(diabetes = na_if(diabetes, "not applicable"),
         sex = na_if(sex, "not provided"),
         host_weight = as.numeric(host_weight),
         host_height = as.numeric(host_height)) |> 
  mutate(bmi = host_weight/((host_height/100)^2))

# number of unique patients excluding NA values
filtered <- unshaped |> 
  filter(!is.na(diabetes) & !is.na(sex)) |> 
  select(-host_weight, -host_height)

length(unique(filtered$Sample))
# we only lose 16 patients by filtering out NAs, which is GREAT


shaped <- filtered |> 
  mutate(presence = 1) |> 
  pivot_wider(names_from = Name, values_from = presence, values_fill = 0) |> 
  group_by(Sample, INSTI_drug_current, sex, diabetes, host_age, bmi) |> 
  summarize(across(starts_with("p_"), sum)) |> 
  ungroup()

# with the filtered code, we get multiples for each species, which is because some species appear multiple times
# how do we want to approach this? we could convert it to presence/absence, but i don't know how many taxa actually have 0s

write_csv(shaped, "preliminary_discussion_data.csv")

