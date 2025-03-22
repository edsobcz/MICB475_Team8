### Loading in necessary packages and objects ###

# loading in necessary packages
library(tidyverse)
library(phyloseq)
library(indicspecies)

# loading in the phyloseq object
load("mpt_rare.RData")
# Uncomment the following line to load Monoinfected mpt final with updated file path
# load("mpt_final_monoinfected.RData")

set.seed(1)

### Perform indicator species analysis ###
Mono_Infected <- subset_samples(mpt_final, hiv_status_clean == "HIV+"& hcv == "NO")
insti_relative <- transform_sample_counts(mpt_rare, fun=function(x) x/sum(x))

isa_output <- multipatt(t(otu_table(insti_relative)), 
                        cluster = sample_data(insti_relative)$INSTI_drug_current)


### Data summary ###
insti_summary <- summary(isa_output)

# full output
taxtable <- tax_table(mpt_rare) %>% as.data.frame() %>% rownames_to_column(var="ASV")

insti_isa <- isa_output$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05)

insti_isa


### Cleaning up output ###
# for INSTI
INSTI_indic_species <- filter(insti_isa, s.YES == 1) |> 
  select(ASV:Species) |> 
  arrange(desc(stat))

I

summary(INSTI_indic_species)

group_by(INSTI_indic_species, Class) |> 
  summarize(count = n()) |> 
  ggplot(aes(x = Class, y = count)) +
    geom_bar(stat = "identity") +
    labs(title = "INSTI Indicator Species") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

# for no INSTI
no_INSTI_indic_species <- filter(insti_isa, s.NO == 1) |> 
  select(stat:Species) |> 
  arrange(desc(stat))

summary(no_INSTI_indic_species)

group_by(no_INSTI_indic_species, Class) |> 
  summarize(count = n()) |> 
  ggplot(aes(x = Class, y = count)) +
  geom_bar(stat = "identity") +
  labs(title = "INSTI Indicator Species") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# together
all_indic_species <- insti_isa |> 
  mutate(insti = case_when(index == 1 ~ "no INSTI",
                           index == 2 ~ "INSTI")) |> 
  select(stat:insti) |>
  arrange(desc(stat)) |> 
  group_by(Class, insti) |> 
  summarize(count = n()) |> 
  ggplot(aes(x = Class, y = count)) +
    geom_bar(stat = "identity") +
    facet_grid(rows = vars(insti)) +
    labs(title = "INSTI Indicator Species") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

all_indic_species

#### Updated: SEBASTIAN Indicator Species/Taxa Analysis ####
# Glom to Species level to prevent loose of significant ASVs
set.seed(1234)
mpt_genus <- tax_glom(Mono_Infected, "Species", NArm = FALSE)
mpt_genus_RA <- transform_sample_counts(mpt_genus, fun=function(x) x/sum(x))

#ISA
isa_mpt <- multipatt(t(otu_table(mpt_genus_RA)), cluster = sample_data(mpt_genus_RA)$`INSTI_drug_current`)
summary(isa_mpt)
taxtable <- tax_table(mpt_final) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# consider that your table is only going to be resolved up to the genus level, be wary of 
# anything beyond the glomed taxa level
inst_isa <- isa_mpt$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) 

# together
all_indic_species <- inst_isa |> 
  mutate(insti = case_when(index == 1 ~ "no INSTI",
                           index == 2 ~ "INSTI")) |> 
  select(stat:insti) |>
  arrange(desc(stat)) |> 
  group_by(Class, insti) |> 
  summarize(count = n()) |> 
  ggplot(aes(x = Class, y = count)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(insti)) +
  labs(title = "INSTI Indicator Species") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

all_indic_species

