### Loading in necessary packages and objects ###

# loading in necessary packages
library(tidyverse)
library(phyloseq)
library(indicspecies)
library(ggplot2)
library(dplyr)

# loading in the phyloseq object
load("mpt_final.RData")
Mono_Infected <- mpt_final %>% subset_samples(hiv_status_clean == "HIV+"& hcv == "NO")

### Indicator Species/Taxa Analysis ####
# Glom to Genus level 
set.seed(1234)
mpt_genus_RA <- Mono_Infected %>%
  tax_glom("Genus", NArm = FALSE) %>%
  transform_sample_counts(fun=function(x) x/sum(x))

#ISA
isa_mpt <- multipatt(t(otu_table(mpt_genus_RA)), cluster = sample_data(mpt_genus_RA)$`INSTI_drug_current`)
taxtable <- tax_table(mpt_genus_RA) %>% as.data.frame() %>% rownames_to_column(var="ASV")

# Uncomment to check summary 
# summary(isa_mpt)

# consider that your table is only going to be resolved up to the genus level, be wary of 
# anything beyond the glomed taxa level
inst_isa <- isa_mpt$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) 

# Checking INSTI+
INSTI_indic_species <- filter(inst_isa, s.YES == 1) |> 
  select(ASV:Species) |> 
  arrange(desc(stat))

group_by(INSTI_indic_species, Class) |> 
  summarize(count = n()) |> 
  ggplot(aes(x = Class, y = count)) +
  geom_bar(stat = "identity") +
  labs(title = "INSTI Indicator Species") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Checking for INSTI-
no_INSTI_indic_species <- filter(inst_isa, s.NO == 1) |> 
  select(stat:Species) |> 
  arrange(desc(stat))

group_by(no_INSTI_indic_species, Class) |> 
  summarize(count = n()) |> 
  ggplot(aes(x = Class, y = count)) +
  geom_bar(stat = "identity") +
  labs(title = "INSTI Indicator Species") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Comparing INSTI+/INSTI-
all_indic_species <- inst_isa |> 
  mutate(insti = case_when(index == 1 ~ "No INSTI",
                           index == 2 ~ "INSTI")) |> 
  select(stat:insti) |>
  arrange(desc(stat)) |> 
  group_by(Genus, insti) |> 
  summarize(count = n()) |> 
  ggplot(aes(x = Genus, y = count)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(insti)) +
  labs(title = "INSTI Indicator Species") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.55, vjust = 0.65))

all_indic_species

# Other possible visualization
# Create visualization with stat on y-axis
genus_indicators <- inst_isa %>%
  mutate(insti_group = case_when(
    s.YES == 1 ~ "INSTI+",
    s.NO == 1 ~ "INSTI -"
  )) %>%
  mutate(Genus = gsub("g__", "", Genus)) %>%  
  filter(!is.na(Genus)) %>%
  mutate(insti_group = factor(insti_group, levels = c("INSTI+", "INSTI -")))

# Create plot with facet_grid using stat values
all_indic_species2 <- ggplot(genus_indicators, aes(x = Genus, y = stat,fill = insti_group)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(insti_group), scales = "free_x") +
  labs(title = NULL,
    x = "Genus",
    y = "Indicator Value (stat)",
    fill = NULL) +
  scale_fill_manual(values = c("INSTI+" = "darkgreen", "INSTI -" = "orange")) +
  theme_bw() +
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 0.95, vjust = 0.95)) +
  scale_x_discrete(labels = c("Eubacterium",
                              "Anaerostignum",
                              "Corynebacterium",
                              "Gluconobacter",
                              "Harryflintia",
                              "Lactonifactor",
                              "Odoribacter",
                              "Pseudoflavonifractor")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(text = element_text(size = 14))

all_indic_species2

# saving indicator species as an object
write_csv(inst_isa, file = "insti_indicspecies.csv")
ggsave("Indicator_species1.png", plot = all_indic_species)
ggsave("Indicator_species2.png", plot = all_indic_species2)





