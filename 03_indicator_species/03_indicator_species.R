### Loading in necessary packages and objects ###

# loading in necessary packages
library(tidyverse)
library(phyloseq)
library(indicspecies)

# loading in the phyloseq object
load("mpt_rare.RData")


### Perform indicator species analysis ###
insti_relative <- transform_sample_counts(mpt_rare, fun=function(x) x/sum(x))
set.seed(123)

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
  select(stat:Species) |> 
  arrange(desc(stat))

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
    xlab("Bacterial Class") +
    ylab("Count") +
    scale_x_discrete(labels = c("c__Actinobacteria" = "Actinobacteria", 
                                "c__Alphaproteobacteria" = "Alphaproteobacteria", 
                                "c__Bacilli" = "Bacilli", "c__Bacteroidia" = "Bacteroidia", 
                                "c__Clostridia" = "Clostridia", "c__Desulfovibrionia" = "Desulfovibrionia", 
                                "c__Gammaproteobacteria" = "Gammaproteobacteria", "c__Negativicutes" = "Negativicutes")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11))

all_indic_species

ggsave(filename = 'indic_plot.png',
       all_indic_species,
       height = 6, width = 5)
