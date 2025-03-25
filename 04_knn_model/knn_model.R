library(phyloseq)
library(tidymodels)
library(microbiome)
library(metagMisc)
library(tidyverse)
set.seed(6666)

# loading in phyloseq objects
load("longmeta_phyloseq_rare.RData")

rarefied <- mpt_rare

indicspecies <- read_csv("../03_indicator_species/insti_indicspecies.csv")

insticore <- read_csv("../03_core_microbiome/insti_core_microbiome.csv") |> pull(Genus)

### CREATING A DATA FRAME WITH METADATA AND INDIC SPECIES TAXA DATA ONLY ###

# pulling a list of indicspecies with all information 
# indicspecies_all <- indicspecies |> 
#   unite(col = "Name", Phylum:Species, remove = FALSE) |> 
#   pull(Name)
# indicspecies_all

# pulling the Genus of the indicator species
indicspecies_genus <- indicspecies |> 
  filter(!is.na(Genus)) |> 
  pull(Genus)

# Normalize microbiome data -- can change subset_taxa to be only the indic species or only the core microbiome
df_clr <- rarefied |>  microbiome::transform('clr') |> # CLR transformation
  subset_taxa(Genus %in% insticore) |>  
  psmelt() |>  
  # Add Z transformation to each species individually
  group_by(Species) |> 
  mutate(Abundance = scale(Abundance)) |> 
  ungroup() |> 
  filter(hiv_status_clean == "HIV+" & hcv == "NO")


### RESHAPING THE DATAFRAME FOR USE IN KNN MODELS ###

# number of unique patients in the dataset after filtering for only HIV+ HCV- that have at least one indicator species
length(unique(df_clr$Sample))

# select relevant columns, mutate non applicable values to NAs
unshaped <- df_clr |> 
  mutate(diabetes = na_if(diabetes, "not applicable"),
         sex = na_if(sex, "not provided"),
         host_weight = as.numeric(host_weight),
         host_height = as.numeric(host_height)) |> 
  mutate(bmi = host_weight/((host_height/100)^2))

# number of unique patients excluding NA values
filtered <- unshaped |> 
  filter(!is.na(diabetes) & !is.na(sex)) |> 
  mutate(is_female = case_when(sex == "female" ~ 1,
                               sex != "female" ~ 0)) |> 
  mutate(is_diabetic = case_when(diabetes == "YES" ~ 1,
                                 diabetes != "YES" ~ 0)) |>
  select(Sample, INSTI_drug_current, is_diabetic, host_age, is_female, bmi, Genus, Abundance)

length(unique(filtered$Sample))
length(unique(filtered$Genus))
# we only lose 16 patients by filtering out NAs, which is GREAT


# we now need to combine the scaled abundance for each of the species within a genus
# this is necessary to make one abundance value for each genus for each person
cleaned <- filtered |>
  group_by(Sample, INSTI_drug_current, is_diabetic, host_age, is_female, bmi, Genus) |> 
  summarize(genus_abundance = sum(Abundance)) |> 
  ungroup() |> 
  pivot_wider(names_from = Genus, values_from = genus_abundance) |> 
  mutate(host_age = as.numeric(host_age)) |> 
  mutate(INSTI_drug_current = as.factor(INSTI_drug_current))

# save final model data
write_csv(cleaned, "cleaned_model_data.csv")


### MODELLING WITH KNN ###

set.seed(7777)

# splitting data
insti_split <- initial_split(cleaned, prop = 0.75, strata = INSTI_drug_current)
train <- training(insti_split)
test <- testing(insti_split)

# creating the recipe. normalizing the numerical variables and imputing the missing bmi with median
knn_recipe <- recipe(INSTI_drug_current ~ is_diabetic + is_female + host_age + bmi + g__Anaerostipes + g__Fusicatenibacter,
                 data = train) |> 
              step_impute_median(bmi) |> 
              step_normalize(all_predictors())

# creating the model spec, using linear distance, tuning for the value of k
knn_tune <- nearest_neighbor(weight_func = "rectangular", neighbors = tune()) |> 
  set_engine("kknn") |> 
  set_mode("classification")

# creating a cross-validation object
knn_vfold <- vfold_cv(train, v = 5, strata = INSTI_drug_current)

# creating a workflow to tune for k
knn_tune_workflow <- workflow() |> 
  add_recipe(knn_recipe) |> 
  add_model(knn_tune) |> 
  tune_grid(resamples = knn_vfold, grid = tibble(neighbors = c(1:20)))

# looking at the best k
best_k <- knn_tune_workflow |> 
  collect_metrics() |> 
  filter(.metric == "roc_auc" | .metric == "accuracy")

best_val <- ggplot(best_k, aes(x = neighbors, y = mean)) +
  geom_line() +
  facet_grid(rows = vars(.metric))

best_val  
# based on the plot, it looks like the best roc/auc is at 5, which is also where the 
# accuracy plateaus a bit. i think that going above 10 is overfitting

# creating a new model spec with the best k
knn_spec <- nearest_neighbor(weight_func = "rectangular", neighbors = 5) |> 
  set_engine("kknn") |> 
  set_mode("classification")

# creating a new workflow with the best k
knn_fit <- workflow() |> 
  add_recipe(knn_recipe) |> 
  add_model(knn_spec) |> 
  fit(data = train)

# predicting
predictions <- predict(knn_fit, test) |> bind_cols(test)

# analyzing model accuracy
class_metrics <- metric_set(accuracy)
metrics <- class_metrics(predictions, truth = INSTI_drug_current, estimate = .pred_class)

confusion_matrix <- predictions |> 
  conf_mat(truth = INSTI_drug_current, estimate = .pred_class)

confusion_matrix

# analyzing the class distribution of the test set
classes <- test |> 
  group_by(INSTI_drug_current) |> 
  summarize(count = n(),
            prop = count/nrow(test))
classes

### CONCLUSION ###

# our classifier has an accuracy of 0.54. just picking the majority class gives us an accuracy
# of 0.66. therefore our knn classifier (unsurprisingly) sucks