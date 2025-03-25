
# Load all packages

library(ggplot2)
library(phyloseq)
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
ps <- phyloseq(OTU, SAMP, TAX, phylotree)

### Random Forest Model
# Part 1: Load and format data

# Load the example dataset.
# This phyloseq object can be created by running Random Forest - Create Dataset.Rmd.
#ps =  %>% 
  #tax_glom('Species') # Does not work because there are NAs
```

For Random Forest, we need the following data characteristics:
  
  -   Our outcome variable: here we will use disease (PD) vs. Control.

-   We should only include taxonomy and variables that we want to use in the
model.

-   Not usually ideal to test every single variable - chances of overfitting
or underfitting is high, especially if the dataset isn't very large

    -   Common to use RF to validate that differentially abundant taxa/indicator
        taxa can collectively predict the outcome variable

-   Our data can't contain missing values.

-   The variables should be normalized. Continuous data should all be on roughly
the same scale.

-   Microbiome data will first be CLR-transformed (similar to log
                                                   transformation).

-   Continuous variables will then be Z-transformed (mean of zero, standard
                                                     deviation of 1).

```{r Format dataset for RF}


# Calculate average abundance:
require(microbiome)
avg_abundance = ps %>%
  microbiome::transform('compositional') %>% # relative abundance
  psmelt() %>% 
  group_by(Species) %>% # Perform the next line separately for each Species value
  summarize(Abundance = mean(Abundance)) %>%
  ungroup() # Stop perfoming subsequent code separately for each Species

keep_taxa <- avg_abundance %>% 
  filter(Abundance>0.001) %>% 
  pull(Species) %>% # Get Species column as a vector

# Normalize microbiome data
df_clr <- ps %>% microbiome::transform('clr') %>%  # CLR transformation
  subset_taxa(Species %in% keep_taxa) %>% 
  psmelt() %>% 
  # Add Z transformation to each species individually
  group_by(Species) %>% 
  mutate(Abundance = scale(Abundance)) %>% 
  ungroup()

### We edited up until this point - but also we don't know if we did the last part right
# Select only the metadata of interest
df_filt = df_clr %>% 
  select(Sample,Case_status,Age_at_collection,Sex,BMI,Genus,Abundance) %>% 
  # Turn each genus into its own column
  pivot_wider(names_from = Genus, values_from = Abundance)

# Remove rows with NA values in the metadata
df_filt = df_filt %>% na.omit()

# Test if the continuous metadata need additional transformations to be normalized
hist(df_filt$Age_at_collection) # Normally distributed
hist(df_filt$BMI) # Reasonably normal

# Our final dataset will have the sample ID, outcome variable of interest, possible predictive variables from the metadata, and then the microbiome data.

# Z-transform any numerical data
for_RF = df_filt %>% mutate_if(is.numeric,~as.numeric(scale(.)))
for_RF_nosampleID = for_RF %>% select(-Sample) # Not necessary for RF
```

# Part 2: Run RF

For your project, think about the following approaches. What unique information
do you get from each of them, or by contrasting different models?
  
  -   Only metadata variables

-   Only microbiome variables

-   Metadata and microbiome variables together

**We will run metadata & taxonomy variables together in a single model.**
  
  ```{r Prep Data}
# Divide the table into predictors and outcome.
X <- for_RF_nosampleID %>% select(-Case_status)
y <- for_RF_nosampleID %>% pull(Case_status)

# Convert y to a factor. If continuous, the order should be c(control,case).
y = factor(y,levels = c('Control','PD'))

# Create empty vectors for storing model information.
# These will be populated with the stats from each cross-validation.
train_auc_scores <- c()
test_auc_scores <- c()

all_labels_train = list()
all_labels_test = list()
feature_importance_values <- list()

# K-fold cross validation:
# Randomly subsets the rows into 10 equal bins.
number_of_folds = 10
set.seed(421) # Keep the randomness reproducible
folds <- createFolds(y, k = number_of_folds, list = TRUE)
```

```{r RUN RF ON EACH FOLD}
for (fold in folds) {
  # Create train and test datasets.
  X_train_fold <- X[-fold, ]
  y_train_fold <- y[-fold]
  X_test_fold <- X[fold, ]
  y_test_fold <- y[fold]
  
  # Hyperparameter tuning setup
  # Give each one a few possible values (look up recommended value ranges)
  # mtry: number of variables that will be used per forest. 
  #       High = overfitting, low = uninformative
  # Splitrule: how decisions are made (see notes below!!)
  # Related to tree complexity. Larger = simpler tree
  tune_grid <- expand.grid(mtry = c(3, 5, 7, 10), 
                           # gini or extratrees for boolean outcomes, 
                           # varaince for continuous.
                           # DELETE THE OTHERS.
                           splitrule = c("gini"),#, "extratrees","variance"),
                           min.node.size = c(10, 15,  20))
  
  # This will tell the RF command how to perform the RF.
  train_control <- trainControl(method = "cv", # K-fold cross validation
                                number = number_of_folds, # 10 folds
                                classProbs = TRUE, # Predicted class probabilities 
                                # are returned instead of just class labels.
                                summaryFunction = twoClassSummary) # compute AUC
  
  # Use hyperparameter tuning to optimize each parameter.
  # Note that optimal settings are chosen based on ROC/AUC - prone to overfitting!
  set.seed(421) # Reproducible randomness
  rf_model <- train(X_train_fold, y_train_fold, # training dataset
                    method = "ranger",
                    trControl = train_control, # Perform tuning
                    tuneGrid = tune_grid,
                    metric = "ROC")
  
  # Finally, run random forest using the optimal settings
  set.seed(421) # Reproducible randomness
  final_model <- randomForest(X_train_fold, y_train_fold, 
                              mtry = rf_model$bestTune$mtry,
                              splitrule = rf_model$bestTune$splitrule, 
                              min.node.size = rf_model$bestTune$min.node.size,
                              importance = TRUE)
  
  # Calculate and save model statistics (TRAINING DATA)
  train_pred_proba <- predict(final_model, X_train_fold, type = "prob")[, 2]
  train_auc <- auc(roc(y_train_fold, train_pred_proba))
  train_auc_scores <- c(train_auc_scores, train_auc)
  # Save predictions
  temp = tibble(row = c(1:nrow(X))[-fold], # Row from the original training dataset
                true_labels = y_train_fold, # Actual outcomes
                predicted_probabilities = train_pred_proba) # Predicted outcomes
  all_labels_train[[length(all_labels_train) + 1]] = temp
  
  # Calculate and save model statistics (TESTING DATA)
  test_pred_proba <- predict(final_model, X_test_fold, type = "prob")[, 2]
  test_auc <- auc(roc(y_test_fold, test_pred_proba))
  test_auc_scores <- c(test_auc_scores, test_auc)
  # Save predictions
  temp = tibble(row = c(1:nrow(X))[fold], # Row from the original testing dataset
                true_labels = y_test_fold, # Actual outcomes
                predicted_probabilities = test_pred_proba) # Predicted outcomes
  all_labels_test[[length(all_labels_test) + 1]] = temp
  
  # Save feature importance values
  # Ctrl, PD: higher number = higher in that group
  # MeanDecrease: how does the model quality decrease if the variable is removed? Higher=more important
  feature_importance_values[[length(feature_importance_values) + 1]] <- final_model$importance
}
```

```{r CALCULATE AVERAGE MODEL FROM K-FOLDS}

# Combine the AUC values using a bootstrap approach
boot_auc_train <- boot(data = train_auc_scores, statistic = function(d, i) mean(d[i]), R = 1000)
boot_auc_test <- boot(data = test_auc_scores, statistic = function(d, i) mean(d[i]), R = 1000)

# Calculate AUC confidence intervals
ci_train <- boot.ci(boot_auc_train, type = "perc")$percent[4:5] 
ci_test <- boot.ci(boot_auc_test, type = "perc")$percent[4:5]

# ci_train: if all AUC values are equal to 1 in the training dataset, the above fails. Interval is c(1,1).
# We'll assign this manually instead:
if(is.null(ci_train)) ci_train = c(1,1)

# Combine the test data labels and predictions (each fold is currently in a separate df)
# Each sample is only included in a test dataset once, so no need to average.
test_labels = bind_rows(all_labels_test)

# Combine the train data labels and predictions (each fold is currently in a separate df)
# Each sample is included 9 times, so need to average.
train_labels = bind_rows(all_labels_train) %>% 
  group_by(row,true_labels) %>% 
  summarize(predicted_probabilities = mean(predicted_probabilities)) %>% 
  ungroup()

# Calculate the average importance values for each variable.
# We will use Reduce to add together all values that are at the same row/column coordinate across all datasets in feature_importance_values, then divide by the number of datasets to get the mean.
mean_feature_importance = Reduce("+", feature_importance_values) / length(feature_importance_values)

# Convert to a data frame and add a Feature column
importance_df = data.frame(Feature = rownames(mean_feature_importance), mean_feature_importance)
# Remove row names
rownames(importance_df) = NULL
# Sort the results from most to least important. We will use MeanDecreaseGini.
importance_df = importance_df %>% arrange(-MeanDecreaseGini)
```

```{r PLOT RESULTS}
# Calculate ROC curve points from RF output
roc_test <- roc(test_labels$true_labels, test_labels$predicted_probabilities)
roc_train <- roc(train_labels$true_labels, train_labels$predicted_probabilities)

# True positive rate = sensitivity
# False positive rate = 1-specificity
ggplot() +
  # Training data: this is a type of control
  geom_line(aes(x = 1 - roc_train$specificities, y = roc_train$sensitivities), 
            color = "red",size=1) +
  # Test data: tells us the strength of the prediction
  geom_line(aes(x = 1 - roc_test$specificities, y = roc_test$sensitivities), 
            color = "black",size=1) +
  geom_abline(slope = 1, intercept = 0, color = "gray", linetype = "dashed",size=1) +
  labs(x = "False Positive Rate", y = "True Positive Rate") +
  annotate("text", x = 0.7, y = 0.2, label = sprintf("Train (red): %.2f (%.2f-%.2f)\nTest (black): %.2f (%.2f-%.2f)",
                                                     auc(roc_train), ci_train[1], ci_train[2],
                                                     auc(roc_test), ci_test[1], ci_test[2]), size = 6) +
  theme_minimal(base_size=18)

# Save plot results so we can easily compare different models
roc_data = data.frame(Dataset = 'RF Tutorial Data',
                      Training_AUC = round(mean(boot_auc_train$t), 2),
                      Training_AUC_CI = paste0(round(ci_train[1], 2), "-", round(ci_train[2], 2)),
                      Testing_AUC = round(mean(boot_auc_test$t), 2),
                      Testing_AUC_CI = paste0(round(ci_test[1], 2), "-", round(ci_test[2], 2)))
```

# Part 3: Assess Results

In the above plot, we can see that although the model is reasonably good at
detecting PD cases from controls in the test dataset (AUC=0.76), the training
dataset has an AUC of 1. This is not realistic and means that the model is
overfitted. The hyperparameter tuning should be adjusted in order to reduce this
effect. Ideally, the training dataset should curve similarly to the test dataset
(though it should still have a higher AUC).

The feature importance data can be portrayed in a bar plot or even as a table.
Below, we can see that Roseburia is the most important feature, which matches
the literature on Parkinson's disease.

```{r}
importance_df
```

```{r}
importance_df %>% 
  # Data are arranged by decreasing importance - turn it into a factor.
  # Otherwise the features will show up alphabetically instead.
  mutate(Feature = factor(.$Feature,levels = .$Feature)) %>% 
  ggplot(aes(Feature,MeanDecreaseGini,fill=MeanDecreaseGini)) +
  geom_col() +
  theme_classic(base_size=18) +
  theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1)) +
  ylab('Importance (Gini)') + xlab(NULL)
  
```