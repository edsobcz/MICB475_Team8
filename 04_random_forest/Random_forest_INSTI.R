library(randomForest)
library(caret)
library(pROC)
library(boot)
library(ggplot2)
library(phyloseq)
library(tidyverse)

# Load the example dataset.
# This phyloseq object can be created by running Random Forest - Create Dataset.Rmd.
ps = readRDS('processed_phyloseq.rds') %>% 
  tax_glom('Genus') # We will only use Genus-level data.

#----- Part 1 --------
# Calculate average abundance:
require(microbiome)
avg_abundance = ps %>%
  microbiome::transform('compositional') %>% # relative abundance
  psmelt() %>% 
  group_by(Genus) %>% # Perform the next line separately for each Genus value
  summarize(Abundance = mean(Abundance)) %>%
  ungroup() # Stop perfoming subsequent code separately for each Genus

keep_taxa = avg_abundance %>% 
  arrange(-Abundance) %>% # Sort high to low
  pull(Genus) %>% # Get Genus column as a vector
  .[c(1:10)] # First 10 values (. = 'everything in keep_taxa up until this line')

# Normalize microbiome data
df_clr = ps %>% microbiome::transform('clr') %>%  # CLR transformation
  subset_taxa(Genus %in% keep_taxa) %>% 
  psmelt() %>% 
  # Add Z transformation to each genus individually
  group_by(Genus) %>% 
  mutate(Abundance = scale(Abundance)) %>% 
  ungroup()

# Select only the metadata of interest
df_filt = df_clr %>% 
  select(Sample,INSTI_drug_current,hiv_status_clean,hcv,bdi_ii,ethnicity,host_age,host_height,host_weight,Genus,Abundance) %>% 
  # Turn each genus into its own column
  pivot_wider(names_from = Genus, values_from = Abundance)

# Remove rows with NA values in the metadata
df_filt = df_filt %>% na.omit()

df_filt <- df_filt %>%
  filter(!is.na(host_age)) %>%
  mutate(host_age = as.numeric(host_age))

df_filt <- df_filt %>%
  filter(!is.na(host_weight)) %>%
  mutate(host_weight = as.numeric(host_weight))

# Test if the continuous metadata need additional transformations to be normalized
hist(df_filt$host_age) # Normally distributed
hist(df_filt$host_weight) # Reasonably normal



# Our final dataset will have the sample ID, outcome variable of interest, possible predictive variables from the metadata, and then the microbiome data.

# Z-transform any numerical data
for_RF = df_filt %>% mutate_if(is.numeric,~as.numeric(scale(.)))
for_RF_nosampleID = for_RF %>% select(-Sample) # Not necessary for RF


#---- Part 2 ------
# Divide the table into predictors and outcome.
X <- for_RF_nosampleID %>% select(-INSTI_drug_current)
y <- for_RF_nosampleID %>% pull(INSTI_drug_current)

# Convert y to a factor. If continuous, the order should be c(control,case).
y = factor(y,levels = c('NO','YES'))

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

#---- Graph results -----
# Calculate ROC curve points from RF output
roc_test <- roc(test_labels$true_labels, test_labels$predicted_probabilities)
roc_train <- roc(train_labels$true_labels, train_labels$predicted_probabilities)

# True positive rate = sensitivity
# False positive rate = 1-specificity
ROC_curve<- ggplot() +
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


Importance_features <- importance_df %>% 
  # Data are arranged by decreasing importance - turn it into a factor.
  # Otherwise the features will show up alphabetically instead.
  arrange(desc(MeanDecreaseGini)) %>%
  mutate(Feature = factor(Feature, levels = Feature)) %>% 
  ggplot(aes(x = Feature, y = MeanDecreaseGini, fill = MeanDecreaseGini)) +
  geom_col() +
  scale_fill_gradient(low = "#9ECAE9", high = "navy", name = "Mean\nDecrease\nGini") +
  labs(y = "Importance (Gini)", x = NULL) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    panel.grid.major.x = element_blank())

Importance_features


# Create a data frame
auc_df <- data.frame(
  Type = c(rep("Train", length(train_auc_scores)), 
           rep("Test", length(test_auc_scores))),
  AUC = c(train_auc_scores, test_auc_scores)
)

# Convert lists to data frames
train_df <- do.call(rbind, all_labels_train)
test_df <- do.call(rbind, all_labels_test)

# Add a column to identify the source
train_df$Source <- "Train"
test_df$Source <- "Test"

# Combine the dataframes
combined_df <- rbind(train_df, test_df)


# Combine all train and test labels
train_labels_df <- do.call(rbind, all_labels_train)
test_labels_df <- do.call(rbind, all_labels_test)

comprehensive_comparison <- data.frame(
  Fold = 1:length(train_auc_scores),
  Train_AUC = train_auc_scores,
  Test_AUC = test_auc_scores,
  Train_Mean_Prob = sapply(all_labels_train, function(x) mean(x$predicted_probabilities)),
  Test_Mean_Prob = sapply(all_labels_test, function(x) mean(x$predicted_probabilities))
)

# Add summary statistics
summary_stats <- data.frame(
  Fold = c("Mean", "SD"),
  Train_AUC = c(mean(comprehensive_comparison$Train_AUC), 
                sd(comprehensive_comparison$Train_AUC)),
  Test_AUC = c(mean(comprehensive_comparison$Test_AUC), 
               sd(comprehensive_comparison$Test_AUC)),
  Train_Mean_Prob = c(mean(comprehensive_comparison$Train_Mean_Prob), 
                      sd(comprehensive_comparison$Train_Mean_Prob)),
  Test_Mean_Prob = c(mean(comprehensive_comparison$Test_Mean_Prob), 
                     sd(comprehensive_comparison$Test_Mean_Prob))
)

# Combine the original data with summary statistics
final_comparison <- rbind(comprehensive_comparison, summary_stats)

# View the table
print(final_comparison)

# Formatted table
library(knitr)
kable(final_comparison, caption = "Comprehensive Model Performance")

ggsave("ROC_curve_random_forest.png", plot = ROC_curve)
ggsave("Importance_features.png", plot = Importance_features)
