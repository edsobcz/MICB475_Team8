
library(tidyverse)
library(dplyr)
library(stringr)

##### Load initial files #####
Lmeta <- read.delim("Long_metadata.txt")
depression_mdata <- read.delim("depression_metadata.txt")
Full_exp_list <- read.delim("Experiment_list.txt")[,c("run_accession","experiment_accession")]

#Proces original metadata file
Ometa <- depression_mdata[!is.na(depression_mdata$sample.id), ]
Ometa$Submitter_Id <- gsub("^qiita_sid_11135:", "", Ometa$Submitter_Id)
Ometa <- rename(Ometa, sample_name = Submitter_Id, run_accession = sample.id)

#Compile list of samples that belong to 16S
S1_16S <- read.delim("16S1.tsv")[,c("sample_name","hcv")] #Note: used hcv to keep dataframe format for later processing. 
S2_16S <- read.delim("16S2.tsv")[,c("sample_name","hcv")]
S3_16S <- read.delim("16S3.tsv")[,c("sample_name","hcv")]
S4_16S <- read.delim("16S4.tsv")[,c("sample_name","hcv")]

compiled_list_16S <- bind_rows(S1_16S, S2_16S, S3_16S, S4_16S) %>%
  distinct(sample_name, .keep_all = TRUE)

#Compile list that map Experiment number to sample name that belong to 16S
E1_16S <- read.delim("E1.tsv")
E2_16S <- read.delim("E2.tsv")
E3_16S <- read.delim("E3.tsv")
E4_16S <- read.delim("E4.tsv")

compiled_experimetn_list_16S <- bind_rows(E1_16S, E2_16S, E3_16S, E4_16S) %>%
  distinct(sample_name, .keep_all = TRUE)

#Cross reference lists with matadata to create final metadata file
cross_exp_16S <- inner_join(compiled_experimetn_list_16S, compiled_list_16S, by = "sample_name")

cross_exp_16S_run <- inner_join(Full_exp_list, cross_exp_16S, by = "experiment_accession")

cross_exp_16S_run_meta <- inner_join(Ometa, cross_exp_16S_run, by = "run_accession")

cross_exp_16S_run_meta <- inner_join(Ometa, cross_exp_16S_run, 
                                     by = "run_accession", 
                                     suffix = c("_remove", ""))
cross_exp_16S_run_meta_Lmeta <- inner_join(Lmeta,cross_exp_16S_run_meta, by = "sample_name",suffix = c("_remove", ""))

check <- cross_exp_16S_run_meta_Lmeta[,c("sample_name", "run_accession")]

cross_exp_16S_run_meta_Lmeta_final <- cross_exp_16S_run_meta_Lmeta %>%
  select(-contains("_remove"))

#Save final csv
write_csv(cross_exp_16S_run_meta_Lmeta, "final_metadata.csv")

#### Creating metadata file for analysis ####

insti_drugs <- c("CAB", "DTG", "RAL", "BIC", "EVG")

sub_metadata <- cross_exp_16S_run_meta_Lmeta_final[,c("sample_name","run_accession","sample_type","regimen_type","combo_regimen","priorvisitreg","current_regimen","hiv_status_clean","hcv","bdi_ii")] %>%
  mutate(INSTI_drug_current = case_when(str_detect(current_regimen, paste(insti_drugs, collapse="|")) ~ "YES", TRUE ~ "NO")) %>%
  mutate(INSTI_drug_prior = case_when(str_detect(priorvisitreg, paste(insti_drugs, collapse="|")) ~ "YES", TRUE ~ "NO")) %>%
  select("run_accession","sample_name","sample_type","INSTI_drug_current","INSTI_drug_prior","current_regimen","hiv_status_clean","hcv","bdi_ii") %>%
  rename(experiment_name = sample_name,`sample-id` = run_accession) #Note: this renaming is done for qiime2 metadata requirments.

#Save modified version of metadata file
write_tsv(sub_metadata, "modified_metadata.tsv")

