##### MICB 475 Team 8

# February 11 2025 Team Meeting

Agenda

- Discuss column selection and organization for metadata file completion.
- Discuss control for our samples: do we need to remove any samples (such as the BLANKS) ? 
- Discuss analysis requirements: are there specific metadata fields we should consider for downstream visualizations?
- Discuss project proposal objectives. 

# February 4 2025 Team Meeting

Depression dataset preliminary wrangling - works well

Confirming which measurements we will be using
- HIV drugs: use drug categories already listed in the metadata, some booster drugs to increase effects of the antivirals they are on (can ignore these in classification as they are combined with a specific drug, or can just say they are using a booster)
- Prediction model predicting microbiome based on what drugs (random forest model machine learning to predict the top taxon)

Need to subset the dataset - 30% v 70%

Confounding variables
- coinfection with HIV
- depression
  --> may need to separate into cohorts to isolate the effects of the drugs
  --> minimal sample size (if ~50 samples would be a big limiation)

--> possible approach } cohort based on coinfection, then use depression status as a second variable

Does the class of HIV drug determine the patients microbiome? If so, we would need to class the microbiome.

- keep 3 metadata variables (depression, coinfection, drug type)
- look at clustering and figure out which one best describes the clustering
- determine which taxa are related to each level of that variable (ex. if depression is the driving variable, see which taxa is most linked to depressed status)
  --> Go through every analysis (alpha and beta) between these three variables (just to see which variable has the biggest effect through taxonomic analysis) then will move forward with machine learning to build a taxonomic prediction model
  --> need to re-label depression score into categorical (says there is a Beck Depression Inventory but we can't find the column, would be best if it was a yes/no but may need to be more nuanced categories)

--> original paper found co-infection has significant differences, and depression was only a factor in co-infection
- since depression did not have affect in single HIV infection, can just look at HIV drugs and not control for depression
  **need to contact authors to ask if the data is only co-infection or mono-HIV infection (or a mix)

Team Proposal } due after reading week, come to team meeting next week with questions
- work through outline of proposal in team meeting
- includes a dataset overview section (need to complete processing to do this)
- will get team server credentials (but keep a copy of qza and qzv files on Github)
- all scripts, qiime files etc. should be on Github
- processing on server, analysis on local R

# 28 January 2025 Team Meeting
No official meeting notetaker. Preliminary discussion about data wrangling will be done in [this Google Doc](https://docs.google.com/document/d/19ViECbRmhkQHRDq6u6QZakAeyeBRO01sz017HTwXPtE/edit?usp=sharing).