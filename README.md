##### MICB 475 Team 8
# March 10 2025 Team Meeting
### Agenda:
- discuss indicator species analysis results
- next steps: knn and random forest; where to start? weekly timeline?
- delegate team tasks for the week

### Meeting Minutes:
- Attendees: Emily W, Emily H, Sebastian
- Eddie absent; no indicator taxa results at the moment
- distribution of roles/tasks
  - Sebastian comfortable with coding
  - everyone good with distribution so far
- after this week: core microbiome, machine learning
- ML will probably mostly be Eddie
- by next week: alpha/beta, indicator taxa, core microbiome; have preliminary results
- last week, mostly worked on rewriting proposal
- Sebastian rerunning alpha/beta with newly subsetted monoinfected depression vs no depression
- Shannon for depression: no change
- want to look at Shannon, Faith, weighted/unweighted unifrac
- saving images: save as pdf, easier to edit later
- 2 more meetings until oral presentation!
- Emily W ran Eddie's indicator taxa code from GitHub
  - INSTI group sees greater richness in Clostridia
  - INSTI group depletion of certain taxa
  - no INSTI has more indicator species
  - look into these species, if there's prior research on the taxa wrt HIV
  - Actinobacteria is INSTI-specific
  - Sebastian ran the code and sees gamma in INSTI but Emily W doesn't; need to set seed
- Sebastian to finish alpha/beta diversity for HIV w/ depression
- next steps:
  - look into indicator species
  - set a seed for indicator analysis
  - core microbiome
  - differential analysis
  - upload figures to GitHub for next week

# March 4 2025 Team Meeting
### Agenda:
- Discuss team proposal and feedback
- Discuss alpha and beta diversity results
- Plan for this week

### Meeting Minutes:

Introduction: what do we know and what are we missing
- mention machine learning approach because it’s an advancement in the field
- move dataset overview to the last paragraph of the introduction
  - "in order to look at this question, the authors used this data…"
- add information about the insti drug classes to introduction
  - hiv —> hiv arvs —> hiv arvs with relation to microbiome —> details of our specific project

Experimental Aims: split first aim into two aims and reshape the others
1. hiv and depression — see if there's an effect
2. alpha and beta diversity for INSTI treatment
3. indicator species
4. machine learning
in between aim 2 and 3 is functional analysis, which will not be added to the proposed approach but which we will do if we have time

Title: revamp to take out depression

Proposed Approach:
- take out depression from aim 2 and add depression comparisons to aim 1
- remove depression comparisons from aim 2

Metadata:
- we want to add sex, age, and bmi to the metadata
- hopefully this will help with our PCoA analysis, which is currently lacking

Current Updates:
- we see a difference in beta diversity between INSTI+ and INSTI- PLWH but no difference for alpha diversity
- this is unusual in the literature but matches up with what these researchers saw with their dataset



# February 25 2025 Team Meeting
### Agenda:
- Review proposal and plan of action
- Discuss our machine learning model and decide on a final model approach
- Discuss whether we should look at depression for all aims

### Meeting Minutes

  Machine Learning Script } Python or R version

  Proposal Overview
  - good experimental aims with clear layout
  - metadata wrangling
    - 4 categories : INSTI or no, depression or no
    - did rarefaction on Qiime2 as R was taking too long
    - Github includes both the whole metadata and the one that only includes the columns we are using
   
    Next Steps
    - show that HIV positive patients have lower alpha diversity than patients without HIV (HIV+ vs HIV-HCV-)
    -  HIV+dep+ vs HIV+dep-
    -   ake a histogram that shows how many people are taking each drug type (or can do pie chart if there are too may columns) and another one that shows how many drugs people take (ex. 10 people take 2 drugs, 20 people take 3 drugs)
         -- general/broad overview of the entire dataset
         -- find a drug as a control that has a decent amount of people taking and not taking and use this as our comparison to INSTI
    - first aim: just look at depression and drug as a preliminary correlation, may be hard to study because we don't know the order of events (are the patients depressed before they started to take the drug)
         -- show this in a pie chart - abandon correlation
         -- if theres no difference than disregard depression here on out
    - generate the alpha and beta diversity matrices with HIV+INSTI+ vs HIV+ INSTI- } shannon, faith and unweighted unifrac
    - statistical test and PCoA plot
    - compare HIV+INSTI+ vs HIV-HCV- to see if alpha diversity is restored or just increased
   

# February 11 2025 Team Meeting

### Agenda

- Discuss column selection and organization for metadata file completion.
- Discuss control for our samples: do we need to remove any samples (such as the BLANKS) ? 
- Discuss analysis requirements: are there specific metadata fields we should consider for downstream visualizations?
- Discuss project proposal objectives.

### Meeting Minutes

How we rebuilt the metadata:
Only a subset of the project was 16S, and the accession numbers in the metdata were not complete (only 385 of 1031 rows have an accession number). By combining multiple data files available for this dataset, we were able to achieve a metadata with 1030 samples, looking at ~1800 variables. Data has been processed and denoised by Sebastian.

- all fastq files on the server are accounted for

We'll leave out the BLANK items, which leaves us with ~800 samples

INSTI: promotes higher alpha diversity
All other drugs (NRTI, NNTRI, PI, etc.) decreased alpha diversity
--> might be worth looking into the physiological mechanisms of the medications, as that may be how it's affecting alpha diversity.
--> we can divide up the drug regimens based on whether they take INSTI or not

variables we're looking at:
- HIV-only (667 patients) & HIV-HCV coinfection (154 patients)
- INSTI & no-INSTI
- depression (potentially, if we don't see anything after looking at the above variables) (we may want to look at if depression overrides the impact of drugs on the microbiome)
    --> we would split based on the Beck depression cutoff (either by the numerical values, or categorical (eg. mild, moderate)

Git commands
- git pull origin main
- git add .
- git commit -m ""
- git push origin main

Proposal discussion
- don't force a research gap
- experimental aims
      - QIIME
      - diversity
      - composition
      - indicator taxa + machine learning (likely a clustering model)
- cite all packages we use, cite QIIME2

ACTION --> ask Evelyn for the machine learning script

ACTION --> Bessie will let us know if we can format our references in a style other than ASM

# February 4 2025 Team Meeting

### Meeting Minutes

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

### Meeting Minutes 

No official meeting notetaker. Preliminary discussion about data wrangling will be done in [this Google Doc](https://docs.google.com/document/d/19ViECbRmhkQHRDq6u6QZakAeyeBRO01sz017HTwXPtE/edit?usp=sharing).
