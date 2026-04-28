# Beluga Whales (Delphinapterus leucas) Discriminate Individual Callers Across Multiple Call Types

# Files Included: 
Beluga_Dataset_et_al_2026.xlsx

Beluga_Code.R

# Summary: 
The repository contains the analysis R code and dataset required to reproduce the statistical analyses, figures, and table for a study examining whether beluga whales discriminate individual callers across multiple call types. The study uses habituation-dishabituation playback experiments and analyzes behavioral responses (primarily number of looks toward the playback speaker) using generalized linear mixed models (GLMMs). The goal is to test whether individual identity is encoded and perceived across distinct call structures.


# File: Beluga_Dataset_et_al_2026.xlsx
•	Dataset used for all analyses 

•	Contains one row per session within a playback trial

# Column Names and Descriptions:
year: year of data collection

tod: time of day 

sod: playback sessions/ trial of day (first (1), second (2), third (3), fourth (4))

trial_number: unique trial identifier

session_order: playback session number (0- control tone, 1- first playback of call type with caller, 2- second playback of same call type with same caller, 3- third playback of same call type with same caller, 4- fourth playback of same call type with different caller)

playback_type: type of call played back during sessions

caller_ID: the whale who produced the call being played back

caller_sex: sex of the whale who produced the call being played back

caller_facility: facility of the whale who produced the call being played back

caller_target_relation: relation (kin, half kin, nonkin) between caller and target

target: the whale who is receiving the playback calls

target_sex: sex of the whale who is receiving the calls

num_app: total number of approaches per session

dur_app: total duration of approaches per session

avg_app: average length of approaches per session

num_look: total number of looks per session

dur_look: total duration of look per session

avg_look: average length of looks per session

# File: Beluga_Code.R

# Required Software: 
R Version 4.5.3 (2026-03-11)
RStudio Version 2025.09.2+418

# Required R Packages:
readxl 

dplyr 

tidyr 

ggplot2 

lme4 

glmmTMB 

emmeans 

performance 

effsize 

writexl

# Analyses were conducted using:
Platform: aarch64-apple-darwin20

Running under: macOS Sequoia 15.7.4

# Code Version and Reproducibility
Script version: v1.0

Manual versioning is used for this submission.

# Overview of Workflow
1.	Load dataset from Excel file. 
2.	Process the data by converting variables to factors and setting the reference level for session_order. 
3.	Subset dataset by call type (CCC, SCC, BiPT, SegW, BiRWS, FW, TBC). 
4.	Fit a Poisson generalized linear mixed model (GLMM) for each call type. 
5.	Assess dispersion using performance::check_overdispersion(). 
6.	If necessary, fit a negative binomial model and compare model fit using AIC. 
7.	Select the best-fitting model based on AIC and dispersion diagnostics. 
8.	Perform pairwise contrasts using emmeans with Tukey adjustment. 
9.	Calculate effect sizes (Cohen’s d) and 95% confidence intervals for key contrasts. 
10. Generate boxplot figures for number of looks by session order. 
11. Export pairwise comparison and effect size results to excel.

# Instructions to Run the Code
1.	Install required packages in R:
install.packages(c("readxl","dplyr","tidyr","ggplot2","lme4","glmmTMB","emmeans","performance","effsize","writexl")) 

2.	Set the working directory to the location of the repository:
setwd("path_to_repository") 

3.	Ensure the dataset file path is correct in the script:
beluga <- read_excel("Beluga_Dataset_et_al_2026.xlsx") 

4.	Run the analysis script:
source("Beluga_Code.R") 

5.	Outputs will include model summaries, pairwise comparisons, effect size tables, and figures displayed in the plotting window. 

# Additional Information:
•	No external configuration files are required.

•	All analysis parameters are defined within the R script.

•	Model structure is consistent across call types.

•	Distribution choice (Poisson or negative binomial) is determined based on dispersion diagnostics and AIC.

# Anonymization Statement
•	This README is anonymized for peer review.

•	Upon acceptance, the following will be added:
   
    - Author names and affiliations
    
    - Corresponding author/s contact details
    
    - Repository DOI (GitHub)
    
    - Link to preprint (if applicable)
    
    - Link to published article

# License
•	An open-source license (MIT) has been selected
