# Code to categorize patients into Moderate (M), Moderate with Hx of Severe (S), Severe (S), and Dead 
# Last Modified: 2022-04-11

# Load required packages ----
library(slider)
library(lubridate)
library(tidyverse)

# Import data ----
dm_stg_clean <- read_csv("/Users/sreynolds2/Documents/GitHub/MS-Markov_Covid_Staging/data/dm_stg_clean_11.08.21.csv")

# Categorize patients into Moderate (M), Moderate with Hx of Severe (S), Severe (S), and Dead ---- 
    # Moderate w/o history of severe (M)
    # Moderate w/ history of severe (MS)
    # Severe (S)
    # Dead (Dead)
df <- dm_stg_clean %>%
    mutate(stage_cat = case_when(stage %in% 4:5 ~ "M",
                                 stage %in% 6:9 ~ "S",
                                 stage == 10 ~ "Dead"),
           severe_yn = ifelse(stage_cat=="S", 1, 0)) %>% 
    group_by(ID) %>% 
    mutate(hx_severe = slider::slide_index_sum(x = severe_yn, 
                                               i = date, 
                                               before = lubridate::days(365)),
           state = case_when(stage_cat=="Dead" ~ "Dead",
                             stage_cat=="S" ~ "S",
                             stage_cat=="M" & hx_severe==0 ~ "M",
                             stage_cat=="M" & hx_severe>0 ~ "MS")) %>% 
    select(ID, date, state, age, sex, BMI, race, ethnicity)
    # Include covariates, e.g., age, sex, BMI, race, ethnicity 

################################################################################

# Code to build and fit Markov Model using msm package 
# Last Modified: 2022-04-11

# Load msm package ----
library(msm)

# Stage transition table ----
statetbl <- statetable.msm(state, ID, data=df)

# Matrix of allowed transitions, where 0=not allowed and 1=allowed ----
qmat <- rbind(c(0,0,0,0),
              c(1,1,0,1),
              c(1,0,1,1),
              c(1,0,1,1))

# Assign row and column names to qmat ----
rownames(qmat) <- c("Dead",rownames(statetbl))
colnames(qmat) <- colnames(statetbl)

# Reformat `state` variable and create variable `DaysSinceEntry` ----
df$state <- as.numeric(as.factor(df$state))

#df$JulianDay <- as.numeric(df$date-min(df$date))

df <- df %>% 
    group_by(ID) %>% 
    mutate(DaysSinceEntry = as.numeric(date-min(date)))

# Specify initial values and assign to `Q.crude` ----
Q.crude <- crudeinits.msm(state ~ DaysSinceEntry, ID, data=df, qmatrix=qmat)

# BASELINE MODEL - Run msm to fit model ----
m0 <- msm(state ~ DaysSinceEntry, subject=ID, data=df, qmatrix=Q.crude)

# COVARIATE MODELS ----
# Age ----
df$age.cen <- df$age/100

m1 <- msm(state ~ DaysSinceEntry, subject=ID, data=df, qmatrix=Q.crude,
          covariates = ~ age.cen)

m1a <- msm(state ~ DaysSinceEntry, subject=ID, data=df, qmatrix=Q.crude,
           covariates = list("2-4" = ~ age.cen,
                             "3-4" = ~ age.cen,
                             "4-3" = ~ age.cen))

# Evaluate log likelihoods and compare 
pvalue.age <- 1 - pchisq((logLik(m1a) - logLik(m0)) * 2, df=3)
    # Tells us that age is important; sig P value
    # Will need data on vacc status - age liked confounded by vacc status 

# Evaluate AIC and compare
AIC(m0, m1, m1a)

# BMI ----
df$BMI.cen <- df$BMI/100

m2 <- msm(state ~ DaysSinceEntry, subject=ID, data=df, qmatrix=Q.crude,
          covariates = ~ BMI.cen)

m2a <- msm(state ~ DaysSinceEntry, subject=ID, data = df, qmatrix = Q.crude,
           covariates = list("2-4" = ~ BMI.cen,
                             "3-4" = ~ BMI.cen,
                             "4-3" = ~ BMI.cen))

# Evaluate AIC and compare
AIC(m0, m2, m2a)

# Sex ---- 
df <- df %>% 
    mutate(Male = ifelse(sex=="Male", 1, 0),
           Female = ifelse(sex=="Female", 1, 0))

#m3 <- msm(state ~ DaysSinceEntry, subject=ID, data=df, qmatrix=Q.crude, covariates = ~ Male + Female)

m3 <- msm(state ~ DaysSinceEntry, subject=ID, data=df, qmatrix=Q.crude,
          covariates = ~ as.factor(sex))

m3 <- msm(state ~ DaysSinceEntry, subject=ID, data=df, qmatrix=Q.crude,
          covariates = ~ Male + Female)

# Race ----
df <- df %>% 
    mutate(Black = ifelse(race=="Black or African American", 1, 0),
           Asian = ifelse(race=="Asian", 1, 0),
           White = ifelse(race=="White or Caucasian", 1, 0),
           Race_Other = ifelse(race %in% c("Native Hawaiian or Other Pacific Islander", "Other", "Unknown"), 1, 0))

#m4 <- msm(state ~ DaysSinceEntry, subject=ID, data=df, qmatrix=Q.crude, covariates = ~ Black + Asian + White + Race_Other)

m4 <- msm(state ~ DaysSinceEntry, subject=ID, data=df, qmatrix=Q.crude,
          covariates = ~ as.factor(race))

# Ethnicity ---- 
df <- df %>% 
    mutate(Hispanic = ifelse(ethnicity=="Hispanic or Latino", 1, 0),
           NotHispanic = ifelse(ethnicity=="Not Hispanic or Latino", 1, 0),
           Ethnicity_Unk = ifelse(ethnicity=="Unknown or Caucasian", 1, 0))

#m5 <- msm(state ~ DaysSinceEntry, subject=ID, data=df, qmatrix=Q.crude, covariates = ~ Hispanic + NotHispanic + Ethnicity_Unk)

m5 <- msm(state ~ DaysSinceEntry, subject=ID, data=df, qmatrix=Q.crude,
          covariates = ~ as.factor(ethnicity))

# Include all covariates: Race, ethnicity, sex, age, BMI ----
m.full <- msm(state ~ DaysSinceEntry, subject=ID, data=df, qmatrix=Q.crude,
          covariates = ~ as.factor(ethnicity) + as.factor(race) + as.factor(sex) + age.cen + BMI.cen)
m.full

# Compare full and baseline models ---- 
# AIC
AIC(m0, m.full)

# Log likelihood
logLik(m.full) - logLik(m0)

# Summary of obs/pred prevalances at each state and hazard ratios for covariate effects ----
summary.msm(m.full)

# Plot of observed and expected prevalences ----
plot.prevalence.msm(m.full)

# Transition probability matrix ----
pmatrix.msm(m.full)

################################################################################

# OTHER FUNCTIONS FROM MSM PACKAGE ----

# Rough estimate of goodness of fit of model ----
    # Estimates observed number of individuals at each state at a series of times
    # Compares these with predicted numbers 
# prevalence.msm()

# Plot of observed and expected prevalences ----
# plot.prevalence.msm()

# Pearson-type goodness-of-fit test for model ---- 
# pearson.msm()

# Summarize observed and expected prevalences at each state for each time ----
# Also prints hazard ratios for models with covariate effects 
# summary.msm(m1a)

# List of tables containing odds ratio estimates, one table per covariate -----
# odds.msm(m1a)

################################################################################

# OLD CODE - DISREGARD 

# Print MSM
printold.msm(m.full)

# Save as R object
saveRDS(m.full, "/Users/sreynolds2/Documents/GitHub/MS-Markov_Covid_Staging/MSM_full.rds")

# Probability matrix  
pmatrix.msm(m.full)

# Estim transition intensity (rate) matrix
qmatrix.msm(m.full)
