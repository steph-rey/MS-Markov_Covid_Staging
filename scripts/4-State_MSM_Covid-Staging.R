# @Title: 4-State Markov Model for Covid Staging Data
# @Project: Mindscape
# @Author: Stephanie.Reynolds@ucsf.edu
# @DateCreated: 01-25-2022
# @DateModified: 01-26-2022

# NOTE: This is a 4-state Markov Model where Moderate=1, Severe=2, Dead=3, and Discharged=4

# Load packages ----
library(tidyverse)
library(msm)

# Import data ----
raw <- read_csv("/Users/sreynolds2/Documents/GitHub/MS-Markov_Covid_Staging/data/dm_stg_clean_11.08.21.csv")

# Transform stages to stage categories ----
    # Stage 4-5 = 1 (Moderate)
    # Stage 6-9 = 2 (Severe)
    # Stage 10 = 3 (Dead)
    # Discharged = 4 
df <- raw %>%
    mutate(stagecat = case_when(stage %in% c(4, 5) ~ "1",
                                stage %in% 6:9 ~ "2",
                                stage == 10 ~ "3")) %>% 
    group_by(ID) %>% 
    mutate(adm_stagecat = stagecat[date==min(date)]) %>% 
    select(ID, date, stagecat) 

# Prep dataset for Markov Model ----
    # For each 'ID', create new row to represent discharge state
    # Bind to dataset and arrange by 'ID' and 'date'
    # Set stagecat to '4' or 'discharged'
    markov <- df %>%
    group_by(ID) %>%
    summarize(date = max(date) + 1) %>%
    mutate(stagecat = NA) %>%
    rbind(df) %>%
    arrange(ID, date) %>%
    mutate(stagecat = case_when(is.na(stagecat) ~ "4",
                                TRUE ~ stagecat))

# Convert stagecat to ordinal factor variable ----
markov$stagecat <- factor(markov$stagecat, 
                          ordered = TRUE, 
                          levels = c("1", "2", "3", "4"))

# Generate table of stage transitions ----
    # Calculate freq table counting # of times each pair of states were observed in succession
    # rows = starting states, and columns = ending states
statetable.msm(state = stagecat, subject = ID, data = markov)

# Create matrix of "allowed" transitions ---- 
    # 1 = allowed transitions
    # 0 = forbidden transitions 
qmat <- matrix(c(1,1,1,1,
                  1,1,1,1,
                  0,0,0,1,
                  0,0,0,1),
                nrow = 4, ncol = 4,
                byrow=TRUE,
                dimnames=list(from=1:4, to=1:4))

# Fit continuous time Markov multi-state model (MSM) by maximum likelihood ----
    # Set obstype = 2 --> Considered an "exact transition time." Even though this 
    # is not necessarily true with our dataset since it's more discrete than continuous, 
    # we are assuming that the pt remains at the state at previous observation until 
    # the current observation. It does not allow for a change in state in between 
    # observation points.
m <- msm(stagecat ~ date, ID, data = markov, qmatrix = qmat, obstype = 2)

# Print MSM ----
printold.msm(m)

# Total LOS for patients at moderate (1) and severe (2) stage categories ----
totlos.msm(m, end = c("1", "2"))

# Probability matrix ----
pmatrix.msm(m)

# Estimated transition intensity (rate) matrix ----
qmatrix.msm(m)

# Log likelihood ratio ----
logLik.msm(m)

