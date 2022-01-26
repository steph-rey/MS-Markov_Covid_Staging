# @Title: 5-State Markov Model for Covid Staging Data
# @Project: Mindscape
# @Author: Stephanie.Reynolds@ucsf.edu
# @DateCreated: 01-25-2022
# @DateModified: 01-25-2022

# NOTE: This is a 5-state Markov Model where LowMod=1, HighMod=2, Severe=3, Dead=4, and Discharged=5
# This model will be compared to the 4-state MM.

# Load packages ----
library(tidyverse)
library(msm)

# Import data ----
raw <- read_csv("/Users/sreynolds2/Documents/GitHub/MS-Markov_Covid_Staging/data/dm_stg_clean_11.08.21.csv")

# Frequency table for stage ----
table(raw$stage)

# Transform stages to stage categories ----
    # Stage 4 -> 1 (Low Moderate)
    # Stage 5 -> 2 (High Moderate)
    # Stage 6-9 -> 3 (Severe)
    # Stage 10 -> 4 (Dead)
    # Discharged -> 5 (Discharged)
df <- raw %>%
    mutate(stagecat = case_when(stage == 4 ~ "1",
                                stage == 5 ~ "2",
                                stage %in% 6:9 ~ "3",
                                stage == 10 ~ "4")) %>%
    select(ID, date, stagecat)

# Prep dataset for Markov Model ----
    # For each 'ID', create new row to represent discharge state
    # Bind to dataset and arrange by 'ID' and 'date'
    # Set stagecat to '4' or 'discharged'
markov_w_5_stages <- df %>%
    group_by(ID) %>%
    summarize(date = max(date) + 1) %>%
    mutate(stagecat = NA) %>%
    rbind(df) %>%
    arrange(ID, date) %>%
    mutate(stagecat = case_when(is.na(stagecat) ~ "5",
                                TRUE ~ stagecat))

# Convert stagecat to ordinal factor variable ----
markov_w_5_stages$stagecat <- factor(markov_w_5_stages$stagecat, 
                          ordered = TRUE, 
                          levels = c("1", "2", "3", "4", "5"))

# Generate table of state transitions ----
    # Calculate freq table counting # of times each pair of states were observed in succession
    # rows = starting states, and columns = ending states
statetable.msm(state = stagecat, subject = ID, data = markov_w_5_stages)

# Create matrix of "allowed" transitions ---- 
    # 1 = allowed transitions
    # 0 = forbidden transitions 
qmat1 <- matrix(c(1,1,1,1,1,
                  1,1,1,1,1,
                  1,1,1,1,1,
                  0,0,0,0,1,
                  0,0,0,0,1),
                nrow = 5, ncol = 5,
                byrow=TRUE,
                dimnames=list(from=1:5, to=1:5))

# Fit continuous time Markov multi-state model (MSM) by maximum likelihood ---- 
    # Set obstype = 2 --> Considered an "exact transition time." Even though this 
    # is not necessarily true with our dataset since it's more discrete than continuous, 
    # we are assuming that the pt remains at the state at previous observation until 
    # the current observation. It does not allow for a change in state in between 
    # observation points.
m_w_5_states <- msm(stagecat ~ date, ID, data = markov_w_5_stages, qmatrix = qmat1, obstype = 2)

# Print MSM ----
printold.msm(m_w_5_states)

# Total LOS for patients at low moderate (1), high moderate (2), and severe (3) stage categories ----
totlos.msm(m_w_5_states, end = c("1", "2", "3"))

# Probability matrix ----
pmatrix.msm(m_w_5_states)

# Estimated transition intensity (rate) matrix ----
qmatrix.msm(m_w_5_states)

# Log likelihood ratio ----
logLik(m_w_5_states)


# Comparing both models ----

# Calculate difference in log likelihood ratios between `m` and `m_w_5_states` ----
diff_logLik <- logLik.msm(m_w_5_states) - logLik.msm(m)

# Difference in Restricted AIC (DRAIC) ---- 
    # Compare predictive ability of two Markov multi-state models (MSM)
    # NOTE: The functions below take forever to run.

# DRAIC_results <- draic.msm(msm.full = m_w_5_states, msm.coarse = m)
# draic_value <- DRAIC_results$draic

# Difference in Restricted AIC (DRAIC) = 0.144 indicating that the coarse model (one w/ fewer parameters) is preferred.
    # Positive value for DRAIC = coarse model is preferred, while negative value = full model is preferred.

# Calculate standard AIC for both models ----
AIC(m, m_w_5_states)






# DISREGARD CODE BELOW ---- 
# Use code below if wanting to add variable for stagecat at ADMISSION

# df <- raw %>%
mutate(stagecat = case_when(stage == 4 ~ "1",
                            stage == 5 ~ "2",
                            stage %in% 6:9 ~ "3",
                            stage == 10 ~ "4")) %>% 
    group_by(ID) %>% 
    mutate(adm_stagecat = stagecat[date==min(date)]) %>% 
    select(ID, date, stagecat, adm_stagecat) 
# table(df$adm_stagecat)
