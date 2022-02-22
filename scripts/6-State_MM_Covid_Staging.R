# Load required packages ----
library(msm)
library(tidyverse)

# Import data ----
dm_stg_clean <- read_csv("/Users/sreynolds2/Documents/GitHub/MS-Markov_Covid_Staging/data/dm_stg_clean_11.08.21.csv")

# Frequency table for stage ----
table(dm_stg_clean$stage)

# Convert stages to stage categories: Moderate (M), Severe (S), and Dead ----
df <- dm_stg_clean %>%
    mutate(stagecat = case_when(stage %in% 4:5 ~ "M",
                                stage %in% 6:9 ~ "S",
                                stage == 10 ~ "Dead"),
           hx_sev = ifelse(stagecat=="S", 1, 0)) %>% 
    group_by(ID) %>% 
    mutate(hx_sev_yn = ifelse(sum(hx_sev>=1), 1, 0)) %>% 
    select(ID, date, stagecat, hx_sev_yn)

# Add new variable to indicate whether pt increased to moderate (M+), decreased to moderate (M-), remained at moderate (M) ----
df2 <- df %>%
    mutate(before = lag(stagecat),
           new_stage = case_when(stagecat == "M" & (is.na(before) | before=="M") ~ "M",
                             stagecat == "M" & before=="S" ~ "MS",
                             stagecat== "S" ~ "S")) %>% 
    mutate(new_stage_2 = if_else(stagecat=="M" & lag(new_stage)=="MS", new_stage))
    #select(ID, date, stage)

df3 <- df %>% 
    mutate(before = lag(stagecat),
           new_stage = case_when(stagecat=="M" & is.na(before) ~ "M",
                                 stagecat=="M" & hx_sev_yn==1 ~ "MS",
                                 stagecat=="M" & hx_sev_yn==0 ~ "M",
                                 stagecat=="S" ~ "S",
                                 TRUE ~ "M"))

#TRUE ~ stagecat
# Bind new row to represent discharge, ONLY if patient does not die ----
# If pt dead, then keep as-is and do not add row for discharge 

# Identify patients who did not die - and filter for max date 
alive <- df %>% 
    group_by(ID) %>% 
    filter(date==max(date),
           stage!="Dead") 

# For loop to create new row and bind row to `alive` 
for (i in length(alive)) {
    new_row <- tibble(ID=alive$ID, date=alive$date+1, stage="Discharge")
    new_df <- rbind(alive, new_row) %>% arrange(ID, date) %>% filter(stage=="Discharge")
}

# Bind newly created discharge rows to `df` and arrange by ID and date; assign to `markov`
markov <- rbind(df, new_df) %>% arrange(ID, date)

# Recode stages to numeric/integer type for Markov model (msm pkg requires numeric vars)
# M-: 1
# M : 2
# M+: 3
# S : 4
# Dead : 5
# Discharge : 6

markov <- markov %>% 
    mutate(state = case_when(stage=="M-" ~ 1,
                             stage=="M" ~ 2,
                             stage=="M+" ~ 3,
                             stage=="S" ~ 4,
                             stage=="Dead" ~ 5,
                             stage=="Discharge" ~ 6,
                             TRUE ~ 0)) %>% 
    select(-stage)

# Generate table of state transitions
statetable.msm(state = state, subject = ID, data = markov)

# Create matrix of "allowed" transitions
# where 1 = allowed transitions and 0 = forbidden transitions 
qmat <- matrix(c(0,1,0,1,1,1,
                 0,1,0,1,1,1,
                 0,1,0,1,1,1,
                 1,0,0,1,1,1,
                 0,0,0,0,0,0,
                 0,0,0,0,0,0),
               nrow = 6, ncol = 6,
               byrow = TRUE,
               dimnames(list(from=1:6, to=1:6)))

# Fit continuous time Markov multi-state model (MSM) by maximum likelihood
m_w_6_states <- msm(state ~ date, ID, data = markov, qmatrix = qmat, obstype = 2)

# Print MSM
printold.msm(m_w_6_states)

# Save as R object
saveRDS(m_w_6_states, "/Users/sreynolds2/Documents/GitHub/MS-Markov_Covid_Staging/6-State_Markov_M+/-")

# Probability matrix  
pmatrix.msm(m_w_6_states)

# Estim transition intensity (rate) matrix
qmatrix.msm(m_w_6_states)

# Log likelihood ratio
logLik.msm(m_w_6_states)

# Load 4-state_MM and 5_state_MM
m_w_4_states <- readRDS("/Users/sreynolds2/Documents/GitHub/MS-Markov_Covid_Staging/4-state_MM")
m_w_5_states <- readRDS("/Users/sreynolds2/Documents/GitHub/MS-Markov_Covid_Staging/5-state_MM")

# Calculate difference in log likelihood ratios between
# `m_w_4_states` and `m_w_6_states` -- higher log likelihood is usually considered the better model
diff_logLik <- logLik.msm(m_w_5_states) - logLik.msm(m_w_4_states)
diff_logLik

# Calculate AIC - lower AIC is usually considered the better model
AIC(m_w_6_states)
AIC(m_w_4_states)
AIC(m_w_5_states)