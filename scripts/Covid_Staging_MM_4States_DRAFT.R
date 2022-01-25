# @Title: Prep Covid Staging Data for Markov Model
# @Project: Mindscape
# @DateCreated: 01-25-2022
# @DateModified: 01-25-2022

# @ModelStates: Moderate, Severe, Dead, Discharged

# Load packages
library(tidyverse)
library(msm)

# Import data
raw <- read_csv("/Users/sreynolds2/Documents/GitHub/MS-Markov_Covid_Staging/data/dm_stg_clean_11.08.21.csv")

# Transform stages to stage categories where:
# stages 4-5 = 1 (moderate)
# stages 6-9 = 2 (severe)
# stage 10 = 3 (dead)
# discharged = 4 ???
df <- raw %>%
    mutate(stagecat = case_when(stage %in% c(4, 5) ~ "1",
                                stage %in% 6:9 ~ "2",
                                stage == 10 ~ "3")) %>% 
    group_by(ID) %>% 
    #mutate(adm_stagecat = stagecat[date==min(date)]) %>% 
    select(ID, date, stagecat) #, adm_stagecat) 

#table(df$adm_stagecat)

# Prep dataset for Markov Model
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

# Convert stagecat to ordinal factor variable 
markov$stagecat <- factor(markov$stagecat, 
                          ordered = TRUE, 
                          levels = c("1", "2", "3", "4"))
# Run statetable.msm fxn
statetable.msm(state = stagecat, subject = ID, data = markov)

# Create matrix to show allowed transitions where:
# 1 = allowed transitions
# 0 = forbidden transitions
qmat0 <- matrix(c(1,1,1,1,
                  1,1,1,1,
                  0,0,0,1,
                  0,0,0,1),
                nrow = 4, ncol = 4,
                byrow=TRUE,
                dimnames=list(from=1:4, to=1:4))
# Print qmatrix
qmat0

# Determine crude estimate of transition rates
crudeinits.msm(formula = stagecat ~ date, subject = ID, data = markov, qmatrix = qmat0)

# Run msm fxn
# Set obstype = 2
m <- msm(stagecat ~ date, ID, data = markov, qmatrix = qmat0, exacttimes = T)

totlos.msm(m, end = c("1", "2"))
pmatrix.msm(m)
logLik.msm(m)

df %>% group_by(ID) %>% 
    mutate(adm_stagecat = stagecat[date==min(date)]) %>%
    filter(adm_stagecat==1) %>% 
    summarize(n_distinct(ID))
10868/914

df %>% group_by(ID) %>% 
    mutate(adm_stagecat = stagecat[date==min(date)]) %>%
    filter(adm_stagecat==2) %>% 
    summarize(n_distinct(ID))
4922/203

# ERROR
#gen.inits = TRUE, 




#markovchain pkg
library(markovchain)

#markovchainFit(data = markov)
states = c("1", "2", "3", "4")
mctest <- new("markov", states = states, byrow = TRUE, generator = qmat0, name = "test")

plot(mctest, package = "diagram")


markovchainFit(data = markov, possibleStates = states)
