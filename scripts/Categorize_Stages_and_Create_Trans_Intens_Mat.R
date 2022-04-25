# This script: (1) Creates vector of initial states (based on stage_adm). 
# (2) Categorizes stages into M, MS, S, D, and R.
# (3) Creates a transition intensity matrix using msm::crudeinits function.

# Load required packages
library(msm)
library(slider) 
library(lubridate)
library(tidyverse)

# Import and preview data
dm_stg_clean <- read_csv("/Users/sreynolds2/Documents/GitHub/MS-Markov_Covid_Staging/data/dm_stg_clean_11.08.21.csv")
glimpse(dm_stg_clean)

# Create vector of initstates from data 
df_for_init_vec <- dm_stg_clean %>%
    distinct(ID, stage_adm) %>% 
    mutate(stage_adm.cat = case_when(stage_adm %in% 4:5 ~ 1,         # Moderate
                                     stage_adm %in% 6:9 ~ 2)) %>%    # Severe
    select(stage_adm.cat)

initstate.vec_covid <- df$stage_adm.cat
initstate.vec_covid

# Categorize patients into:
    # Moderate w/o hx of Severe (M=1)
    # Moderate WITH hx of Severe (MS=2)
    # Severe (S=3)
    # Dead (D=4)
df <- dm_stg_clean %>%
    mutate(stage_cat = case_when(stage %in% 4:5 ~ 1,  # Moderate
                                 stage %in% 6:9 ~ 3,  # Severe
                                 stage == 10 ~ 4),    # Dead 
           severe_yn = ifelse(stage_cat==3, 1, 0)) %>% 
    group_by(ID) %>% 
    mutate(hx_severe = slider::slide_index_sum(x = severe_yn, 
                                               i = date, 
                                               before = lubridate::days(365)),
           state = case_when(stage_cat==4 ~ 4,
                             stage_cat==3 ~ 3,
                             stage_cat==1 & hx_severe==0 ~ 1,
                             stage_cat==1 & hx_severe>0 ~ 2)) %>% 
    select(ID, date, state)   # Choose whether to include additional covariates (i.e. age) here

# If patient doesn't die, bind new row to represent recovered/discharge state
    # Recovered (R=5) 
df <- df %>% 
    do ({                        # If pt dies, leave as-is
        if (max(.$state) == 4) {    
            df <- .
        } else {                 # Else... create new row where state==5 to represent recovered 
            new <- tibble(ID = .$ID[1],
                          date = max(.$date) + 1,
                          state = 5)
            df <- rbind(., new)  # Bind to last row for each grouped ID
        }
        df
    })

df

# Now, we have a table of ID, date, state where state categorized into: 
    # 1 = Moderate (M)
    # 2 = Moderate WITH history of severe (MS)
    # 3 = Severe (S)
    # 4 = Dead (D)
    # 5 = Recovered (R)

# Create stage transition table 
statetbl <- statetable.msm(state, ID, data=df)

statetbl

# Assign matrix of allowed transitions, where 0=not allowed and 1=allowed
qmat <- matrix(c(1, 0, 1, 1, 1,
                 0, 1, 1, 1, 1,
                 0, 1, 1, 1, 1,
                 0, 0, 0, 1, 0,
                 0, 0, 0, 0, 1),
               nrow = 5, ncol = 5,
               byrow = T,
               dimnames = list(c("M", "MS", "S", "D", "R"),
                               c("M", "MS", "S", "D", "R")))

qmat

# Specify initial values and assign to `Q.crude`
Q.crude <- crudeinits.msm(state ~ date, ID, data=df, qmatrix=qmat)
Q.crude

# Build function to check row sums of transition intensity matrix 
Check.Row.Sum <- function(intensmat, eps=1e-8) {
    all(abs(apply(intensmat,1,sum) - 0) < eps) 
}

# Run function to check row sums
Check.Row.Sum(Q.crude)  # TRUE

# Using msm package, fit Markov model by maximum likelihood 
mc_fit <- msm(state ~ date, data = df, subject = ID, qmatrix = Q.crude)
mc_fit
