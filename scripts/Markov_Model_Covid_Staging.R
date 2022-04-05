

# Load required packages
library(tidyverse)
library(markovchain)

# Import and preview data
raw <-
    read_csv(
        "/Users/sreynolds2/Documents/GitHub/MS-Markov_Covid_Staging/data/dm_stg_clean_11.08.21.csv"
    )

# Transform stages to stage categories
# Stage 4-5 = 1 (Moderate)
# Stage 6-9 = 2 (Severe)
# Stage 10 = 3 (Dead)
# Discharged = 4 (Discharged)
df <- raw %>%
    mutate(stage.cat = case_when(stage %in% c(4, 5) ~ 1,
                                 stage %in% 6:9 ~ 2,
                                 stage == 10 ~ 3)) %>%
    select(ID, date, stage.cat) %>%
    group_by(ID) %>%
    do ({
        if (max(.$stage.cat) == 3) {
            df <- .
        } else {
            new <- tibble(
                ID = .$ID[1],
                date = max(.$date) + 1,
                stage.cat = 4
            )
            df <- rbind(., new)
        }
        df
    })

df_factor <- df %>%
    mutate(
        stage.cat = case_when(
            stage.cat == 1 ~ "Moderate",
            stage.cat == 2 ~ "Severe",
            stage.cat == 3 ~ "Dead",
            stage.cat == 4 ~ "Discharge"
        )
    )

# Determine stage transition matrix probabilities
stages_today_and_next <- df_factor %>%
    group_by(ID) %>%
    mutate(
        stage.cat.tmw = lead(stage.cat),
        stage.cat.tmw = ifelse(is.na(stage.cat.tmw), "NULL", stage.cat.tmw)
    )

# Count number of pt-days in each stage.cat
total_ptdays_per_stagecat <- stages_today_and_next %>%
    group_by(stage.cat) %>%
    count()

n_dead = total_ptdays_per_stagecat[1, 2]
n_disc = total_ptdays_per_stagecat[2, 2]
n_mod = total_ptdays_per_stagecat[3, 2]
n_sev = total_ptdays_per_stagecat[4, 2]

# Count number of pt-days categorized by each possible transition
# FROM stage.cat 1 TO 1, 2, 3, 4
# FROM stage.cat 2 TO 1, 2, 3, 4
# FROM stage.cat 3 TO 1, 2, 3, 4 -- This will be 0, since Dead is an absorbing state.
# FROM stage.cat 4 TO 1, 2, 3, 4 -- This will be 0, since Discharge is an absorbing state.

# Table showing type of stage.transition for each pt-day
# Categories include: Mod_to_Mod, Mod_to_Sev, Mod_to_Dead, Mod_to_Disc, Sev_to_Mod, Sev_to_Sev, Sev_to_Dead, Sev_to_Disc, Dead, Disc.
stage_transitions <- stages_today_and_next %>%
    mutate(stage.transition =
               as.factor(
                   case_when(
                       stage.cat == "Moderate" & stage.cat.tmw == "Moderate" ~ "Mod_to_Mod",
                       stage.cat == "Moderate" &
                           stage.cat.tmw == "Severe" ~ "Mod_to_Sev",
                       stage.cat == "Moderate" &
                           stage.cat.tmw == "Dead" ~ "Mod_to_Dead",
                       stage.cat == "Moderate" &
                           stage.cat.tmw == "Discharge" ~ "Mod_to_Disc",
                       stage.cat == "Severe" &
                           stage.cat.tmw == "Moderate" ~ "Sev_to_Mod",
                       stage.cat == "Severe" &
                           stage.cat.tmw == "Severe" ~ "Sev_to_Sev",
                       stage.cat == "Severe" &
                           stage.cat.tmw == "Dead" ~ "Sev_to_Dead",
                       stage.cat == "Severe" &
                           stage.cat.tmw == "Discharge" ~ "Sev_to_Disc",
                       stage.cat == "Dead" ~ "Dead",
                       stage.cat == "Discharge" ~ "Disc"
                   )
               ))

# Table showing number of pt-days for each possible stage transition
counts <- stage_transitions %>%
    group_by(stage.transition) %>%
    summarize(n = n()) %>%
    arrange(sub("Mod_to", "", stage.transition))

probs <- counts %>%
    summarize()

probs = c()

# prob = round(n/length(df$ID),4))
n_mod <- 34 + 987 + 8717 + 332
probs %>%
    summarize(sum(starts_with("Mod", stage.transition)))

stage_transitions %>%
    group_by(stage.transition) %>%
    summarize(n) %>%
    summarize(ifelse(starts_with("Mod", stage.transition), n(), NA))

# Frequency counts of each type of stage.transition
table(stage_transitions$stage.transition)

# Total number of patient-days
n_ptdays <- length(df$ID)
cat("Total number of patient-days =", sum(table(stage_transitions$stage.transition)))

# Alternative method would be to use ifelse and create new var/column for each type of stage.transition
stage_transitions_2 <- stages_today_and_next %>%
    mutate(
        Mod_to_Mod = ifelse(stage.cat == "Moderate" &
                                stage.cat.tmw == "Moderate", 1, 0),
        Mod_to_Sev = ifelse(stage.cat == "Moderate" &
                                stage.cat.tmw == "Severe", 1, 0),
        Mod_to_Dead = ifelse(stage.cat == "Moderate" &
                                 stage.cat.tmw == "Dead", 1, 0),
        Mod_to_Disc = ifelse(stage.cat == "Moderate" &
                                 stage.cat.tmw == "Discharge", 1, 0),
        Sev_to_Mod = ifelse(stage.cat == "Severe" &
                                stage.cat.tmw == "Moderate", 1, 0),
        Sev_to_Sev = ifelse(stage.cat == "Severe" &
                                stage.cat.tmw == "Severe", 1, 0),
        Sev_to_Dead = ifelse(stage.cat == "Severe" &
                                 stage.cat.tmw == "Dead", 1, 0),
        Sev_to_Disc = ifelse(stage.cat == "Severe" &
                                 stage.cat.tmw == "Discharge", 1, 0),
        Dead = ifelse(stage.cat == "Dead", 1, 0),
        Disc = ifelse(stage.cat == "Discharge", 1, 0),
        Moderate = ifelse(stage.cat == "Moderate", 1, 0),
        Severe = ifelse(stage.cat == "Severe", 1, 0)
    )

# Then, you can basically summarize all numeric columns -- you get the same results as above
stage_trans_counts <- stage_transitions_2 %>%
    summarize_if(is.numeric, sum) %>%
    select(-ID)

#write.csv(stage_trans_counts, "/Users/sreynolds2/Documents/GitHub/MS-Markov_Covid_Staging/stage_trans_counts.csv")

stage_trans_probs <- stage_trans_counts %>%
    summarize(
        prob_mm = sum(Mod_to_Mod) / sum(Moderate),
        prob_ms = sum(Mod_to_Sev) / sum(Moderate),
        prob_mdead = sum(Mod_to_Dead) / sum(Moderate),
        prob_mdisc = sum(Mod_to_Disc) / sum(Moderate),
        prob_sm = sum(Sev_to_Mod) / sum(Severe),
        prob_ss = sum(Sev_to_Sev) / sum(Severe),
        prob_sdead = sum(Sev_to_Dead) / sum(Severe),
        prob_sdisc = sum(Sev_to_Disc) / sum(Severe),
        prob_dead = sum(Dead) / sum(Dead),
        prob_disc = sum(Disc) / sum(Disc)
    )

# Moving along, we can load the `markovchain` package
library(markovchain)

# Print probs
stage_trans_probs

# Assign state names and probability matrix
states <- c("Moderate", "Severe", "Dead", "Discharge")
prob_mat <-
    matrix(
        data = c(
            stage_trans_probs[1],
            stage_trans_probs[2],
            stage_trans_probs[3],
            stage_trans_probs[4],
            stage_trans_probs[5],
            stage_trans_probs[6],
            stage_trans_probs[7],
            stage_trans_probs[8],
            0,
            0,
            stage_trans_probs[9],
            0,
            0,
            0,
            0,
            stage_trans_probs[10]
        ),
        byrow = TRUE,
        nrow = 4,
        dimnames = list(from = states,
                        to = states)
    )
# Attempt 2 - I think the first didn't work b/c need numeric where each row sums to 1
prob_mat_2 <- matrix(
    data = c(.866, .033, .003, .098,
             .078, .905, .011, .006,
             0, 0, 1, 0,
             0, 0, 0, 1),
    byrow = TRUE,
    nrow = 4,
    dimnames = list(from = states,
                    to = states)
)

# Create Markov Chain object
mc <- new(
    "markovchain",
    states = states,
    byrow = T,
    transitionMatrix = prob_mat_2,
    name = "Covid Markov"
)
mc

# Plot Markov chain
plot(mc, main = "Covid Stages Markov Chain")

# Define probabilities of initial state
init_state <- c(.5, .5, 0, 0)

# Calculate likelihood of each state after 2 days
after2Days <- init_state * (mc ^ 2)
after2Days

mc ^ 2

# Print states of Markov chain
states(mc)

# Number of states of dimensions
dim(mc)

# Print transition matrix
t(mc)

# Obtain transition probabilities
# Must provide initial state (t0) and subsequent state (t1)
transitionProbability(mc, t0 = "Moderate", t1 = "Discharge")

# Convert Markov chain to a data frame
as(mc, "data.frame")

# Shows the closed, recurrent, and transient classes.
# Show absorbing states
summary(mc)

# Transient states
transientStates(mc)

# Absorbing states
absorbingStates(mc)


################################################################################
#### igraph package ############################################################
################################################################################

# Load package
library(igraph)

mc.igraph <- as(mc, "igraph")

#finding and formatting the clusters
SCC <- clusters(mc.igraph, mode = "strong")
V(mc.igraph)$color <- rainbow(SCC$no)[SCC$membership]

# Shows communicating classes or strongly connected components
plot(mc.igraph,
     mark.groups = split(1:vcount(mc.igraph), SCC$membership),
     main = "Communicating classes - strongly connected components")

# Plot Markov Chain using igraph package
library(igraph)
g <- as(mc, "igraph")

V(g)["Moderate"]$color <- "red"
V(g)["Severe"]$color <- "green"
V(g)["Dead"]$color <- "blue"
V(g)["Discharge"]$color <- "purple"

plot(
    g,
    vertex.color = V(g)$color,
    vertex.size = 25,
    vertex.label.cex = 1.25,
    vertex.label.family = "Helvetica",
    vertex.label.dist = 3.5,
    edge.width = 2,
    edge.arrow.size = 1,
    edge.label = E(g)$prob,
    main = "Markov Chain for Predicting
     the Trajectory of COVID-19 Patients"
)
