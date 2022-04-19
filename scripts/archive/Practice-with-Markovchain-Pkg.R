# Practice with `markovchain` package

library(markovchain)

# Set matrix
    # All rows or all columns sum to 1 
mat <- matrix(c(.25, .25, .25, .25,
                .25, .25, .25, .25,
                0, 0, 1, 0,
                0, 0, 0, 1),
              byrow=TRUE, 
              nrow=4, 
              ncol=4)

# Create Markov chain  
mc <- as(mat, "markovchain")

# Plot Markov chain
plot(mc, main="Covid Stages Markov Chain")

# Define probabilities of initial state 
init_state <- c(.5, .5, 0, 0)

# Calculate likelihood of each state after 2 days 
after2Days <- init_state * (mc^2)
after2Days

mc^2

# Print states of Markov chain
states(mc)

# Number of states of dimensions
dim(mc)

# Print transition matrix
t(mc)

# Obtain transition probabilities
    # Must provide initial state (t0) and subsequent state (t1)
transitionProbability(mc, t0 = "s1", t1 = "s2")

# Convert Markov chain to a data frame
as(mc, "data.frame")

# Shows the closed, recurrent, and transient classes. 
# Show absorbing states 
summary(mc)

# Vector of steady states in matrix form
steadyStates(mc)

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
SCC <- clusters(mc.igraph, mode="strong")
V(mc.igraph)$color <- rainbow(SCC$no)[SCC$membership]

# Shows communicating classes or strongly connected components 
plot(mc.igraph, mark.groups = split(1:vcount(mc.igraph), SCC$membership),
     main="Communicating classes - strongly connected components")

# Preview first 20 rows of stagecat variable
df$stagecat[1:20]

# Before next step, need to split data into lists by ID 
df <- split(df, f = df$ID)
df 

# Obtain empirical transition matrix 
createSequenceMatrix(stringchar = df$stagecat)  # this doesnt work on lists 
#lapply(df, createSequenceMatrix)
#try for loop
for (i in length(df)) {
    mat <- createSequenceMatrix(stringchar = df[[i]]$stagecat)
    print(mat)
}
    
df[1]

#fitting the DTMC by bootstrap
McFitBootstrap <- markovchainFit(data = df$stagecat, method = "bootstrap", name = "Covid Staging")
McFitBootstrap

#fitting the DTMC by MLE
McFitMLE <- markovchainFit(data = df$stagecat, method = "mle", name = "Covid Staging")
McFitMLE


