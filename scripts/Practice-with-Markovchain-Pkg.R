# Practice with `markovchain` package

library(markovchain)

mat <- matrix(c(.25, .25, .25, .25,
                .25, .25, .25, .25,
                0, 0, 1, 0,
                0, 0, 0, 1),
              byrow=TRUE, 
              nrow=4, 
              ncol=4)

mc <- as(mat, "markovchain")

plot(mc, main="Covid Stages Markov Chain")

init_state <- c(.5, .5, 0, 0)

after2Days <- init_state * (mc^2)
after2Days

mc^2

states(mc)
dim(mc)
t(mc)
transitionProbability(mc)

as(mc, "data.frame")

summary(mc)

steadyStates(mc)
transientStates(mc)
absorbingStates(mc)

library(igraph)

mc.igraph <- as(mc, "igraph")

#finding and formatting the clusters
SCC <- clusters(mc.igraph, mode="strong")
V(mc.igraph)$color <- rainbow(SCC$no)[SCC$membership]

#plotting
plot(mc.igraph, mark.groups = split(1:vcount(mc.igraph), SCC$membership),
     main="Communicating classes - strongly connected components")

#load df
df$stagecat[1:20]

#before next step, need to split data into lists by ID 
df <- split(df, f = df$ID)
df 

#obtaining the empirical transition matrix 
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


