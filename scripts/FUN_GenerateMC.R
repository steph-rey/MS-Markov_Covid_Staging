# Script to build function to create Markov chains 
# Author(s): Steph Reynolds (Stephanie.Reynolds@ucsf.edu)
#            Travis Porco (Travis.Porco@ucsf.edu)
#            Seth Blumberg (Seth.Blumberg@ucsf.edu)
# Date Created: 2022-04-15
# Date Modified: 2022-04-15 at 1:00PM PT 

# Load required packages 
library(tidyverse)

# Transition Intensity Matrix: p_ji = rate(to j from i) where i!=j 
    # p_11  p_12  p_13  
    # p_21  p_22  p_23 ...
    # p_31  p_32  p_33  
    #     ... 
    # Each col corresponds to indiv leaving that state

# 3-state Model 
    # A -> B, A -> C

# Assign values to transition intensity matrix 
TransIntensityMat1 <- matrix(c(-0.15,    0,    0,
                               0.1,      0,    0,
                               0.05,     0,    0),
                             nrow=3, 
                             ncol=3,
                             byrow=T)

# Build function to check column sums of transition intensity matrix 
Check.Col.Sum <- function(intensmat, eps=1e-8) {
    all(abs(apply(intensmat,2,sum) - 0) < eps) 
}

# Run function to check col sums 
Check.Col.Sum(TransIntensityMat1)

# Build function to check if integer
is.integer.value <- function(nn) {
    all((round(nn) - nn) == 0)
}

# Build function to generate Markov chain 
GenerateMC <- function(initstate, endtime, intensmat, maxit=65536, eps=1e-10) {
    if (!is.numeric(initstate) || initstate <= 0 || initstate > dim(intensmat)[2] || is.integer.value(initstate)) {
        stop("Invalid state")
    }
    if (!is.numeric(endtime) || endtime <= 0) {
        stop("Invalid end time")
    }
    etime <- 0
    #event_list <- tibble(etime = 0, fromstate = current_state)
    curit <- 0
    curstate <- initstate
    iplus <- intensmat[intensmat > 0]
    if (length(iplus) == 0) {
        stop("Bad intensity matrix")
    }
    tinyrate <- min(iplus) * 1e-10
    tinyrate <- ifelse(tinyrate <= 0, 1e-300, tinyrate)
    rec <- list()
    while(etime <= endtime) {
        if (curit > maxit) {
            stop("Max iterations exceeded")
        } else {
            curit <- curit + 1
        }
        if (curstate <= 0 || curstate > dim(intensmat)[2]) {
            stop("Invalid state")
        }
        flows <- intensmat[ , curstate]
        if (!all(is.finite(flows))) {
            stop("Bad intensities!")
        } 
        flows[curstate] <- 0    # Note: Can move this to beginning to make cleaner 
        if (all(abs(flows) < eps)) {
            break
        }
        flows[flows==0] <- tinyrate
        event.times <- rexp(length(flows), rate = flows)
        if (!all(is.finite(event.times))) {
            stop("Bad event times")
        }
        which.min(event.times) -> new.state
        etime + event.times[new.state] -> new.time
        rec[[length(rec) + 1]] <- c(from = curstate, to = new.state, event.time = new.time)
        if (new.time <= etime) {
            stop("Didn't go forward in time")
        }
        etime <- new.time
        curstate <- new.state
    }
    return(rec)
}

#GenerateMC <- function(initstate, endtime, intensmat, maxit=65536, eps=1e-10)
GenerateMC(initstate = , endtime = 500, intensmat = TransIntensityMat1)
