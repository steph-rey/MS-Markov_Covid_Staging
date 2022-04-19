# Script to build function to create Markov chains 
# Author(s): Steph Reynolds (Stephanie.Reynolds@ucsf.edu)
#            Travis Porco (Travis.Porco@ucsf.edu)
#            Seth Blumberg (Seth.Blumberg@ucsf.edu)
# Date Created: 2022-04-15
# Date Modified: 2022-04-18 at 3:30PM PT 

# Load required packages 
library(tidyverse)

# Transition Intensity Matrix: p_ji = rate(to j from i) where i!=j 
    # p_11  p_12  p_13  
    # p_21  p_22  p_23 ...
    # p_31  p_32  p_33  
    #     ... 
    # Each col corresponds to indiv leaving that state

# 3-State Model 
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

# Build function to fix column sums of transition intensity matrix 
Fix.Col.Sum <- function(intensmat, eps=1e-8) {
    if (nrow(intensmat) != ncol(intensmat)) {
        stop("Fix.Col.Sum: Matrix not square")
    }
    diag(intensmat) <- rep(0, dim(intensmat)[1])
    colsum <- apply(intensmat,2,sum)  
    diag(intensmat) <- -colsum
    intensmat
}

# Build function to check if value is an integer
is.integer.value <- function(nn) {
    all((round(nn) - nn) == 0)
}

# Build function to generate Markov chain 
GenerateMC <- function(initstate, endtime, intensmat, maxit=65536, eps=1e-10) {
    if (!is.numeric(initstate) || initstate <= 0 || initstate > dim(intensmat)[2] || !is.integer.value(initstate)) {
        stop(paste0("At top of GenerateMC: Invalid state ", as.character(initstate)))
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
            stop(paste0("Invalid state ", as.character(curstate)))
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

# Assign values to transition intensity matrix 
TransIntensityMat2 <- structure(c(0, 0, 0, 0, 0, 0.00229847315711706, -0.130701034312921, 
                                  0, 0.0284025611558036, 0.1, 0.00666666666666667, 0, -0.159666666666667, 
                                  0.053, 0.1, 0.0107334525939177, 0, 0.0785330948121646, -0.0892665474060823, 
                                  0, 0, 0, 0, 0, 0), 
                                .Dim = c(5L, 5L), 
                                .Dimnames = list(c("D", "M", "MS", "S", "R"), c("D", "M", "MS", "S", "R")))

# Run function to check column sums of transition intensity matrix
Check.Col.Sum(TransIntensityMat2)  # Returns TRUE 

GenerateMC(initstate = 2, endtime = 500, intensmat = TransIntensityMat2)

# Build function to check event list element
Good.Event.List.Element <- function(elt) {
    class(elt) == "numeric" && length(elt) == 3 && names(elt)[1] == "from" && names(elt)[2] == "to" && names(elt)[3] == "event.time"
}

# Build function to check event list
Bad.Event.List <- function(event.list) {
    if (class(event.list) != "list") {
        return(FALSE)
    }
    if (length(event.list) == 0) {
        return(TRUE)
    }
    !all(sapply(event.list, Good.Event.List.Element))
}

# Build function to generate patient trajectory - returns initstate, endtime, and event.list
GenerateTrajectory <- function(initstate, endtime, intensmat) {
    list(initstate = initstate,
         endtime = endtime,
         event.list = GenerateMC(initstate = initstate, endtime = endtime, intensmat = intensmat))
}

# Build function to convert event.list to patient-day table - returns dataframe containing day and state 
GeneratePtDayTbl <- function(traj) {
    if (class(traj) != "list" || length(traj) != 3) {
        stop("Bad trajectory object")
    }
    if (names(traj)[1] != "initstate" || class(traj[[1]]) != "numeric") {
        stop("Bad trajectory object")
    }
    if (names(traj)[2] != "endtime" || class(traj[[2]]) != "numeric" || traj$endtime <= 0) {
        stop("Bad trajectory object")
    }
    if (names(traj)[3] != "event.list" || Bad.Event.List(traj$event.list)) {
        stop("Bad event list")
    }
    days <- 0:round(traj$endtime)
    if (length(traj$event.list) == 0) {
        data.frame(time = days,
                   state = rep(traj$initstate, length(days)))
    } else {
        inittrans <- data.frame(from = 0, 
                                to = traj$initstate, 
                                event.time = 0)
        elist <- rbind(inittrans, as.data.frame(do.call(rbind, traj$event.list)))
        tmptbl <- merge(elist, data.frame(event.time = days, keep = rep(TRUE, length(days))), by = "event.time", all.x = T, all.y = T)
        tmptbl$keep[is.na(tmptbl$keep)] <- FALSE
        tmptbl$from <- NULL
        names(tmptbl)[names(tmptbl) == "to"] <- "state"
        tmptbl <- tmptbl[order(tmptbl$event.time), ]
        if (nrow(tmptbl) == 1) {
            stop("Cannot happen")
        }
        for (ii in 2:nrow(tmptbl)) {
            if (is.na(tmptbl$state[ii])) {
                if (is.na(tmptbl$state[ii-1])) {
                    stop("Also cannot happen")
                }
                tmptbl$state[ii] <- tmptbl$state[ii-1]
            }
        }
        tmptbl <- tmptbl[tmptbl$keep, ]
        tmptbl$keep <- NULL
        names(tmptbl)[names(tmptbl) == "event.time"] <- "time"
        tmptbl[ , c("time", "state")]
    }
}

# Run function and view output 
GeneratePtDayTbl(xtemp)









