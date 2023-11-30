rm(list=ls())
# packages ----------------------------------------------------------------
library(goldfish)
library(broom)
library(pixiedust)
library(survival)
library(mlogit)
library(RSiena)
library(matrixStats)
library(stats)
library(foreach)
library(doParallel)

library(ggraph)
library(migraph)

set.seed(1)

# Single simulation ------------------------
endTime <- 200

parmsChoiceCrea <- c(0.6,0,0,0) #indegree, outdegree,reciprocity, transitivity of alter (reciprocity positiva para creacion/ negativa eliminacion)
parmsChoiceDel <- c(-1,0,0,0)     

  
nSimCrea <- 59 # from s501 and s502 data, number of ties created between t_0 and t_1
nSimDel <- 56 # from s501 and s502 data, number of ties deleted between t_0 and t_1

actDf <-dimnames(s501)[[2]]
nAct <- length(actDf)


# parameters = c(log(nSim/endTime/avgActors),...), avgActors = actors that can create/delete ties in t_0 
parmsCrea <- c(log(nSimCrea / 20 / nAct )) #indegree and outdegree of ego
parmsDel <- c(log(nSimDel / 20 / (nAct-4))) #indegree and outdegree of ego
# there are 4 actors with no ties, so no deletion is possible


# simulate DyNAM ------------------------------------------------------
# Xstate = matrix(0L, nrow=nAct, ncol=nAct)#Xstate=s501

Xstate = s501
rownames(Xstate) = dimnames(s501)[[2]]
initState = Xstate

XAux <- cbind(colSums(Xstate),rowSums(Xstate))
dimnames(XAux) <- list(actDf,c("indeg","outdeg"))

XAuxInit = XAux

X <- cbind("intercept"=rep(1,nrow(XAux)))

# Xchoice[sender, receiver, effect]
Xchoice <- array(
  c(matrix(XAux[,"indeg"],nrow=nAct,ncol=nAct,byrow=TRUE),
    matrix(XAux[,"outdeg"],nrow=nAct,ncol=nAct,byrow=TRUE),
    t(Xstate),
    Xstate%*%Xstate),
  dim = c(nAct, nAct, 4),
  dimnames = list(actDf, actDf, c( "recip", "indeg","outdeg","trans"))
)

# initial info
expXbCrea <- exp(X %*% parmsCrea)
sumRateCrea <- sum(expXbCrea)
expXbDel <- exp(X %*% parmsDel)
sumRateDel <- sum(expXbDel)

# init simulation:
time <- 0
events <- NULL
supportConstrain <- list()
event <- 0L
receiver <- 0L
num_event <- 0
#events = data.frame(time=numeric(200),label=character(200))|>as.list()

while (time < endTime) {
# while (num_event < tot_events + 1){
  stop = FALSE
  
  cat("event:", event, "time:", time, "\r")
  timeEventCrea <- rexp(1, sumRateCrea)
  timeEventDel <- rexp(1, sumRateDel)
  
  if(timeEventCrea < timeEventDel){
    # In this case, we have a creation event
    timeEvent = timeEventCrea
    
    # The nodes that are available to create a tie are those with an outdeg < nAct
    isAvailable <-XAux[,"outdeg"]<(nAct-1) & XAuxInit[,"outdeg"]<(nAct-1)
    if(sum(isAvailable)<1) next
    sender <- sample(sum(isAvailable), 1, prob = expXbCrea[isAvailable] / sumRateCrea)
    labSender <- actDf[isAvailable][sender]
    
    # The nodes that are available to receive a tie are: 
    isAvailableChoice = which(Xstate[labSender,]==0 & actDf != labSender & initState[labSender,]==0)
    counter=0
    
    while(length(isAvailableChoice)<1){ #redraw sender
      counter = counter +1
      if(counter>sum(isAvailable)){stop =TRUE; break}
      isAvailable <-XAux[,"outdeg"]>0 & XAuxInit[,"outdeg"]>0 & !actDf==labSender
      sender <- sample(sum(isAvailable), 1, prob = expXbDel[!(rownames(XAux) %in% labSender)][isAvailable] / sumRateDel)
      labSender <- actDf[isAvailable][sender]
      isAvailableChoice = which(Xstate[labSender,]==0 & actDf != labSender & initState[labSender,]==0)
    }
    if(stop) next
    
    
    choiceSet <- Xchoice[labSender,isAvailableChoice,]
    expUtil <- exp(choiceSet %*% parmsChoiceCrea) 
    receiver <- sample(length(isAvailableChoice), 1, prob = expUtil / sum(expUtil))
    labReceiver <- names(isAvailableChoice)[receiver]
    
    time <- time + timeEvent
    event <- event + 1L
    events <- rbind(
      events,
      data.frame(
        time = time,
        sender = labSender,
        receiver = labReceiver,
        replace = 1,
        timeElapse = timeEvent
      )
    )
    supportConstrain[[event]]= isAvailableChoice
    Xstate[labSender,labReceiver]=1
    XAux[labSender,"outdeg"] = XAux[labSender,"outdeg"]+1
    XAux[labReceiver,"indeg"] = XAux[labReceiver,"indeg"]+1
    # Xchoice[labReceiver,labSender,"recip"] = 1
    # Xchoice[,labReceiver,"indeg"] = Xchoice[,labReceiver,"indeg"]+1
    # Xchoice[,labSender,"outdeg"] = Xchoice[,labSender,"outdeg"]+1
    
    Xchoice <- array(
      c(matrix(XAux[,"indeg"],nrow=nAct,ncol=nAct,byrow=TRUE),
        matrix(XAux[,"outdeg"],nrow=nAct,ncol=nAct,byrow=TRUE),
        t(Xstate),
        Xstate%*%Xstate),
      dim = c(nAct, nAct, 4),
      dimnames = list(actDf, actDf, c("indeg","outdeg","recip","trans"))
    )
    num_event = num_event+1
  }else{
    # In this case, we have a deletion event
    timeEvent = timeEventDel
    # The nodes that are available to delete a tie are those with an outdeg > 0
    isAvailable <-XAux[,"outdeg"]>0 & XAuxInit[,"outdeg"]>0 
    if(sum(isAvailable)<1) next
    sender <- sample(sum(isAvailable), 1, prob = expXbDel[isAvailable] / sumRateDel)
    labSender <- actDf[isAvailable][sender]
    
    # The nodes that are available to break a tie with sender are: 
    isAvailableChoice = which(Xstate[labSender,]==1 & initState[labSender,]==1)
    
    counter=0
    while(length(isAvailableChoice)<1){ #redraw sender
      counter = counter +1
      if(counter>sum(isAvailable)){stop =TRUE; break}
      isAvailable <-XAux[,"outdeg"]>0 & XAuxInit[,"outdeg"]>0 & !actDf==labSender
      sender <- sample(sum(isAvailable), 1, prob = expXbDel[!(rownames(XAux) %in% labSender)][isAvailable] / sumRateDel)
      labSender <- actDf[isAvailable][sender]
      isAvailableChoice = which(Xstate[labSender,]==1 & initState[labSender,]==1)
    }
    if(stop) next
    
    
    choiceSet <- Xchoice[labSender,isAvailableChoice,]
    expUtil <- exp(choiceSet %*% parmsChoiceDel) 
    receiver <- sample(length(isAvailableChoice), 1, prob = expUtil / sum(expUtil))
    labReceiver <- names(isAvailableChoice)[receiver]
    
    time <- time + timeEvent
    event <- event + 1L
    events <- rbind(
      events,
      data.frame(
        time = time,
        sender = labSender,
        receiver = labReceiver,
        replace = 0,
        timeElapse = timeEvent
      )
    )
    supportConstrain[[event]]= isAvailableChoice
    Xstate[labSender,labReceiver]=0
    XAux[labSender,"outdeg"] = XAux[labSender,"outdeg"]-1
    XAux[labReceiver,"indeg"] = XAux[labReceiver,"indeg"]-1
    # Xchoice[labReceiver,labSender,"recip"] = 0
    # Xchoice[,labReceiver,"indeg"] = Xchoice[,labReceiver,"indeg"]-1
    # Xchoice[,labSender,"outdeg"] = Xchoice[,labSender,"outdeg"]-1
    
    Xchoice <- array(
      c(matrix(XAux[,"indeg"],nrow=nAct,ncol=nAct,byrow=TRUE),
        matrix(XAux[,"outdeg"],nrow=nAct,ncol=nAct,byrow=TRUE),
        t(Xstate),
        Xstate%*%Xstate),
      dim = c(nAct, nAct, 4),
      dimnames = list(actDf, actDf, c( "indeg","outdeg","recip","trans"))
    )
    num_event = num_event+1
  }
}

events
actDfnodes <- data.frame(label=actDf,present = TRUE)
actDfnodes <- defineNodes(actDfnodes) 

eventsMod <- events[, -5]

initState <- s501
dimnames(initState) <- list(actDfnodes$label, actDfnodes$label)

initNet <- defineNetwork(nodes = actDfnodes, matrix = initState)

netEvents <- defineNetwork(nodes = actDfnodes,matrix=initState) |>
  linkEvents(changeEvents = eventsMod, nodes = actDfnodes)

depEventsCrea <- defineDependentEvents(
  eventsMod[eventsMod$replace == 1,],
  nodes = actDfnodes, 
  defaultNetwork = netEvents
)
depEventsDel <- defineDependentEvents(
  eventsMod[eventsMod$replace == 0,],
  nodes = actDfnodes, 
  defaultNetwork = netEvents
)

supportConstrainCrea <- supportConstrain[eventsMod$replace == 1]
supportConstrainDel <- supportConstrain[eventsMod$replace == 0]


modTest2 <- estimate(
  depEventsCrea ~ indeg + outdeg + recip + trans + inertia + tie(initNet),
  model = "DyNAM", subModel = "choice",
  estimationInit = list(
    startTime = 0,
    fixedParameters = c(NA, NA, NA, NA, -20, -20),
    returnIntervalLogL = TRUE,
    returnEventProbabilities = TRUE
  ),
  progress = FALSE
)


modTest4 <- estimate(
  depEventsDel~ indeg + outdeg + recip + trans + inertia + tie(initNet),
  model = "DyNAM", subModel = "choice",
  estimationInit = list(
    startTime = 0,
    fixedParameters = c(NA, NA, NA, NA, 20, 20),
    returnIntervalLogL = TRUE,
    returnEventProbabilities = TRUE
  ),
  progress = FALSE
)





# The network at the beginning
callNetworkBgn <- as.matrix(netEvents)
autographr(callNetworkBgn, labels = FALSE, layout = "fr")

callNetworkHlf <- as.matrix(netEvents,time = 100)
autographr(callNetworkHlf, labels = FALSE, layout = "fr")


callNetworkEnd <- as.matrix(netEvents,time = 200)
autographr(callNetworkEnd, labels = FALSE, layout = "fr")

table(as.matrix(netEvents))
table(as.matrix(netEvents, time = 200))





coefOriginal = list(creation = parmsChoiceCrea ,deletion = parmsChoiceDel)

save(modTest2, modTest4, coefOriginal, initState, Xstate, events,file = "simulationData.RData")
