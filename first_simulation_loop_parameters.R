rm(list=ls())
# packages ----------------------------------------------------------------
library(goldfish)
library(broom)
library(pixiedust)
library(survival)
library(mlogit)
library(RSiena)

# Start of the loop to store the values of parameters
# Simulation parameters
endTime <- 30
nSimCrea <- 59 # from s501 and s502 data, number of ties created between t_0 and t_1
nSimDel <- 56 # from s501 and s502 data, number of ties deleted between t_0 and t_1

actDf <-dimnames(s501)[[2]]
nAct <- length(actDf)

# parameters = c(log(nSim/endTime/avgActors),...), avgActors = actors that can create/delete ties in t_0 
parmsCrea <- c(log(nSimCrea / endTime / nAct )) #indegree and outdegree of ego
parmsDel <- c(log(nSimDel / endTime / (nAct-4))) #indegree and outdegree of ego
# there are 4 actors with no ties, so no deletion is possible

parmsChoiceCrea <- c(1,0.05,0.1) #reciprocity, indegree, outdegree of alter (reciprocity positiva para creacion/ negativa eliminacion)
parmsChoiceDel <- c(-1,-0.05,-0.1)     

n = 100
storeParam = list(
  "mod1Crea"=data.frame("indeg"=numeric(n), "outdeg"=numeric(n), "recip"=numeric(n)),
  "mod2Crea"=data.frame("indeg"=numeric(n), "outdeg"=numeric(n), "recip"=numeric(n)),
  "mod1Del"=data.frame("indeg"=numeric(n), "outdeg"=numeric(n), "recip"=numeric(n)),
  "mod2Del"=data.frame("indeg"=numeric(n), "outdeg"=numeric(n), "recip"=numeric(n))
)

repeatedAction = list()

for(i in 1:n){
    # simulate DyNAM ------------------------------------------------------
    # Xstate = matrix(0L, nrow=nAct, ncol=nAct)#Xstate=s501
    Xstate = s501
    rownames(Xstate) = dimnames(s501)[[2]]
    
    XAux <- cbind(colSums(Xstate),rowSums(Xstate))
    dimnames(XAux) <- list(actDf,c("indeg","outdeg"))
    
    X <- cbind("intercept"=rep(1,nrow(XAux)))
    
    # Xchoice[sender, receiver, effect]
    Xchoice <- array(
      c(t(Xstate),
        matrix(XAux[,"indeg"],nrow=nAct,ncol=nAct,byrow=TRUE),
        matrix(XAux[,"outdeg"],nrow=nAct,ncol=nAct,byrow=TRUE)),
      dim = c(nAct, nAct, 3),
      dimnames = list(actDf, actDf, c( "recip", "indeg","outdeg"))
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
    #events = data.frame(time=numeric(200),label=character(200))|>as.list()
    
    while (time < endTime) {
      cat("event:", event, "time:", time, "\r")
      timeEventCrea <- rexp(1, sumRateCrea)
      timeEventDel <- rexp(1, sumRateDel)
    
      if(timeEventCrea < timeEventDel){
        # In this case, we have a creation event
        timeEvent = timeEventCrea
        
        # The nodes that are available to create a tie are those with an outdeg < nAct
        isAvailable <-XAux[,"outdeg"]<(nAct-1)
        sender <- sample(sum(isAvailable), 1, prob = expXbCrea[isAvailable] / sumRateCrea)
        labSender <- actDf[isAvailable][sender]
        
        # The nodes that are available to receive a tie are: 
        isAvailableChoice = which(Xstate[labSender,]==0 & actDf != labSender)
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
        Xchoice[labReceiver,labSender,"recip"] = 1
        Xchoice[,labReceiver,"indeg"] = Xchoice[,labReceiver,"indeg"]+1
        Xchoice[,labSender,"outdeg"] = Xchoice[,labSender,"outdeg"]+1
        
      }else{
        # In this case, we have a deletion event
        timeEvent = timeEventDel
        # The nodes that are available to delete a tie are those with an outdeg > 0
        isAvailable <-XAux[,"outdeg"]>0
        sender <- sample(sum(isAvailable), 1, prob = expXbDel[isAvailable] / sumRateDel)
        labSender <- actDf[isAvailable][sender]
        
        # The nodes that are available to break a tie with sender are: 
        isAvailableChoice = which(Xstate[labSender,]==1) 
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
        Xchoice[labReceiver,labSender,"recip"] = 0
        Xchoice[,labReceiver,"indeg"] = Xchoice[,labReceiver,"indeg"]-1
        Xchoice[,labSender,"outdeg"] = Xchoice[,labSender,"outdeg"]-1
        
      }
    }
    
    ### repeated events -----------------------------------------------------
    tt = table(events$sender,events$receiver)
    t = which(tt>1,arr.ind=TRUE)
    repeatedAction[[i]] = data.frame("sender"=rownames(tt)[t[,1]],
                                     "receiver"=colnames(tt)[t[,2]],
                                     "n"=tt[t])
    
    ### goldfish --------------------------------------------------------------
    actDfnodes <- data.frame(label=actDf,present = TRUE)
    actDfnodes <- defineNodes(actDfnodes) 
    
    eventsMod <- events[, -5]
    
    initState <- s501
    dimnames(initState) <- list(actDfnodes$label, actDfnodes$label)
    
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
  
    modTest1 <- estimate(
      depEventsCrea ~ indeg(netEvents, weighted = FALSE) +
        outdeg(netEvents, weighted = FALSE) +
        recip(netEvents),
      model = "DyNAM", subModel = "choice",
      estimationInit = list(
        startTime = 0,
        opportunitiesList = supportConstrainCrea,
        returnIntervalLogL = TRUE,
        returnEventProbabilities = TRUE
      ),
      progress = FALSE
    )
    storeParam$mod1Crea[i,] = coef(modTest1)
    
    modTest2 <- estimate(
      depEventsCrea ~ indeg + outdeg + recip + inertia,
      model = "DyNAM", subModel = "choice",
      estimationInit = list(
        startTime = 0,
        fixedParameters = c(NA, NA, NA, -20),
        returnIntervalLogL = TRUE,
        returnEventProbabilities = TRUE
      ),
      progress = FALSE
    )
    storeParam$mod2Crea[i,] = coef(modTest2)
    
    
    modTest3 <- estimate(
      depEventsDel ~ indeg(netEvents, weighted = FALSE) +
        outdeg(netEvents, weighted = FALSE) +
        recip(netEvents),
      model = "DyNAM", subModel = "choice",
      estimationInit = list(
        startTime = 0,
        opportunitiesList = supportConstrainDel,
        returnIntervalLogL = TRUE,
        returnEventProbabilities = TRUE
      ),
      progress = FALSE
    )
    storeParam$mod1Del[i,] = coef(modTest3)
    
    modTest4 <- estimate(
      depEventsDel ~ indeg + outdeg + recip + inertia,
      model = "DyNAM", subModel = "choice",
      estimationInit = list(
        startTime = 0,
        fixedParameters = c(NA, NA, NA, 20),
        returnIntervalLogL = TRUE,
        returnEventProbabilities = TRUE
      ),
      progress = FALSE
    )
    storeParam$mod2Del[i,] = coef(modTest4)
}


str(storeParam)
str(repeatedAction)
table(sapply(repeatedAction,nrow))
sapply(repeatedAction,function(x) table(x$n))



par(mfrow=c(1,2))
hist(storeParam$mod1Crea$indeg,main="Mod1 creation, indeg")
hist(storeParam$mod2Crea$indeg,main="Mod2 creation,indeg")
hist(storeParam$mod1Del$indeg,main="Mod1 deletion, indeg")
hist(storeParam$mod2Del$indeg,main="Mod2 deletion, indeg")
hist(storeParam$mod1Crea$outdeg,main="Mod1 creation, outdeg")
hist(storeParam$mod2Crea$outdeg,main="Mod2 creation,outdeg")
hist(storeParam$mod1Del$outdeg,main="Mod1 deletion, outdeg")
hist(storeParam$mod2Del$outdeg,main="Mod2 deletion, outdeg")
hist(storeParam$mod1Crea$recip,main="Mod1 creation, recip")
hist(storeParam$mod2Crea$recip,main="Mod2 creation,recip")
hist(storeParam$mod1Del$recip,main="Mod1 deletion, recip")
hist(storeParam$mod2Del$recip,main="Mod2 deletion, recip")

sapply(storeParam$mod1Crea,median)
sapply(storeParam$mod2Crea,mean)
sapply(storeParam$mod1Del,median)
sapply(storeParam$mod2Del,mean)


which(storeParam$mod1Crea$recip< (-6))
which(storeParam$mod2Crea$recip< (-6))

which(storeParam$mod1Del$recip< (-3))
which(storeParam$mod2Del$recip< (-3))
