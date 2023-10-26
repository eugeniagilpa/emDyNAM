rm(list=ls())
# packages ----------------------------------------------------------------
library(goldfish)
library(broom)
library(pixiedust)
library(survival)
library(mlogit)
library(RSiena)
library(matrixStats)

set.seed(123)

# Single simulation ------------------------
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

  # simulate DyNAM ------------------------------------------------------
  # Xstate = matrix(0L, nrow=nAct, ncol=nAct)#Xstate=s501
  
  Xstate = s501
  rownames(Xstate) = dimnames(s501)[[2]]
  initState = Xstate
  
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
      isAvailable <-XAux[,"outdeg"]<(nAct-1) # añadir condiciones para tener en cuenta repeticion de eventos
      sender <- sample(sum(isAvailable), 1, prob = expXbCrea[isAvailable] / sumRateCrea)
      labSender <- actDf[isAvailable][sender]
      
      # The nodes that are available to receibe a tie are: 
      isAvailableChoice = which(Xstate[labSender,]==0  & actDf != labSender) #acutalizar dependiendo de events
      #isAvailableChoice = which(Xstate[labSender,]==0 & initState[labSender,]==0 & actDf != labSender) #acutalizar dependiendo de events
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
      isAvailable <-XAux[,"outdeg"]>0 # meter condicion para tener en cuenta repeticiones (out/ind initial states)
      sender <- sample(sum(isAvailable), 1, prob = expXbDel[isAvailable] / sumRateDel)
      labSender <- actDf[isAvailable][sender]
      
      # The nodes that are available to break a tie with sender are: 
      isAvailableChoice = which(Xstate[labSender,]==1)
      #isAvailableChoice = which(Xstate[labSender,]==1 & initState[labSender,]==1) 
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

  
  
# Setup of the EM algorithm ---------------
  
# We will start from two matrix or networks (initial and final) and a series 
# of events.
  
# Xstate and initState
# actors = rownames(Xstate)  
# sequence = which(initState!=Xstate,arr.ind = TRUE)
# colnames(sequence)  = c("sender","receiver")
# sequence = as.data.frame(sequence)
# sequence$replace = ifelse(initState[initState!=Xstate]==1,0,1)
# sequence$sender = actors[sequence$sender]
# sequence$receiver = actors[sequence$receiver]
# rownames(sequence) = 1:nrow(sequence)


EMPreprocessing = function(X0, X1){
  # Function that preprocess the network matrixes
  # X0 : initial network matrix (row, col names must be actor names)
  # X1 : final network matrix (row, col names must be actor names)
  actors = colnames(X0)  
  sequence = which(X0!=X1,arr.ind = TRUE)
  colnames(sequence)  = c("sender","receiver")
  sequence = as.data.frame(sequence)
  sequence$replace = ifelse(X0[X0!=X1]==1,0,1)
  sequence$sender = actors[sequence$sender]
  sequence$receiver = actors[sequence$receiver]
  rownames(sequence) = 1:nrow(sequence)
  return(sequence)
}

# sequence = EMPreprocessing(initState,Xstate)  
  
# Permutations -------------------------------------

  #nmax =1000, 10000 (to be revisited)
permute = function(x,nmax = 5){
  # x: original data frame with a column of senders, a column of receivers and a column of replace
  #    (replace = 0 means deletion of tie, replace = 1 means creation of tie)
  # nmax: maximum number of permutations performed. The number of permutations
  #       computed will be the min{length(x)!, nmax}
  # out: output of the function is a list of dataframes. Each dataframe is a permutation
  #      of the rows of the original list x
  n = nrow(x)
  if(factorial(n)<nmax){
    out = lapply(1:factorial(n),function(t) x[sample(n,n),])
  }else{
    out = lapply(1:nmax,function(t) x[sample(n,n),])
  }
  return(out)

  # if(factorial(n)<nmax){
  #   out = replicate(factorial(n), x[sample(n,n),],simplify = "list")
  # }else{
  #   out = replicate(nmax,x[sample(n,n),],simplify = "list")
  # }
  # return(out)
}

# permut = permute(sequence)
  
# Function to create time difference of events ---------

timeGenerator = function(seq,nAct,theta){
  # seq : sequence of events
  # nAct : number of actors
  # theta : list with Crea and Del objects with parameters for the sender
    
  X_aux <- cbind("intercept"=rep(1,nAct))
    
  expXbCrea <- exp(X_aux %*% theta$Crea)
  sumRateCrea <- sum(expXbCrea)
  expXbDel <- exp(X_aux %*% theta$Del)
  sumRateDel <- sum(expXbDel)
    
  time = numeric(nrow(seq))
  if(seq$replace[1]==0){
    time[1] =  rexp(1, sumRateDel)
  }else{
    time[1] = rexp(1, sumRateCrea)
  }
    
  for(i in 2:nrow(seq)){
    if(seq$replace[i]==0){
      time[i] = time[i-1] + rexp(1, sumRateDel)
    }else{
      time[i] = time[i-1] + rexp(1, sumRateCrea)
    }
  }
  return(time)
}
  
 
GatherPreprocessingDF = function(formula,seqTime,nodes=actDfnodes,matrix = net0){
  # formula : string formula with the effects for the model
  
  netEvents = defineNetwork(nodes = actDfnodes,matrix=net0) |>
      linkEvents(changeEvents = seqTime, nodes = actDfnodes)
    
  depEvents = defineDependentEvents(seqTime, nodes = actDfnodes, defaultNetwork = netEvents)
  
  formGather = paste("depEvents ~",formula,sep="")
  
  dataProcessed <- GatherPreprocessing(
    as.formula(formGather),
    model = "DyNAM", subModel = "choice",
    progress = FALSE
  )
  namesEffects <- dataProcessed$namesEffects
  nEvents <- length(dataProcessed$n_candidates)
  
  expandedDF <- cbind(
    setNames(
      as.data.frame(dataProcessed$stat_all_events),
      namesEffects
    ),
    data.frame(
      event = rep(seq.int(nEvents), dataProcessed$n_candidates),
      selected = sequence(dataProcessed$n_candidates) ==
        rep(dataProcessed$selected, dataProcessed$n_candidates),
      sender = rep(dataProcessed$sender, dataProcessed$n_candidates)
    )
  
  )
  
  nEvents = expandedDF$event[length(expandedDF$event)]
  indexEvents = data.frame("start"=seq(1,length(expandedDF$event),nAct-1),
                           "end"=seq(nAct-1,length(expandedDF$event),nAct-1),
                           "event"=seq(1,nEvents))
  creation_events = 0
  deletion_events = 0
  for(i in 1:nEvents){
  # Select rows of expandedDF corresponding to the event
  # Check if inertia == 0 (creation of ties) or inertia == 1 (delition of ties) for (receiver == TRUE)
   df_aux = expandedDF[seq(indexEvents$start[i],indexEvents$end[i]),]
      index_aux = df_aux[,"selected"]==TRUE
      if (df_aux[index_aux,"inertia_netEvents"]==0){ #creation of ties
         creation_events=c(creation_events,i)
      }else{
         deletion_events=c(deletion_events,i)
       }
    }
    creation_events=creation_events[-1]
    deletion_events=deletion_events[-1]
    
    expandedDFCreation = subset(expandedDF,event %in% creation_events & inertia_netEvents == 0)
    expandedDFDeletion = subset(expandedDF,event %in% deletion_events & inertia_netEvents == 1)
  
    listExpandedDF = list(
    expandedDFCreation = expandedDFCreation[,!names(expandedDFCreation) %in% c("inertia_netEvents", "sender")],
    expandedDFDeletion = expandedDFDeletion[,!names(expandedDFDeletion) %in% c("inertia_netEvents", "sender")])
    
  return(listExpandedDF)
}  
  
# Log-likelihood computation function ----------------

logLikelihood = function(listExpandedDF, beta){
  # listExpandedDF : list returned from function GatherePreprocessingDF
  # beta : choice parameters 
  
  # Computations of loglikelihood for creation events
  xbCrea = by(listExpandedDF$expandedDFCreation[,!names(listExpandedDF$expandedDFCreation) %in% c("event","selected")],
                  listExpandedDF$expandedDFCreation[,"event"],function(x) as.matrix(x)%*%beta$Crea)
  xb_auxCrea = by(listExpandedDF$expandedDFCreation[listExpandedDF$expandedDFCreation$selected==TRUE,!names(listExpandedDF$expandedDFCreation) %in% c("event","selected")],
              listExpandedDF$expandedDFCreation[listExpandedDF$expandedDFCreation$selected==TRUE,"event"],function(x) as.matrix(x)%*%beta$Crea)
  loglikCrea =  sum(unlist(xb_auxCrea) - sapply(xbCrea,colLogSumExps))

  # Computation of loglikelihood for deletion events
  xbDel = by(listExpandedDF$expandedDFDeletion[,!names(listExpandedDF$expandedDFDeletion) %in% c("event","selected")],
          listExpandedDF$expandedDFDeletion[,"event"],function(x) as.matrix(x)%*%beta$Del)
  xb_auxDel = by(listExpandedDF$expandedDFDeletion[listExpandedDF$expandedDFDeletion$selected==TRUE,!names(listExpandedDF$expandedDFDeletion) %in% c("event","selected")],
              listExpandedDF$expandedDFDeletion[listExpandedDF$expandedDFDeletion$selected==TRUE,"event"],function(x) as.matrix(x)%*%beta$Del)
  loglikDel = sum(unlist(xb_auxDel) - sapply(xbDel,colLogSumExps))
  
  return(loglikCrea+loglikDel)
  
}



# EM algorithm --------------------


EMAlgorithm = function(net0,net1,theta0,beta0,formula){
  # net0 : network matrix in initial state
  # net1: network matrix in final state
  # theta0 : initial parameters for sender (list of creation and deletion with 
  #          labels Crea and Del)
  # beta0: initial parameters for receiver (list of creation and deletion with 
  #        labels Crea and Del)
  # formula: string with formula for the model (only rhs of formula)

  actDfnodes <- defineNodes(data.frame(label=colnames(net0)))
  
  # Initialization of parameters
  theta = theta0
  beta = beta0
  
  
  for (i in 10){
  logLik = 0
  
  # Creation of sequence of events from initial data
  sequence = EMPreprocessing(net0,net1)
  
  # Creation of permutations
  permut = permute(sequence, nmax = 5)
  
  for(seq in permut){
   #seq=permut[[1]]
    # GatherPreprocessing ----------------

    seq$time = seq(1,nrow(seq))
    seqTime=seq
    #seqTime = data.frame("time"=timeGenerator(seq,nAct,theta),seq)

    # netEvents = defineNetwork(nodes = actDfnodes,matrix=net0) |>
    #  linkEvents(changeEvents = seqTime, nodes = actDfnodes)

    #depEvents = defineDependentEvents(seqTime, nodes = actDfnodes, defaultNetwork = netEvents)
    
    
    #formGather = paste("depEvents ~",formula,sep="")
    listExpandedDF = GatherPreprocessingDF(formula,seqTime,nodes=actDfnodes,matrix = net0)

    # LogLikelihood computation -------------
    logLik = c(logLik,logLikelihood(listExpandedDF,beta))
  }
  logLik = logLik[-1]
  
  
  # For the maximum value, we perform goldfish to estimate the parameters and update beta
  seqMax = permut[[which.max(logLik)]]
  seqMax$time = seq(1,nrow(seqMax))
  netEvents = defineNetwork(nodes = actDfnodes,matrix=net0) |>
    linkEvents(changeEvents = seqMax, nodes = actDfnodes)
  
  depEventsCrea <- defineDependentEvents(
    seqMax[seqMax$replace == 1,],
    nodes = actDfnodes, 
    defaultNetwork = netEvents
  )
  depEventsDel <- defineDependentEvents(
    seqMax[seqMax$replace == 0,],
    nodes = actDfnodes, 
    defaultNetwork = netEvents
  )
  
  supportConstrainCrea <- supportConstrain[seqMax$replace == 1]
  supportConstrainDel <- supportConstrain[seqMax$replace == 0]
  
  
  formCrea = paste("depEventsCrea ~",formula,sep="")
  modCrea <- estimate(
    as.formula(formCrea),
    model = "DyNAM", subModel = "choice",
    estimationInit = list(
      startTime = 0,
      fixedParameters = c(NA, NA, NA, -20),
      returnIntervalLogL = TRUE,
      returnEventProbabilities = TRUE
    ),
    progress = FALSE
  )
  beta$Crea = coef(modCrea)
  
  formDel = paste("depEventsDel ~",formula,sep="")
  modDel <- estimate(
      as.formula(formDel),
      model = "DyNAM", subModel = "choice",
      estimationInit = list(
        startTime = 0,
        fixedParameters = c(NA, NA, NA, 20),
        returnIntervalLogL = TRUE,
        returnEventProbabilities = TRUE
      ),
      progress = FALSE
    )
  beta$Del = coef(modDel)
  }
 
  return(beta)
   
}
  
net0=s501
rownames(net0)=colnames(net0)
net1=Xstate
theta0=data.frame("Crea"=parmsCrea,"Del"=parmsDel)
beta0=data.frame("Crea"=parmsChoiceCrea,"Del"=parmsChoiceDel)
formula = "indeg + outdeg + recip + inertia"

likelihood = EMAlgorithm(net0,net1,theta0,beta0,formula)




