#' Function that generates simulations of the process of competition between
#' creation and deletion of ties.
#' Indegree, outdegree, reciprocity and transitivity in the choice model.
#' Constant rate.
#'
#' @param nSimCrea number of Creation events to be created (in average)
#' @param nSimDel number of Deletion events to be created (in average)
#' @param parmsChoiceCrea values of the Creation choice parameters (indeg,
#'                        outdeg, recip, trans)
#' @param parmsChoiceDel values of the Deletion choice parameters (indeg,
#'                        outdeg, recip, trans)
#' @param endTime time between networks, set to 1 by definition
#' @param rep allow for repetition of events with same sender-receiver pair
#' @param init_matrix adjacency matrix of initial network
#'
#' @return List with results of estimation given by goldfish, initial and final
#'          networks (adjacency matrix), sequence of events (eventsMod)
#'
#' @example library(RSiena)
#' parmsChoiceCrea = c(0.2,0.2,0.6,0.6)
#' parmsChoiceDel = c(-0.6,-0.6,-2,-2)
#' nSimCrea = 50
#' nSimDel = 50
#' simulation = simulationDyNAM(nSimCrea = nSimCrea,
#'                              nSimDel = nSimDel,
#'                              parmsChoiceCrea = parmsChoiceCrea,
#'                              parmsChoiceDel = parmsChoiceDel,
#'                              endTime = 1, rep = FALSE, init_matrix = s501)
#'
#' @export
#'
#'
simulationDyNAM = function(nSimCrea,nSimDel,parmsChoiceCrea,parmsChoiceDel,
                           endTime = 1, rep = FALSE,init_matrix = s501){

  actDf <-dimnames(init_matrix)[[2]]
  nAct <- length(actDf)

  # Rate parameters
  parmsCrea <- c(log(nSimCrea / endTime / nAct ))
  parmsDel <- c(log(nSimDel / endTime / (nAct-4)))


  # simulate DyNAM ------------------------------------------------------

  Xstate = init_matrix
  rownames(Xstate) = dimnames(init_matrix)[[2]]
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
  events <- NULL
  supportConstrain <- list()
  event <- 0L
  receiver <- 0L
  num_event <- 0
  time = 0L

  while (time < endTime) {
    # browser()
    if(!rep){
      stop = FALSE
      # browser()
      # cat("event:", event, "time:", time, "\r")
      timeEventCrea <- rexp(1, sumRateCrea)
      timeEventDel <- rexp(1, sumRateDel)

      if(timeEventCrea < timeEventDel){
        # In this case, we have a creation event
        timeEvent = timeEventCrea

        # The nodes that are available to create a tie are those with an outdeg < nAct
        isAvailable <-XAux[,"outdeg"]<(nAct-1) & XAuxInit[,"outdeg"]<(nAct-1)
        # cat('sum is avarilable crea',sum(isAvailable),"\n")
        if(sum(isAvailable)<1){next}
        sender <- sample(sum(isAvailable), 1, prob = expXbCrea[isAvailable] / sumRateCrea)
        labSender <- actDf[isAvailable][sender]

        # The nodes that are available to receive a tie are:
        isAvailableChoice = which(Xstate[labSender,]==0 & actDf != labSender & initState[labSender,]==0)
        # cat('sum is avarilable choice',sum(isAvailableChoice),"\n")
        counter=0
        labSenderPast = labSender
        while(length(isAvailableChoice)<1){ #redraw sender
          counter = counter +1
          if(counter>sum(isAvailable)){stop =TRUE; break}
          isAvailable <-XAux[,"outdeg"]<(nAct-1) & XAuxInit[,"outdeg"]<(nAct-1) & !actDf==labSenderPast
          # cat('sum is avarilable redraw',sum(isAvailable),"\n")
          if(sum(isAvailable)<1){next}
          sender <- sample(sum(isAvailable), 1, prob = expXbCrea[isAvailable] / sumRateCrea)
          labSender <- actDf[isAvailable][sender]
          labSenderPast <- c(labSender,labSenderPast)
          isAvailableChoice = which(Xstate[labSender,]==0 & actDf != labSender & initState[labSender,]==0)
          # cat('sum is avarilable choice',sum(isAvailableChoice),"\n")
        }
        if(stop){next}


        choiceSet <- Xchoice[labSender,isAvailableChoice,]
        expUtil <- exp(choiceSet %*% parmsChoiceCrea)
        # cat(expUtil,"exputil Crea \n")
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
        # cat('sum is avarilable del',sum(isAvailable),"\n")
        if(sum(isAvailable)<1) next
        sender <- sample(sum(isAvailable), 1, prob = expXbDel[isAvailable] / sumRateDel)
        labSender <- actDf[isAvailable][sender]

        # The nodes that are available to break a tie with sender are:
        isAvailableChoice = which(Xstate[labSender,]==1 & initState[labSender,]==1)
        # cat('sum is avarilable choice del',sum(isAvailableChoice),"\n")
        counter=0
        labSenderPast = labSender
        while(length(isAvailableChoice)<1){ #redraw sender
          counter = counter +1
          if(counter>sum(isAvailable)){stop =TRUE; break}
          isAvailable <-XAux[,"outdeg"]>0 & XAuxInit[,"outdeg"]>0 & !actDf==labSenderPast
          # cat('sum is avarilable redraw del',sum(isAvailable),"\n")
          # browser()
          # cat("prob", expXbDel[isAvailable] / sumRateDel,"\n")
          if(sum(isAvailable)<1){next}
          sender <- sample(sum(isAvailable), 1, prob = expXbDel[isAvailable] / sumRateDel)
          labSender <- actDf[isAvailable][sender]
          labSenderPast = c(labSender,labSenderPast)
          isAvailableChoice = which(Xstate[labSender,]==1 & initState[labSender,]==1)
          # cat('sum is avarilable choice',sum(isAvailableChoice),"\n")
        }
        if(stop){next}


        choiceSet <- Xchoice[labSender,isAvailableChoice,]
        expUtil <- exp(choiceSet %*% parmsChoiceDel)
        # cat(expUtil," exp util \n")
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
    }else{

      stop = FALSE
      # browser()
      # cat("event:", event, "time:", time, "\r")
      timeEventCrea <- rexp(1, sumRateCrea)
      timeEventDel <- rexp(1, sumRateDel)

      if(timeEventCrea < timeEventDel){
        # In this case, we have a creation event
        timeEvent = timeEventCrea

        # The nodes that are available to create a tie are those with an outdeg < nAct
        isAvailable <-XAux[,"outdeg"]<(nAct-1)
        if(sum(isAvailable)<1) next
        sender <- sample(sum(isAvailable), 1, prob = expXbCrea[isAvailable] / sumRateCrea)
        labSender <- actDf[isAvailable][sender]

        # The nodes that are available to receive a tie are:
        isAvailableChoice = which(Xstate[labSender,]==0 & actDf != labSender)
        counter=0

        while(length(isAvailableChoice)<1){ #redraw sender
          counter = counter +1
          if(counter>sum(isAvailable)){stop =TRUE; break}
          isAvailable <-XAux[,"outdeg"]<(nAct-1) & !actDf==labSender
          sender <- sample(sum(isAvailable), 1, prob = expXbCrea[!(rownames(XAux) %in% labSender)][isAvailable] / sumRateCrea)
          labSender <- actDf[isAvailable][sender]
          isAvailableChoice = which(Xstate[labSender,]==0 & actDf != labSender )
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
        isAvailable <-XAux[,"outdeg"]>0
        if(sum(isAvailable)<1) next
        sender <- sample(sum(isAvailable), 1, prob = expXbDel[isAvailable] / sumRateDel)
        labSender <- actDf[isAvailable][sender]

        # The nodes that are available to break a tie with sender are:
        isAvailableChoice = which(Xstate[labSender,]==1 )

        counter=0
        while(length(isAvailableChoice)<1){ #redraw sender
          counter = counter +1
          if(counter>sum(isAvailable)){stop =TRUE; break}
          isAvailable <-XAux[,"outdeg"]>0  & !actDf==labSender
          sender <- sample(sum(isAvailable), 1, prob = expXbDel[!(rownames(XAux) %in% labSender)][isAvailable] / sumRateDel)
          labSender <- actDf[isAvailable][sender]
          isAvailableChoice = which(Xstate[labSender,]==1 )
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

  }


  actDfnodes <- data.frame(label=actDf,present = TRUE)
  actDfnodes <- defineNodes(actDfnodes)

  eventsMod <- events[, -5]

  initState <- init_matrix
  dimnames(initState) <- list(actDfnodes$label, actDfnodes$label)

  initNet <- defineNetwork(nodes = actDfnodes, matrix = initState)


  envirPreproCrea <- new.env()
  envirPreproDel <- new.env()


  envirPreproCrea$eventsMod <- eventsMod
  envirPreproCrea$actDfnodes <- actDfnodes
  envirPreproCrea$initNet <- initNet

  envirPreproDel$eventsMod <- eventsMod
  envirPreproDel$actDfnodes <- actDfnodes
  envirPreproDel$initNet <- initNet


  local(
    {
      netEvents <- defineNetwork(nodes = actDfnodes,matrix=initState) |>
        linkEvents(changeEvents = eventsMod, nodes = actDfnodes)

      depEventsCrea <- defineDependentEvents(
        eventsMod[eventsMod$replace == 1,],
        nodes = actDfnodes,
        defaultNetwork = netEvents
      )
    },
    envirPreproCrea
  )

  local(
    {
      netEvents <- defineNetwork(nodes = actDfnodes,matrix=initState) |>
        linkEvents(changeEvents = eventsMod, nodes = actDfnodes)

      depEventsDel <- defineDependentEvents(
        eventsMod[eventsMod$replace == 0,],
        nodes = actDfnodes,
        defaultNetwork = netEvents
      )
    },
    envirPreproDel
  )


  modCrea <- estimate(
    depEventsCrea ~ indeg + outdeg + recip + trans + inertia + tie(initNet),
    model = "DyNAM", subModel = "choice",
    estimationInit = list(
      startTime = 0,
      fixedParameters = c(NA, NA, NA, NA, -20, -20),
      returnIntervalLogL = TRUE,
      returnEventProbabilities = TRUE
    ),
    progress = FALSE,
    envir = envirPreproCrea
  )

  modDel <- estimate(
    depEventsDel~ indeg + outdeg + recip + trans + inertia + tie(initNet),
    model = "DyNAM", subModel = "choice",
    estimationInit = list(
      startTime = 0,
      fixedParameters = c(NA, NA, NA, NA, 20, 20),
      returnIntervalLogL = TRUE,
      returnEventProbabilities = TRUE
    ),
    progress = FALSE,
    envir = envirPreproDel
  )

  eventsMod <- events[, -5]

  return(list("modCrea" = modCrea, "modDel" = modDel,"net0" = initState, "net1" = Xstate , "eventsMod" = eventsMod ))

}

