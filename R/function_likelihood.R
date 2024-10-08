# Log-likelihood computation function ----------------

#' log-likelihood computation
#'
#' @param listExpandedDF list returned from function GatherPreprocessingDF
#' @param beta list of choice parameters with elements Crea and Del
#'
#' @return loglikelihood (creation + deletion of events models)
#' @export
#'
logLikelihood <- function(listExpandedDFChoice, listExpandedDFRate, beta, theta,
                          initTime, endTime) {
  logLikChoice <- logLikChoice(listExpandedDFChoice, beta)
  logLikRate <- logLikRate(listExpandedDFRate, theta)
  logLikTime <- logLikTime(logLikRate, initTime, endTime)

  return(logLikChoice$loglik + logLikRate$loglik + logLikTime)
}

#' log-likelihood computation of the choice model
#'
#' @param listExpandedDF list returned from function GatherPreprocessingDF
#' @param beta list of choice parameters with elements Crea and Del
#'
#' @return loglikelihood (creation + deletion of events models)
#' @export
#'
logLikChoice <- function(listExpandedDF, beta) {
  # Computations of loglikelihood for creation events

  xbCrea <- by(
    listExpandedDF$expandedDFCreation[, !names(listExpandedDF$expandedDFCreation) %in% c("event", "selected")],
    listExpandedDF$expandedDFCreation[, "event"], function(x) as.matrix(x) %*% beta$Crea
  )
  xb_auxCrea <- by(
    listExpandedDF$expandedDFCreation[
      listExpandedDF$expandedDFCreation$selected == TRUE,
      !names(listExpandedDF$expandedDFCreation) %in% c("event", "selected")
    ],
    listExpandedDF$expandedDFCreation[listExpandedDF$expandedDFCreation$selected == TRUE, "event"],
    function(x) as.matrix(x) %*% beta$Crea
  )
  loglikCrea <- sum(unlist(xb_auxCrea) - sapply(xbCrea, colLogSumExps))

  # Computation of loglikelihood for deletion events
  xbDel <- by(
    listExpandedDF$expandedDFDeletion[, !names(listExpandedDF$expandedDFDeletion) %in% c("event", "selected")],
    listExpandedDF$expandedDFDeletion[, "event"], function(x) as.matrix(x) %*% beta$Del
  )
  xb_auxDel <- by(
    listExpandedDF$expandedDFDeletion[
      listExpandedDF$expandedDFDeletion$selected == TRUE,
      !names(listExpandedDF$expandedDFDeletion) %in% c("event", "selected")
    ],
    listExpandedDF$expandedDFDeletion[listExpandedDF$expandedDFDeletion$selected == TRUE, "event"],
    function(x) as.matrix(x) %*% beta$Del
  )
  loglikDel <- sum(unlist(xb_auxDel) - sapply(xbDel, colLogSumExps))


  return(list(
    "loglik" = loglikCrea + loglikDel, "xb_auxCrea" = xb_auxCrea,
    "xbCrea" = xbCrea, "xb_auxDel" = xb_auxDel, "xbDel" = xbDel
  ))
}






#' log-likelihood computation with time information
#'
#' @param listExpandedDF list returned from function GatherPreprocessingDF
#' @param theta list of rate parameters with elements Crea and Del
#'
#' @return loglikelihood (creation + deletion of events models)
#' @export
#'
logLikRate <- function(expandedDF, theta) {
  # Computations of loglikelihood for creation events
  xbCrea <- by(
    expandedDF$expandedDFCreation[, !names(expandedDF$expandedDFCreation) %in% c("event", "selected")],
    expandedDF$expandedDFCreation[, "event"], function(x) as.matrix(x) %*% theta$Crea
  )

  xbCrea <- as.data.frame(do.call(cbind, xbCrea))

  xb_auxCrea <- by(
    expandedDF$expandedDFCreation[
      expandedDF$expandedDFCreation$selected == TRUE,
      !names(expandedDF$expandedDFCreation) %in% c("event", "selected")
    ],
    expandedDF$expandedDFCreation[expandedDF$expandedDFCreation$selected == TRUE, "event"],
    function(x) as.matrix(x) %*% theta$Crea
  )
  loglikCrea <- sum(unlist(xb_auxCrea) - apply(xbCrea, 2, sum))


  # Computation of loglikelihood for deletion events
  xbDel <- by(
    expandedDF$expandedDFDeletion[, !names(expandedDF$expandedDFDeletion) %in% c("event", "selected")],
    expandedDF$expandedDFDeletion[, "event"], function(x) as.matrix(x) %*% theta$Del
  )

  xbDel <- as.data.frame(do.call(cbind, xbDel))

  xb_auxDel <- by(
    expandedDF$expandedDFDeletion[
      expandedDF$expandedDFDeletion$selected == TRUE,
      !names(expandedDF$expandedDFDeletion) %in% c("event", "selected")
    ],
    expandedDF$expandedDFDeletion[expandedDF$expandedDFDeletion$selected == TRUE, "event"],
    function(x) as.matrix(x) %*% theta$Del
  )
  loglikDel <- sum(unlist(xb_auxDel) - apply(xbDel, 2, sum))

  return(list(
    "loglik" = loglikCrea + loglikDel, "xb_auxCrea" = xb_auxCrea,
    "xbCrea" = xbCrea, "xb_auxDel" = xb_auxDel, "xbDel" = xbDel
  ))
}



#' log-likelihood computation with time information
#'
#' @param logLikRate list returned from logLik
#' @param initTime initial time of observation
#' @param endTime end time of observation
#'
#' @return loglikelihood (creation + deletion of events models)
#' @export
#'
logLikTime <- function(logLikRate, initTime, endTime) {
  # colLogSumExps!!! Todo en escala log
  lambdaCreaDF <- exp(logLikRate$xb_auxCrea - apply(logLikRate$xbCrea, 2, sum))
  lambdaDelDF <- exp(logLikRate$xb_auxDel - apply(logLikRate$xbDel, 2, sum))

  mu_alphaCrea <- sum(1 / lambdaCreaDF[-length(lambdaCreaDF)])
  sigma_alphaCrea <- sum(1 / (lambdaCreaDF[-length(lambdaCreaDF)])^2)
  logkappaCrea <- -log(lambdaCreaDF[length(lambdaCreaDF)]) -
    0.5 * log(2 * pi * sigma_alphaCrea) -
    (endTime - initTime - mu_alphaCrea)^2 / (2 * sigma_alphaCrea)

  mu_alphaDel <- sum(1 / lambdaDelDF[-length(lambdaDelDF)])
  sigma_alphaDel <- sum(1 / (lambdaDelDF[-length(lambdaDelDF)])^2)
  logkappaDel <- -log(lambdaDelDF[length(lambdaDelDF)]) -
    0.5 * log(2 * pi * sigma_alphaDel) -
    (endTime - initTime - mu_alphaDel)^2 / (2 * sigma_alphaDel)

  return(logkappaCrea + logkappaDel)
}



#' log-likelihood computation with time information
#'
#' TO DO : change value of \eqn{\pi_{i_r}} from 1/n to corresponding one.
#'
#' @param rate dataframe with constant rate for creation and deletion
#' @param initTime initial time of observation
#' @param endTime end time of observation
#' @param nactors number of actors
#' @param R lenth of the sequence

#'
#' @return loglikelihood (creation + deletion of events models)
#' (log-kappa with constant rate)
#' @export
#'
logLikConstantRate <- function(rate, initTime, endTime, nactors, R) {
  logkappaCrea <- nactors * rate$Crea * (endTime - initTime) +
    R * (log(nactors) + log(-rate$Crea) + log(endTime - initTime)) -
    log(factorial(R))

  logkappaDel <- nactors * rate$Del * (endTime - initTime) +
    R * (log(nactors) + log(-rate$Del) + log(endTime - initTime)) -
    log(factorial(R))

  logpiCrea <- -R * log(nactors)
  logpiDel <- -R * log(nactors)

  return(logkappaCrea + logkappaDel + logpiCrea + logpiDel)
}








#' log-likelihood computation: Multi-Core
#'
#' @description given parameters, sequences, formula, initial network and
#' multi-core elements, computes log-likelihood of all the sequences.
#'
#' @return vector of loglikelihoods
#' @export
#'
logLikelihoodMC <- function(indexCore, permut, beta, theta, initTime, endTime,
                            splitIndicesPerCore, actDfnodes = actDfnodes,
                            net0 = net0, formula = formula) {
  indicesCore <- splitIndicesPerCore[[indexCore]]
  resLikelihood <- vector("list", length(indicesCore))

  for (i in seq_along(indicesCore)) {
    seq <- permut[[indicesCore[[i]]]]
    envirPrepro <- new.env()
    seq$time <- seq(1, nrow(seq))
    envirPrepro$seqTime <- seq
    envirPrepro$actDfnodes <- actDfnodes
    envirPrepro$net0 <- net0

    local(
      {
        netEvents <- defineNetwork(nodes = actDfnodes, matrix = net0) |>
          linkEvents(changeEvents = seqTime, nodes = actDfnodes)

        depEvents <- defineDependentEvents(
          seqTime,
          nodes = actDfnodes, defaultNetwork = netEvents
        )
      },
      envirPrepro
    )
    listExpandedDF <- GatherPreprocessingDF(formula, envir = envirPrepro)$listExpandedDF

    resLikelihood[[i]] <- logLikChoice(listExpandedDF, beta)$loglik +
      logLikConstantRate(
        theta, initTime, endTime,
        length(actDfnodes$label), nrow(seq)
      )
  }
  return(unlist(resLikelihood))
}


#' log-likelihood computation: Multi-Core with temperatures (parallel tempering)
#'
#' @description given parameters, sequences, formula, initial network and
#' multi-core elements, computes log-likelihood of all the sequences.
#'
#' @return vector of loglikelihoods
#' @export
#'
logLikelihoodMCTemp <- function(indexCore, permut, beta, theta, initTime,
                                endTime, splitIndicesPerCore,
                                actDfnodes = actDfnodes, net0 = net0,
                                formula = formula, temp) {
  indicesCore <- splitIndicesPerCore[[indexCore]]
  resLikelihood <- vector("list", length(indicesCore))

  for (i in seq_along(indicesCore)) {
    seq <- permut[[indicesCore[[i]]]]
    envirPrepro <- new.env()
    seq$time <- seq(1, nrow(seq))
    envirPrepro$seqTime <- seq
    envirPrepro$actDfnodes <- actDfnodes
    envirPrepro$net0 <- net0

    local(
      {
        netEvents <- defineNetwork(nodes = actDfnodes, matrix = net0) |>
          linkEvents(changeEvents = seqTime, nodes = actDfnodes)

        depEvents <- defineDependentEvents(
          seqTime,
          nodes = actDfnodes, defaultNetwork = netEvents
        )
      },
      envirPrepro
    )
    listExpandedDF <- GatherPreprocessingDF(formula, envir = envirPrepro)

    resLikelihood[[i]] <- logLikChoice(listExpandedDF$listExpandedDF, beta)$loglik / temp +
      logLikConstantRate(
        theta, initTime, endTime,
        length(actDfnodes$label),
        nrow(seq)
      ) / temp
  }
  return(unlist(resLikelihood))
}











#' log-likelihood computation: Multi-Core with time
#'
#' @description given parameters, sequences, formula, initial network and
#' multi-core elements, computes log-likelihood of all the sequences.
#'
#' @return vector of loglikelihoods
#' @export
#'
logLikelihoodTimeMC <- function(indexCore, permut, beta, theta, initTime,
                                endTime, splitIndicesPerCore,
                                actDfnodes = actDfnodes, net0 = net0,
                                formulaChoice = formulaChoice,
                                formulaRate = formulaRate) {
  indicesCore <- splitIndicesPerCore[[indexCore]]
  resLikelihood <- vector("list", length(indicesCore))

  for (i in seq_along(indicesCore)) {
    seq <- permut[[indicesCore[[i]]]]
    envirPrepro <- new.env()
    seq$time <- seq(1, nrow(seq))
    envirPrepro$seqTime <- seq
    envirPrepro$actDfnodes <- actDfnodes
    envirPrepro$net0 <- net0

    local(
      {
        netEvents <- defineNetwork(nodes = actDfnodes, matrix = net0) |>
          linkEvents(changeEvents = seqTime, nodes = actDfnodes)

        depEvents <- defineDependentEvents(
          seqTime,
          nodes = actDfnodes, defaultNetwork = netEvents
        )
      },
      envirPrepro
    )
    listExpandedDFChoice <- GatherPreprocessingDF(formulaChoice, envir = envirPrepro)$listExpandedDF

    listExpandedDFRate <- GatherPreprocessingDF(formulaRate, envir = envirPrepro)$expandedDF

    resLikelihood[[i]] <- logLikelihood(
      listExpandedDFChoice, listExpandedDFRate,
      beta, theta, initTime, endTime
    )
  }
  return(unlist(resLikelihood))
}




#' log-likelihood computation: Multi-Core with time
#'
#' @description given parameters, sequences, formula, initial network and
#' multi-core elements, computes log-likelihood of all the sequences.
#'
#' @return vector of loglikelihoods
#' @export
#'
getlogLikelihood <- function(seq, actDfnodes, net0, fixedparameters,
                             parameters, initTime, endTime, formula) {
   # browser()
  envirPrepro <- new.env()
  if ("row" %in% colnames(seq)) {
    seqTime <- seq[, -which(colnames(seq) == "row")]
  } else {
    seqTime <- seq
  }
  seqTime <- cbind(seqTime, "time" = 1:nrow(seqTime))
  envirPrepro$seqTime <- seqTime
  envirPrepro$actDfnodes <- actDfnodes
  envirPrepro$net0 <- net0
  envirPrepro$fixedparameters <- fixedparameters
  # envirPrepro$initTime <- initTime
  # envirPrepro$endTime <- endTime
  envirPrepro$parameters <- parameters


  # CREATION

  envirPrepro$replaceIndex <- 1
  if ((length(unlist(strsplit(formula, "~"))) < 2)) {
    formula <- paste("depEvents ~", formula, sep = "")
  }

  local(
    {
      netEvents <- defineNetwork(nodes = actDfnodes, matrix = net0) |>
        linkEvents(changeEvents = seqTime, nodes = actDfnodes)

      depEvents <- defineDependentEvents(
        seqTime[seqTime$replace == replaceIndex, ],
        nodes = actDfnodes,
        defaultNetwork = netEvents,
      )
    },
    envirPrepro
  )

  resCrea <- estimate(
    as.formula(formula),
    estimationInit = list(
      engine = "default",
      initialParameters = parameters$Crea,
      fixedParameters = fixedparameters$Crea,
      maxIterations = 0,
      initialDamping = 1, dampingIncreaseFactor = 1, dampingDecreaseFactor = 1,
      # startTime = initTime,
      # endTime = endTime,
      onlyScore = TRUE
    ),
    verbose = FALSE,
    progress = FALSE,
    envir = envirPrepro
  )


  # DELETION
  envirPrepro$replaceIndex <- 0
  # formula <- paste("depEvents ~", formula, sep = "")

  local(
    {
      netEvents <- defineNetwork(nodes = actDfnodes, matrix = net0) |>
        linkEvents(changeEvents = seqTime, nodes = actDfnodes)

      depEvents <- defineDependentEvents(
        seqTime[seqTime$replace == replaceIndex, ],
        nodes = actDfnodes,
        defaultNetwork = netEvents,
      )
    },
    envirPrepro
  )

  resDel <- estimate(
    as.formula(formula),
    estimationInit = list(
      engine = "default",
      initialParameters = parameters$Del,
      fixedParameters = fixedparameters$Del,
      maxIterations = 0,
      initialDamping = 1, dampingIncreaseFactor = 1, dampingDecreaseFactor = 1,
      # startTime = initTime,
      # endTime = endTime,
      onlyScore = TRUE
    ),
    verbose = FALSE,
    progress = FALSE,
    envir = envirPrepro
  )


  return(list("resCrea" = resCrea, "resDel" = resDel))
}


#' log-likelihood computation: Multi-Core with time
#'
#' @description given parameters, sequences, formula, initial network and
#' multi-core elements, computes log-likelihood of all the sequences.
#'
#' @return vector of loglikelihoods
#' @export
#'
getlogLikelihoodMC <- function(indexCore, seqsEM, beta, fixedparameters,
                               splitIndicesPerCore, initTime, endTime,
                               actDfnodes = actDfnodes, net0 = net0,
                               formula = formula, temp) {
  indicesCore <- splitIndicesPerCore[[indexCore]]
  resLikelihood <- vector("list", length(indicesCore))
      # browser()
  for (i in seq_along(indicesCore)) {
    seq <- seqsEM[[indicesCore[i]]]$newseq
    if ("row" %in% colnames(seq)) {
      seq <- seq[, -which(colnames(seq) == "row")]
    }

    auxLik <- getlogLikelihood(
      seq = seq, actDfnodes = actDfnodes,
      net0 = net0,
      fixedparameters = fixedparameters,
      parameters = beta,
      initTime = initTime, endTime = endTime,
      formula = formula
    )
    resLikelihood[[i]] <- c("Crea" = auxLik$resCrea$logLikelihood, "Del" = auxLik$resDel$logLikelihood) /
      temp[indicesCore[i]]
  }

  return(resLikelihood)
}




#' log-likelihood computation: Multi-Core with time
#'
#' @description given parameters, sequences, formula, initial network and
#' multi-core elements, computes log-likelihood of all the sequences.
#'
#' @return vector of loglikelihoods
#' @export
#'
getlogLikelihoodRate <- function(seq, actDfnodes, parameters,
                                 initTime, endTime, temp) {

  n = length(actDfnodes$label)
  t = endTime - initTime
 # browser()

  # CREATION
    R = nrow(seq[seq$replace==1,])
    aux = n  *parameters$Crea * t
    logLikelihood = -aux+R*log(aux)-lfactorial(R)
    score = -n*t + R/parameters$Crea
    informationMatrix = R/(parameters$Crea)^2

    resCrea = list("logLikelihood" = logLikelihood/temp,
                "score" = score/temp,
                "informationMatrix" = informationMatrix/temp^2)
  # DELETION
    R = nrow(seq[seq$replace==0,])
    aux = n*parameters$Del*t
    logLikelihood = -aux+R*log(aux)-lfactorial(R)
    score = -n*t + R/parameters$Del
    informationMatrix = R/(parameters$Del)^2

    resDel = list("logLikelihood" = logLikelihood/temp,
                "score" = score/temp,
                "informationMatrix" = informationMatrix/temp^2)

  return(list("resCrea" = resCrea, "resDel" = resDel))
}



#' log-likelihood computation: Multi-Core with time
#'
#' @description given parameters, sequences, formula, initial network and
#' multi-core elements, computes log-likelihood of all the sequences.
#'
#' @return vector of loglikelihoods
#' @export
#'
getlogLikelihoodRateMC <- function(seqs, theta, initTime, endTime,
                               actDfnodes = actDfnodes, temp) {
    # browser()
resLikelihood=vector("list",length=length(seqs))
    for (i in 1:length(seqs)) {
    seq <- seqs[[i]]$newseq
    if ("row" %in% colnames(seq)) {
      seq <- seq[, -which(colnames(seq) == "row")]
    }

    auxLik <- getlogLikelihoodRate(
      seq = seq, actDfnodes = actDfnodes,
      parameters = theta,
      initTime = initTime, endTime = endTime,
      temp = temp[i]
    )
    resLikelihood[[i]] <- c("Crea" = auxLik$resCrea$logLikelihood,
                            "Del" = auxLik$resDel$logLikelihood)
  }

  return(resLikelihood)
}




