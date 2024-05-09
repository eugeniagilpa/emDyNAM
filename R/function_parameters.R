# Parameter computation from goldfish ---------------------


#' Function that computes goldfish estimation of parameters for a given sequence
#'
#' @param seq sequence of values
#' @param actDfnodes. nodes
#' @param net0. initial network
#' @param formula. formula of the model
#'
#' @return list with goldfish estimators and their standard errors
#' @export
#'
parameters <- function(seq, actDfnodes. = actDfnodes, net0. = net0, formula. = formula) {
  envirPreproCrea <- new.env()
  envirPreproDel <- new.env()

  seq$time <- seq(1, nrow(seq))

  envirPreproCrea$seqTime <- seq
  envirPreproCrea$actDfnodes <- actDfnodes.
  envirPreproCrea$net0 <- net0.

  envirPreproDel$seqTime <- seq
  envirPreproDel$actDfnodes <- actDfnodes.
  envirPreproDel$net0 <- net0.


  local(
    {
      netEvents <- defineNetwork(nodes = actDfnodes, matrix = net0) |>
        linkEvents(changeEvents = seqTime, nodes = actDfnodes)

      depEventsCrea <- defineDependentEvents(
        seqTime[seqTime$replace == 1, ],
        nodes = actDfnodes,
        defaultNetwork = netEvents,
      )
    },
    envirPreproCrea
  )
  local(
    {
      netEvents <- defineNetwork(nodes = actDfnodes, matrix = net0) |>
        linkEvents(changeEvents = seqTime, nodes = actDfnodes)

      depEventsDel <- defineDependentEvents(
        seqTime[seqTime$replace == 0, ],
        nodes = actDfnodes,
        defaultNetwork = netEvents,
      )
    },
    envirPreproDel
  )

  formCrea <- paste("depEventsCrea ~", formula., sep = "")
  modCrea <- estimate(
    as.formula(formCrea),
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

  formDel <- paste("depEventsDel ~", formula., sep = "")
  modDel <- estimate(
    as.formula(formDel),
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
  return(list(
    betaCreaDF = coef(modCrea), betaDelDF = coef(modDel),
    seCreaDF = modCrea$standardErrors[modCrea$standardErrors > 0],
    seDelDF = modDel$standardErrors[modDel$standardErrors > 0]
  ))
}


#' Estimation of parameters: Multi-Core
#'
#' @description function that paralelizes computation of parameters on all the sequences using goldfish.
#'
#'
#' @return list with goldfish estimators and their standard errors for each permuation
#' @export
#'
parametersMC <- function(indexCore, splitIndicesPerCore, permut = permut,
                         actDfnodes. = actDfnodes, net0. = net0,
                         formula. = formula) {
  indicesCore <- splitIndicesPerCore[[indexCore]]
  resParameters <- vector("list", length(indicesCore))
  for (i in seq_along(indicesCore)) {
    resParameters[[i]] <- parameters(permut[[indicesCore[[i]]]], actDfnodes., net0., formula.)
  }

  return(resParameters)
}




# Rubin's rule for weighted data -------------

#' Rubin´s rule for weighted data
#'
#' @param x values of the parameters
#' @param se standard errors of the parameters
#' @param w values of the weights
#'
#' @return value of the averaged mean and the standard error
#' @export
#'
rubinsRule <- function(x, se, w) {
  x_mean <- apply(x, 2, weighted.mean, w)
  v_within <- apply(se^2, 2, weighted.mean, w)
  v_between <- apply((sweep(x, 2, x_mean))^2, 2, weighted.mean, w)
  se_total <- sqrt(v_within + (1 + 1 / length(w)) * v_between)
  return(list(mean = x_mean, se = se_total))
}







# Newton-Raphson --------------------


auxEstimateNR <- function(seqTime, envirPrepro, formula) {
  # for(i in 1:length(seqsEM)){ # TO DO: PARALLELIZE THIS

  envirPrepro$seqTime <- seqTime

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


  res <- estimate(
    as.formula(formula),
    estimationInit = list(
      engine = "default_c",
      initialParameters = parameters,
      fixedParameters = fixedparameters,
      maxIterations = 0,
      initialDamping = 1, dampingIncreaseFactor = 1, dampingDecreaseFactor = 1,
      startTime = initTime,
      endTime = endTime
    ),
    verbose = FALSE,
    progress = FALSE,
    envir = envirPrepro
  )


  logLikelihood <- res$logLikelihood
  score <- res$finalScore
  informationMatrix <- res$finalInformationMatrix

  return(list(
    logLikelihood = logLikelihood,
    score = score,
    informationMatrix = informationMatrix
  ))
}


auxEstimateNR_MC <- function(indexCore, splitIndicesPerCore, seqsTime = seqsTime,
                             envirPrepro = envirPrepro, formula = formula) {
  indicesCore <- splitIndicesPerCore[[indexCore]]
  resParameters <- vector("list", length(indicesCore))
  for (i in seq_along(indicesCore)) {
    res[[i]] <- auxEstimateNR(seqsTime[[indicesCore[[i]]]], envirPrepro, formula)
  }

  return(res)
}





estimateNR <- function(parameters, fixedparameters = c(NA, NA, NA, NA, -20, -20),
                       initTime, endTime, actDfnodes, net0,
                       formula, seqsEM, modelType = "Crea",
                       num_cores_parameters = 1) {
  initialDamping <- 1
  minDampingFactor <- initialDamping
  dampingIncreaseFactor <- 1
  dampingDecreaseFactor <- 1 # is it better 2 and 3 (defaults in goldfish) QUESTION
  iIteration <- 1
  isConverged <- FALSE
  isInitialEstimation <- TRUE
  logLikelihood.old <- -Inf
  parameters.old <- parameters
  score.old <- NULL
  informationMatrix.old <- NULL


  idUnfixedCompnents <- which(is.na(fixedparameters))

  envirPrepro <- new.env()

  seqsTime <- seqsEM

  if (!"time" %in% colnames(seqsTime[[1]])) {
    seqsTime <- lapply(seqsTime, function(x) cbind(x, "time" = 1:nrow(x)))
  }

  # envirPrepro$seqsTime <- seqsEM
  envirPrepro$actDfnodes <- actDfnodes
  envirPrepro$net0 <- net0

  envirPrepro$fixedparameters <- fixedparameters
  envirPrepro$initTime <- initTime
  envirPrepro$endTime <- endTime


  if (modelType == "Crea") {
    envirPrepro$replaceIndex <- 1
  } else {
    envirPrepro$replaceIndex <- 0
  }


  formula <- paste("depEvents ~", formula, sep = "")

  logLikelihood <- vector("numeric", length = length(seqsEM))
  score <- vector("list", length = length(seqsEM))
  informationMatrix <- vector("list", length = length(seqsEM))

  while (TRUE) {
    envirPrepro$parameters <- parameters


    if (iIteration == 1) {
      splitIndicesPerCore2 <- splitIndices(length(seqsTime), num_cores_parameters)
      cl2 <- makeCluster(num_cores_parameters)
      on.exit(stopCluster(cl2))
      clusterEvalQ(cl2, {
        library(goldfish)
        library(matrixStats)
        NULL
      })


      clusterExport(cl2, list(
        "auxEstimateNR"
      ))
    }

    res <- clusterApply(cl2, seq_along(splitIndicesPerCore2), auxEstimateNR_MC,
      seqsTime = seqsTime,
      splitIndicesPerCore = splitIndicesPerCore2,
      envirPrepro = envirPrepro, formula = formula
    )


    logLikelihood <- sapply(res, "[[", 1)
    score <- lapply(res, "[[", 2)
    informationMatrix <- lapply(res, "[[", 3)

    if (sum(logLikelihood) <= sum(logLikelihood.old) ||
      any(is.na(logLikelihood)) ||
      any(is.na(unlist(score))) ||
      any(is.na(unlist(informationMatrix)))) {
      # reset values
      logLikelihood <- logLikelihood.old
      parameters <- parameters.old
      score <- score.old
      informationMatrix <- informationMatrix.old
      minDampingFactor <- minDampingFactor * dampingIncreaseFactor
    } else {
      logLikelihood.old <- logLikelihood
      parameters.old <- parameters
      score.old <- score
      informationMatrix.old <- informationMatrix
      minDampingFactor <- max(
        1,
        minDampingFactor / ifelse(isInitialEstimation, 1, dampingDecreaseFactor) # QUESTIONS
      )
    }

    isInitialEstimation <- FALSE

    # Calculate the UPDATE distance taking into account the DAMPING
    dampingFactor <- minDampingFactor



    informationMatrixUnfixed <- Reduce("+", lapply(
      informationMatrix, "[",
      idUnfixedCompnents, idUnfixedCompnents
    )) /
      length(idUnfixedCompnents)

    scoreUnfixed <- Reduce("+", lapply(score, "[", idUnfixedCompnents)) /
      length(idUnfixedCompnents)

    inverseInformationUnfixed <- try(
      solve(informationMatrixUnfixed),
      silent = TRUE
    )
    if (inherits(inverseInformationUnfixed, "try-error")) {
      stop(
        "Matrix cannot be inverted;",
        " probably due to collinearity between parameters."
      )
    }

    update <- rep(0, nParams)
    update[idUnfixedCompnents] <-
      (inverseInformationUnfixed %*% scoreUnfixed) / dampingFactor



    # check for stop criteria
    if (max(abs(score)) <= maxScoreStopCriterion) {
      isConverged <- TRUE
      if (progress) {
        cat(
          "\nStopping as maximum absolute score is below ",
          maxScoreStopCriterion, ".\n",
          sep = ""
        )
      }

      break
    }
    if (iIteration > maxIterations) {
      if (progress) {
        cat(
          "\nStopping as maximum of ",
          maxIterations,
          " iterations have been reached. No convergence.\n"
        )
      }
      break
    }


    parameters <- parameters + update

    iIteration <- iIteration + 1
  } # end of while


  stdErrors <- rep(0, nParams)
  stdErrors[idUnfixedCompnents] <- sqrt(diag(inverseInformationUnfixed))

  return(list(parameters = parameters, stdErrors = stdErrors))
}
