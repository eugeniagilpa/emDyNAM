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







# Newton-Raphson -------------------


# auxEstimateNR <- function(seqTime, envirPrepro, formula) {
#   # for(i in 1:length(seqsEM)){ # TO DO: PARALLELIZE THIS
#
#   envirPrepro$seqTime <- seqTime
#
#   local(
#     {
#       netEvents <- defineNetwork(nodes = actDfnodes, matrix = net0) |>
#         linkEvents(changeEvents = seqTime, nodes = actDfnodes)
#
#       depEvents <- defineDependentEvents(
#         seqTime[seqTime$replace == replaceIndex, ],
#         nodes = actDfnodes,
#         defaultNetwork = netEvents,
#       )
#     },
#     envirPrepro
#   )
#
#
#   res <- estimate(
#     as.formula(formula),
#     estimationInit = list(
#       engine = "default_c",
#       initialParameters = parameters,
#       fixedParameters = fixedparameters,
#       maxIterations = 0,
#       initialDamping = 1, dampingIncreaseFactor = 1, dampingDecreaseFactor = 1,
#       startTime = initTime,
#       endTime = endTime
#     ),
#     verbose = FALSE,
#     progress = FALSE,
#     envir = envirPrepro
#   )
#
#
#   logLikelihood <- res$logLikelihood
#   score <- res$finalScore
#   informationMatrix <- res$finalInformationMatrix
#
#   return(list(
#     logLikelihood = logLikelihood,
#     score = score,
#     informationMatrix = informationMatrix
#   ))
# }


# auxEstimateNR_MC <- function(indexCore, splitIndicesPerCore, seqsTime = seqsTime,
#                              envirPrepro = envirPrepro, formula = formula) {
#   indicesCore <- splitIndicesPerCore[[indexCore]]
#   resParameters <- vector("list", length(indicesCore))
#   for (i in seq_along(indicesCore)) {
#     res[[i]] <- auxEstimateNR(seqsTime[[indicesCore[[i]]]], envirPrepro, formula)
#   }
#
#   return(res)
# }
#



#' newtonraphsonStepChoice
#'
#' @param x values of the parameters
#' @param se standard errors of the parameters
#' @param w values of the weights
#'
#' @return value of the averaged mean and the standard error
#' @export
#'
newtonraphsonStepChoice <- function(parameters, fixedparameters,seqsEM) {

  stats.new <- lapply(seqsEM, "[[", "newlogLikelihoodStats")

# browser()
  logLikelihoodCrea.new <- sapply(unlist(lapply(stats.new, "[", "resCrea"),
    recursive = FALSE
  ), "[", "logLikelihood")
  scoreCrea.new <- unlist(
    lapply(unlist(lapply(stats.new, "[", "resCrea"),
      recursive = FALSE
    ), "[", "finalScore"),
    recursive = FALSE
  )
  informationMatrixCrea.new <- unlist(
    lapply(
      unlist(lapply(stats.new, "[", "resCrea"),
        recursive = FALSE
      ), "[",
      "finalInformationMatrix"
    ),
    recursive = FALSE
  )
  logLikelihoodDel.new <- sapply(unlist(lapply(stats.new, "[", "resDel"),
    recursive = FALSE
  ), "[", "logLikelihood")
  scoreDel.new <- unlist(
    lapply(unlist(lapply(stats.new, "[", "resDel"),
      recursive = FALSE
    ), "[", "finalScore"),
    recursive = FALSE
  )
  informationMatrixDel.new <- unlist(
    lapply(
      unlist(lapply(stats.new, "[", "resDel"),
        recursive = FALSE
      ), "[",
      "finalInformationMatrix"
    ),
    recursive = FALSE
  )

  idunfixedComponentsCrea <- which(is.na(fixedparameters$Crea))
  idunfixedComponentsDel <- which(is.na(fixedparameters$Del))

  # Calculate the UPDATE distance taking into account the DAMPING
  scoreUnfixedCrea <- Reduce("+", lapply(scoreCrea.new, "[", idunfixedComponentsCrea)) /
    length(seqsEM)
  scoreUnfixedDel <- Reduce("+", lapply(scoreDel.new, "[", idunfixedComponentsDel)) /
    length(seqsEM)

  informationMatrixUnfixedCrea <- Reduce("+", lapply(
    informationMatrixCrea.new, "[",
    idunfixedComponentsCrea)) / length(seqsEM) -
    scoreUnfixedCrea %*% t(scoreUnfixedCrea) / (length(seqsEM)^2)

  informationMatrixUnfixedDel <- Reduce("+", lapply(
    informationMatrixDel.new, "[",
    idunfixedComponentsDel)) / length(seqsEM) -
    scoreUnfixedDel %*% t(scoreUnfixedDel) / (length(seqsEM)^2)

  inverseInformationUnfixedCrea <- try(
    solve(informationMatrixUnfixedCrea),
    silent = TRUE
  )

  if (inherits(inverseInformationUnfixedCrea, "try-error")) {
    # stop(
    #   "Matrix cannot be inverted;",
    #   " probably due to collinearity between parameters. CREA"
    # )
    lambdaCrea = scoreUnfixedCrea%*%scoreUnfixedCrea
    inverseInformationUnfixedCrea <- try(
      solve(dkCreaAux),
      silent = TRUE
    )
    dkCreaAux = informationMatrixUnfixedCrea%*% informationMatrixUnfixedCrea +
      lambdaCrea[1,1] * diag(1,nrow(informationMatrixUnfixedCrea))
    inverseInformationUnfixedCrea <- try(
      solve(dkCreaAux),
      silent = TRUE
    )
    updateCrea <- rep(0, length(parameters$Crea))
    updateCrea[idunfixedComponentsCrea] <-
      (inverseInformationUnfixedCrea %*% t(informationMatrixUnfixedCrea) %*% scoreUnfixedCrea)
  }else{
    updateCrea <- rep(0, length(parameters$Crea))
    updateCrea[idunfixedComponentsCrea] <-
      (inverseInformationUnfixedCrea %*% scoreUnfixedCrea)
  }

  inverseInformationUnfixedDel <- try(
    solve(informationMatrixUnfixedDel),
    silent = TRUE
  )

  if (inherits(inverseInformationUnfixedDel, "try-error")) {
    # stop(
    #   "Matrix cannot be inverted;",
    #   " probably due to collinearity between parameters. DEL"
    # )
    lambdaDel = scoreUnfixedDel%*%scoreUnfixedDel
    dkDelAux = informationMatrixUnfixedDel%*% informationMatrixUnfixedDel +
      lambdaDel[1,1] * diag(1,nrow(informationMatrixUnfixedDel))
    inverseInformationUnfixedDel <- try(
      solve(dkDelAux),
      silent = TRUE
    )
    updateDel <- rep(0, length(parameters$Del))
    updateDel[idunfixedComponentsDel] <-
      (inverseInformationUnfixedDel %*% t(informationMatrixUnfixedDel) %*% scoreUnfixedDel)
  }else{
    updateDel <- rep(0, length(parameters$Del))
    updateDel[idunfixedComponentsDel] <-
      (inverseInformationUnfixedDel %*% scoreUnfixedDel)
  }

  parameters$Crea <- parameters$Crea + updateCrea
  parameters$Del <- parameters$Del + updateDel

  return(list(parameters = parameters,inverseInformationUnfixedCrea = inverseInformationUnfixedCrea,
              inverseInformationUnfixedDel = inverseInformationUnfixedDel))

  # if (sum(logLikelihoodCrea.new) <= sum(logLikelihoodCrea.old) ||
  #     any(is.na(logLikelihoodCrea.new)) ||
  #     any(is.na(unlist(scoreCrea.new))) ||
  #     any(is.na(unlist(informationMatrixCrea.new)))) {
  #   # reset values
  #   logLikelihoodCrea.new <- logLikelihoodCrea.old
  #   parameters$Crea <- parameters.old$Crea
  #   scoreCrea.new <- scoreCrea.old
  #   informationMatrixCrea.new <- informationMatrixCrea.old
  #   minDampingFactorCrea <- minDampingFactorCrea * dampingIncreaseFactor
  # } else {
  #   logLikelihoodCrea.old <- logLikelihoodCrea.new
  #   parameters.old$Crea <- parameters$Crea
  #   scoreCrea.old <- scoreCrea.new
  #   informationMatrixCrea.old <- informationMatrixCrea.new
  #   dampingFactorCrea <- max(
  #     1,
  #     minDampingFactorCrea / dampingDecreaseFactor
  #   )
  # # }


  # if (sum(logLikelihoodDel.new) <= sum(logLikelihoodDel.old) ||
  #     any(is.na(logLikelihoodDel.new)) ||
  #     any(is.na(unlist(scoreDel.new))) ||
  #     any(is.na(unlist(informationMatrixDel.new)))) {
  #   # reset values
  #   logLikelihoodDel.new <- logLikelihoodDel.old
  #   parameters$Del <- parameters.old$Del
  #   scoreDel.new <- scoreDel.old
  #   informationMatrixDel.new <- informationMatrixDel.old
  #   minDampingFactorDel <- minDampingFactorDel * dampingInDecreaseFactor
  # } else {
  #   logLikelihoodDel.old <- logLikelihoodDel.new
  #   parameters.old$Del <- parameters$Del
  #   scoreDel.old <- scoreDel.new
  #   informationMatrixDel.old <- informationMatrixDel.new
  #   minDampingFactorDel <- max(
  #     1,
  #     minDampingFactorDel / dampingDecreaseFactor
  #   )
  # # }






  # envirPrepro <- new.env()

  # seqsTime <- seqsEM

  # if (!"time" %in% colnames(seqsTime[[1]])) {
  #   seqsTime <- lapply(seqsTime, function(x) cbind(x, "time" = 1:nrow(x)))
  # }

  # # envirPrepro$seqsTime <- seqsEM
  # envirPrepro$actDfnodes <- actDfnodes
  # envirPrepro$net0 <- net0
  #
  # envirPrepro$fixedparameters <- fixedparameters
  # envirPrepro$initTime <- initTime
  # envirPrepro$endTime <- endTime
  #
  #
  # if (modelType == "Crea") {
  #   envirPrepro$replaceIndex <- 1
  # } else {
  #   envirPrepro$replaceIndex <- 0
  # }


  # formula <- paste("depEvents ~", formula, sep = "")
  #
  # logLikelihood <- vector("numeric", length = length(seqsEM))
  # score <- vector("list", length = length(seqsEM))
  # informationMatrix <- vector("list", length = length(seqsEM))

  # while (TRUE) {
  # envirPrepro$parameters <- parameters
  #
  #
  # if (iIteration == 1) {
  #   splitIndicesPerCore2 <- splitIndices(length(seqsTime), num_cores_parameters)
  #   cl2 <- makeCluster(num_cores_parameters)
  #   on.exit(stopCluster(cl2))
  #   clusterEvalQ(cl2, {
  #     library(goldfish)
  #     library(matrixStats)
  #     NULL
  #   })
  #
  #
  #   clusterExport(cl2, list(
  #     "auxEstimateNR"
  #   ))
  # }
  #
  # res <- clusterApply(cl2, seq_along(splitIndicesPerCore2), auxEstimateNR_MC,
  #   seqsTime = seqsTime,
  #   splitIndicesPerCore = splitIndicesPerCore2,
  #   envirPrepro = envirPrepro, formula = formula
  # )


  # logLikelihood <- sapply(res, "[[", 1)
  # score <- lapply(res, "[[", 2)
  # informationMatrix <- lapply(res, "[[", 3)



  # isInitialEstimation <- FALSE


  # # check for stop criteria
  # if (max(abs(score)) <= maxScoreStopCriterion) {
  #   isConverged <- TRUE
  #   if (progress) {
  #     cat(
  #       "\nStopping as maximum absolute score is below ",
  #       maxScoreStopCriterion, ".\n",
  #       sep = ""
  #     )
  #   }
  #
  #   break
  # }
  # if (iIteration > maxIterations) {
  #   if (progress) {
  #     cat(
  #       "\nStopping as maximum of ",
  #       maxIterations,
  #       " iterations have been reached. No convergence.\n"
  #     )
  #   }
  #   break
  # }


  # } # end of while
}





#' newtonraphsonStepRate
#'
#' @param x values of the parameters
#' @param se standard errors of the parameters
#' @param w values of the weights
#'
#' @return value of the averaged mean and the standard error
#' @export
#'
newtonraphsonStepRate <- function(parameters, seqsEM) {
  initialDamping <- 1
# browser()
  stats.new <- lapply(seqsEM, "[[", "newloglikRate")


  logLikelihoodCrea.new <- sapply(unlist(lapply(stats.new, "[", "resCrea"),
                                         recursive = FALSE
  ), "[", "logLikelihood")
  scoreCrea.new <- unlist(
    lapply(unlist(lapply(stats.new, "[", "resCrea"),
                  recursive = FALSE
    ), "[", "score"),
    recursive = FALSE
  )
  informationMatrixCrea.new <- unlist(
    lapply(
      unlist(lapply(stats.new, "[", "resCrea"),
             recursive = FALSE
      ), "[",
      "informationMatrix"
    ),
    recursive = FALSE
  )
  logLikelihoodDel.new <- sapply(unlist(lapply(stats.new, "[", "resDel"),
                                        recursive = FALSE
  ), "[", "logLikelihood")
  scoreDel.new <- unlist(
    lapply(unlist(lapply(stats.new, "[", "resDel"),
                  recursive = FALSE
    ), "[", "score"),
    recursive = FALSE
  )
  informationMatrixDel.new <- unlist(
    lapply(
      unlist(lapply(stats.new, "[", "resDel"),
             recursive = FALSE
      ), "[",
      "informationMatrix"
    ),
    recursive = FALSE
  )



  # Calculate the UPDATE distance taking into account the DAMPING

  scoreUnfixedCrea <- Reduce("+", scoreCrea.new) /
    length(seqsEM)
  scoreUnfixedDel <- Reduce("+", scoreDel.new) /
    length(seqsEM)

  informationMatrixUnfixedCrea <- Reduce("+",
    informationMatrixCrea.new) / length(seqsEM) -
    scoreUnfixedCrea %*% t(scoreUnfixedCrea) / (length(seqsEM)^2)

  informationMatrixUnfixedDel <- Reduce("+",
    informationMatrixDel.new) / length(seqsEM) -
    scoreUnfixedDel %*% t(scoreUnfixedDel) / (length(seqsEM)^2)


  inverseInformationUnfixedCrea <- try(
    solve(informationMatrixUnfixedCrea),
    silent = TRUE
  )
  if (inherits(inverseInformationUnfixedCrea, "try-error")) {
    lambdaCrea = scoreUnfixedCrea%*%scoreUnfixedCrea
    dkCreaAux = informationMatrixUnfixedCrea%*% informationMatrixUnfixedCrea +
      lambdaCrea[1,1] * diag(1,nrow(informationMatrixUnfixedCrea))
    inverseInformationUnfixedCrea <- try(
      solve(dkCreaAux),
      silent = TRUE
    )
    updateCrea <- rep(0, length(parameters$Crea))
    updateCrea<-
      (inverseInformationUnfixedCrea %*% t(informationMatrixUnfixedCrea) %*% scoreUnfixedCrea)
  }else{
    updateCrea <- rep(0, length(parameters$Crea))
    updateCrea <-
      (inverseInformationUnfixedCrea %*% scoreUnfixedCrea)
  }

  inverseInformationUnfixedDel <- try(
    solve(informationMatrixUnfixedDel),
    silent = TRUE
  )
  if (inherits(inverseInformationUnfixedDel, "try-error")) {
    lambdaDel = scoreUnfixedDel%*%scoreUnfixedDel
    dkDelAux = informationMatrixUnfixedDel%*% informationMatrixUnfixedDel +
      lambdaDel[1,1] * diag(1,nrow(informationMatrixUnfixedDel))
    inverseInformationUnfixedDel <- try(
      solve(dkDelAux),
      silent = TRUE
    )
    updateDel <- rep(0, length(parameters$Del))
    updateDel <-
      (inverseInformationUnfixedDel %*% t(informationMatrixUnfixedDel) %*% scoreUnfixedDel)
  }else{
    updateDel <- rep(0, length(parameters$Del))
    updateDel <-
      (inverseInformationUnfixedDel %*% scoreUnfixedDel)
  }


  parameters$Crea <- parameters$Crea + updateCrea
  parameters$Del <- parameters$Del + updateDel


  return(list(parameters = parameters,inverseInformationUnfixedCrea = inverseInformationUnfixedCrea,
              inverseInformationUnfixedDel = inverseInformationUnfixedDel))


}


