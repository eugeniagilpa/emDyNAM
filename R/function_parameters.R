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
