# EM preprocesing -------------------------

#' Auxiliary function for EM algorithm
#'
#' @param X0 initial network matrix (row, col names must be actor names)
#' @param X1 final network matrix (row, col names must be actor names)
#'
#' @return data frame with sequence of events (sender, receiver and replacement)
#' @export
#'
#' @examples library(RSiena)
#' X0 <- s501
#' X1 <- s502
#' EMPreprocessing(X0, X1)
EMPreprocessing <- function(X0, X1) {
  actors <- colnames(X0)
  sequence <- which(X0 != X1, arr.ind = TRUE)
  colnames(sequence) <- c("sender", "receiver")
  sequence <- as.data.frame(sequence)
  sequence$replace <- ifelse(X0[X0 != X1] == 1, 0, 1)
  sequence$sender <- actors[sequence$sender]
  sequence$receiver <- actors[sequence$receiver]
  rownames(sequence) <- 1:nrow(sequence)
  return(sequence)
}


# Gather preprocesing -------------------------

#' Gather preprocessing function
#'
#' TO DO: include Rate model / REM specification
#'
#' @param formula string formula with the effects for the model
#' @param envir environment were dependentEvents, nodes and net0 are located.
#'
#' @return listExpandedDF (see [goldfish::estimate()] for more)
#' @export
#'
GatherPreprocessingDF <- function(formula, envir = new.env(), submodel = "choice") {
  # formula : string formula with the effects for the model
  # envir : environment were dependentEvents, nodes and net0 are located.
  # return: listExpandedDF

  ##################
  # model = rem

  # browser()
  if ((length(unlist(strsplit(formula, "~"))) < 2)) {
    formGather <- paste("depEvents ~", formula, sep = "")
  } else {
    formGather <- formula
  }

  dataProcessed <- GatherPreprocessing(
    as.formula(formGather),
    model = "DyNAM", subModel = submodel,
    progress = FALSE,
    envir = envir
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

  nEvents <- expandedDF$event[length(expandedDF$event)]
  indexEvents <- data.frame(
    "start" = seq(1, length(expandedDF$event), nAct - 1),
    "end" = seq(nAct - 1, length(expandedDF$event), nAct - 1),
    "event" = seq(1, nEvents)
  )
  creation_events <- 0
  deletion_events <- 0
  for (i in 1:nEvents) {
    # Select rows of expandedDF corresponding to the event
    # Check if inertia == 0 (creation of ties) or inertia == 1 (delition of ties) for (receiver == TRUE)
    df_aux <- expandedDF[seq(indexEvents$start[i], indexEvents$end[i]), ]
    index_aux <- df_aux[, "selected"] == TRUE
    if (df_aux[index_aux, "inertia_netEvents"] == 0) { # creation of ties
      creation_events <- c(creation_events, i)
    } else {
      deletion_events <- c(deletion_events, i)
    }
  }
  creation_events <- creation_events[-1]
  deletion_events <- deletion_events[-1]

  expandedDFCreation <- subset(expandedDF, event %in% creation_events & inertia_netEvents == 0)
  expandedDFDeletion <- subset(expandedDF, event %in% deletion_events & inertia_netEvents == 1)

  listExpandedDF <- list(
    expandedDFCreation = expandedDFCreation[, !names(expandedDFCreation) %in% c("inertia_netEvents", "tie_net0", "sender")],
    expandedDFDeletion = expandedDFDeletion[, !names(expandedDFDeletion) %in% c("inertia_netEvents", "tie_net0", "sender")]
  )
  expandedDF <- list(
    expandedDFCreation = subset(expandedDF, event %in% creation_events)[, !names(expandedDFCreation) %in% c("inertia_netEvents", "tie_net0", "sender")],
    expandedDFDeletion = subset(expandedDF, event %in% deletion_events)[, !names(expandedDFCreation) %in% c("inertia_netEvents", "tie_net0", "sender")]
  )

  return(list(listExpandedDF = listExpandedDF, expandedDF = expandedDF))
}
