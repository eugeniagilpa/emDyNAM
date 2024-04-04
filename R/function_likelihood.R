# Log-likelihood computation function ----------------

#' log-likelihood computation
#'
#' @param listExpandedDF list returned from function GatherPreprocessingDF
#' @param beta list of choice parameters with elements Crea and Del
#'
#' @return loglikelihood (creation + deletion of events models)
#' @export
#'
logLikelihood = function(listExpandedDF, beta){

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

#' log-likelihood computation: Multi-Core
#'
#' @description given parameters, sequences, formula, initial network and multi-core elements, computes log-likelihood of all the sequences.
#'
#' @return vector of loglikelihoods
#' @export
#'
logLikelihoodMC = function(indexCore,permut,beta,splitIndicesPerCore,actDfnodes=actDfnodes,net0=net0,formula=formula){

  indicesCore = splitIndicesPerCore[[indexCore]]
  resLikelihood = vector("list",length(indicesCore))

  for(i in seq_along(indicesCore)){
    seq = permut[[indicesCore[[i]]]]
    envirPrepro <- new.env()
    seq$time <- seq(1,nrow(seq))
    envirPrepro$seqTime <- seq
    envirPrepro$actDfnodes <- actDfnodes
    envirPrepro$net0 <- net0

    local(
      {
        netEvents = defineNetwork(nodes = actDfnodes,matrix=net0) |>
          linkEvents(changeEvents = seqTime, nodes = actDfnodes)

        depEvents = defineDependentEvents(
          seqTime, nodes = actDfnodes, defaultNetwork = netEvents
        )
      },
      envirPrepro
    )
    listExpandedDF = GatherPreprocessingDF(formula, envir = envirPrepro)

    resLikelihood[[i]] = logLikelihood(listExpandedDF,beta)
  }
  return(unlist(resLikelihood))
}


