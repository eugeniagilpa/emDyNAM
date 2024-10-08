% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_MCMC.R
\name{stepPT}
\alias{stepPT}
\title{MCMC step reversible jump proposals}
\usage{
stepPT(
  seq,
  type,
  actDfnodesLab,
  actDfnodes,
  tieNames,
  formula,
  net0,
  beta,
  theta,
  fixedparameters,
  initTime,
  endTime,
  k = 5,
  temp = 1,
  pAug,
  pShort,
  pPerm,
  logLikelihoodStats
)
}
\arguments{
\item{seq}{data frame, sequence of events.}

\item{type}{integer, type of step.}

\item{actDfnodesLab}{vector, node labs.}

\item{tieNames}{vector, tie labs (e.g. "12" is a tie from sender 1
to receiver 2).}

\item{formula}{formula, the one used in goldfish model.}

\item{net0}{matrix, initial observed network.}

\item{beta}{list, estimator of choice parameters of the model
(creation and deletion)}

\item{theta}{dataframe, estimator of rate parameters of the model
(creation and deletion) (constant)}

\item{initTime}{initial time}

\item{endTime}{end time}

\item{k}{number of permutations after augmentations/shortening}
}
\value{
list of
\item{newseq}{new sequence.}
\item{loglikSeq}{log-likelihood from original sequence.}
\item{newloglikSeq}{log-likelihood from new sequence.}
\item{step}{depending on the type of step, result of \code{\link[=stepAugment]{stepAugment()}},
\code{\link[=stepShort]{stepShort()}} or \code{\link[=stepPerm]{stepPerm()}}.}
\item{pUndoStep}{probability of undoing the step of MH acceptance rate.}
}
\description{
Performs an step of MCMC given a type:
\item Type 1: augmentation of sequence.
\item Type 2: shortening of sequence.
followed by k permutation steps.
}
