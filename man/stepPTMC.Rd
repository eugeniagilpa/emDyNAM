% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_MCMC.R
\name{stepPTMC}
\alias{stepPTMC}
\title{Parallel tempering step (multi-core)}
\usage{
stepPTMC(
  indexCore,
  splitIndicesPerCore,
  seqs,
  H,
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
  k,
  temp,
  nStepExch,
  pAug,
  pShort,
  pPerm
)
}
\arguments{
\item{seqs}{list of data frames, sequences of events for each temperature.}

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

\item{type}{vector of integers, type of step for each chain.}

\item{T0}{maximum temperature}
}
\value{
list of result of stepMCMC for each sequence with different temperatures.
}
\description{
Performs an step of MCMC given a type:
\item Type 1: augmentation of sequence.
\item Type 2: shortening of sequence.
}
