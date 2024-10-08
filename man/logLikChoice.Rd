% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_likelihood.R
\name{logLikChoice}
\alias{logLikChoice}
\title{log-likelihood computation of the choice model}
\usage{
logLikChoice(listExpandedDF, beta)
}
\arguments{
\item{listExpandedDF}{list returned from function GatherPreprocessingDF}

\item{beta}{list of choice parameters with elements Crea and Del}
}
\value{
loglikelihood (creation + deletion of events models)
}
\description{
log-likelihood computation of the choice model
}
