% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_likelihood.R
\name{logLikConstantRate}
\alias{logLikConstantRate}
\title{log-likelihood computation with time information}
\usage{
logLikConstantRate(rate, initTime, endTime, nactors, R)
}
\arguments{
\item{rate}{dataframe with constant rate for creation and deletion}

\item{initTime}{initial time of observation}

\item{endTime}{end time of observation}

\item{nactors}{number of actors}

\item{R}{lenth of the sequence}
}
\value{
loglikelihood (creation + deletion of events models)
(log-kappa with constant rate)
}
\description{
TO DO : change value of \eqn{\pi_{i_r}} from 1/n to corresponding one.
}
