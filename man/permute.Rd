% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_permute.R
\name{permute}
\alias{permute}
\title{Permutation of a sequence of events}
\usage{
permute(x, nmax = 1000)
}
\arguments{
\item{x}{original data frame with a column of senders, a column of receivers and a column of replace (replace = 0 means deletion of tie, replace = 1 means creation of tie)}

\item{nmax}{maximum number of permutations performed. The number of permutations computed will be the min{length(x)!, nmax}}
}
\value{
list of dataframes. Each dataframe is a permutation of the rows of the original list x
}
\description{
Permutation of a sequence of events
}
\examples{
library(RSiena)
X0 <- s501
X1 <- s502
seq <- EMPreprocessing(X0, X1)
permute(seq, nmax = 5)
}
