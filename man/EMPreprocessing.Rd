% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/function_preprocessing.R
\name{EMPreprocessing}
\alias{EMPreprocessing}
\title{Auxiliary function for EM algorithm}
\usage{
EMPreprocessing(X0, X1)
}
\arguments{
\item{X0}{initial network matrix (row, col names must be actor names)}

\item{X1}{final network matrix (row, col names must be actor names)}
}
\value{
data frame with sequence of events (sender, receiver and replacement)
}
\description{
Auxiliary function for EM algorithm
}
\examples{
library(RSiena)
X0 <- s501
X1 <- s502
EMPreprocessing(X0, X1)
}
