# Permutations -------------------------------------

#' Permutation of a sequence of events
#'
#' @param x original data frame with a column of senders, a column of receivers and a column of replace (replace = 0 means deletion of tie, replace = 1 means creation of tie)
#' @param nmax maximum number of permutations performed. The number of permutations computed will be the min{length(x)!, nmax}
#'
#' @return list of dataframes. Each dataframe is a permutation of the rows of the original list x
#' @export
#'
#' @examples library(RSiena)
#' X0 <- s501
#' X1 <- s502
#' seq <- EMPreprocessing(X0, X1)
#' permute(seq, nmax = 5)
permute <- function(x, nmax = 1000) {
  n <- nrow(x)
  if (factorial(n) < nmax) {
    out <- lapply(1:factorial(n), function(t) x[sample(n, n), ])
  } else {
    out <- lapply(1:nmax, function(t) x[sample(n, n), ])
  }
  return(out)
}
