# Function to create time difference of events ---------

#' Time generator function.
#'
#' @description Generate possible times for a sequence of events following the DyNAM model.
#'
#' @param seq sequence of events
#' @param nAct number of actors
#' @param theta list with Crea and Del objects with rate parameters
#'
#' @return vector of estimated time values for each event in the sequence
#' @export
#'
#' @examples library(RSiena)
#' X0 <- s501
#' X1 <- s502
#' seq <- EMPreprocessing(X0, X1)
#' timeGenerator(seq, ncol(X0), list(c("Crea" = rep(1, ncol(X0)), "Del" = rep(1, ncol(X0)))))
timeGenerator <- function(seq, nAct, theta) {
  X_aux <- cbind("intercept" = rep(1, nAct))

  expXbCrea <- exp(X_aux %*% theta$Crea)
  sumRateCrea <- sum(expXbCrea)
  expXbDel <- exp(X_aux %*% theta$Del)
  sumRateDel <- sum(expXbDel)

  time <- numeric(nrow(seq))
  if (seq$replace[1] == 0) {
    time[1] <- rexp(1, sumRateDel)
  } else {
    time[1] <- rexp(1, sumRateCrea)
  }

  for (i in 2:nrow(seq)) {
    if (seq$replace[i] == 0) {
      time[i] <- time[i - 1] + rexp(1, sumRateDel)
    } else {
      time[i] <- time[i - 1] + rexp(1, sumRateCrea)
    }
  }
  return(time)
}
