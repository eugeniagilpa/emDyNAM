#' Step for permuting two elements in a given sequence.
#'
#' @param x vector
#'
#' @return sample of elements of the vector
#'
#' @export
#'
sampleVec <- function(x, ...) x[sample(length(x), ...)]
