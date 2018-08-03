#' Generate variates from the OIPPMM.
#'
#' @param n the desired sample size
#' @param l lambda, a vector of positive Poisson mixture components
#' @param p pi, a vector of weights for the mixture components
#' @return returns a list containing \code{y} (the generated variates) and \code{n0}
#' (the number of 0s from the Poisson process that were discarded)
#' @seealso \code{\link{inflmix}}, \code{\link{sumOIPPMM}}
#' @examples
#' y <- rinflmix(n=50, l=c(1,4), p=c(.4,.4))$y
#' sumOIPPMM(inflmix(y))
#' @import stats
#' @import utils
#' @export
rinflmix <- function(n, l, p) {
  # Warning:- p must sum to 1
  y <- numeric(0)
  n0 <- 0
  K <- length(l)

  while(length(y) < n) {
    # Choose a component
    lambda_k <- sample(l, 1, prob = p[1:K])
    # Generate a variate
    z <- rpois(1, lambda_k)
    # If the generated variate was a 0, keep generating until non-zero
    while(z == 0) {
      # Note that a 0 was generated
      n0 <- n0 + 1
      # Try again
      z <- rpois(1, lambda_k)
    }
    y <- c(y, z)
  }
  # One-inflate the data
  change1s <- runif(n,0,1)
  y[change1s < tail(p, 1)] <- 1
  results <- list(y = y, n0 = n0)
  return(results)
}
