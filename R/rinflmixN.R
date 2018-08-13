#' Generate variates from the OIPPMM, beginning with unknown population size, \code{N}.
#'
#' @param N the unknown population size
#' @param l lambda, a vector of positive Poisson mixture components
#' @param p pi, a vector of weights for the mixture components
#' @return a vector of generated OIPPMM variates
#' @seealso \code{\link{inflmix}}, \code{\link{sumOIPPMM}}, \code{\link{rinflmix}}
#' @seealso \code{\link{inflmix}}, \code{\link{sumOIPPMM}}
#' @examples
#' y <- rinflmixN(N=500, l=c(1,4), p=c(.4,.4))
#' sumOIPPMM(inflmix(y))
#' @import stats
#' @import utils
#' @export
rinflmixN <- function(N, l, p) {
  # Use a transformation similar to one from Boehning and Kuhnert (2006)
  q <- sapply(1:length(p), function(j) {
    p[j] / (1 - dpois(0, l[j])) / (sum(p / (1 - dpois(0, l))) + (1 - sum(p)))
  })
  lmixed <- sample(l, N, replace = T, prob = q)
  y <- rpois(N, lmixed)
  y <- y[y != 0]
  change1s <- runif(length(y), 0, 1)
  y[change1s <= (1 - sum(p))] <- 1
  return(y)
}
