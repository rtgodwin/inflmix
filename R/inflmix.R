#' Estimate the one-inflated positive Poisson mixture model (OIPPMM)
#'
#' @param y A vector of positive integers.
#' @param l lambda, a vector of starting values for the positive Poisson
#'   components. If \code{NULL}, starting values will be found via grid search,
#'   and mixture models with successively more components will be estimated
#'   until the non-parametric MLE is found.
#' @param p pi, a vector of starting values for the mixture weights.
#'   \code{l} and \code{p} must be initialized together, or not at all. If
#'   \code{NULL}, grid search and estimation for successive numbers of mixture
#'   components will commence until the non-parametric MLE is found.
#' @param tol Tolerance of the EM algorithm. The EM algorithm proceeds until the
#'   proportional difference between all successive parameter estimates for
#'   lambda and pi are less than \code{tol}. Default is 0.001\%.
#' @param maxLikmethod Maximization method passed to \pkg{maxLik}. Default is
#'   Newton-Raphson.
#' @param maxiters Maximum number of EM iterations.
#' @param minlam The minimum value that a mixture component can have before it
#'   is considered to be a redundant one-inflating component. If any value in
#'   lambda is less than \code{minlam}, the algorithm stops and the
#'   non-parametric MLE is declared to be found. Only relevant if \code{l} and
#'   \code{p} are \code{NULL}, so that \code{inflmix} is searching for the
#'   non-parametric MLE.
#' @param reduntol After the EM algorithm converges, the estimation process will
#'   begin again (including a grid search for new starting values), unless any
#'   two components in lambda are within \code{reduntol} of each other.
#'   The non-parametric MLE is then declared to be found.  Only relevant if
#'   \code{l} and \code{p} are \code{NULL}.
#' @param maxk The maximum number of positive Poisson components to be attempted
#'   in the search for the non-parametric MLE.
#' @return If \code{inflmix} is called with starting values for \code{l} and
#'   \code{p}, returns a list containing:
#'   \tabular{ll}{
#'   \code{termreas} \tab the reason that the EM algorithm terminated (either
#'   convergence or iteration limit) \cr
#'   \code{iterations} \tab the number of iterations until convergence \cr
#'   \code{lambda} \tab the estimated values for the positive Poisson parameters \cr
#'   \code{pi} \tab the estimated values for the component weights \cr
#'   \code{logl} \tab the value of the log-likelihood function evaluated at the
#'   parameter estimates for lambda and pi \cr
#'   \code{n} \tab the sample size, the length of the vector \code{y} \cr
#'   \code{predicted} \tab the predicted counts obtained by evaluting the
#'   probability mass function of the OIPPMM model at the parameter estimates for
#'   lambda and pi, and for \eqn{y = 1,\dots,max(y)} \cr
#'   \code{chisq} \tab the Pearson chi-square distance statistic obtained by
#'   comparing the actual and predicted counts \cr
#'   \code{HTn0} \tab the Horvitz-Thompson estimator for the number of missing
#'   zeros \cr
#'   }
#'   If \code{inflmix} is called without starting values for \code{l} and
#'   \code{p} (\code{l=NULL} and \code{p=NULL}), then \code{inflmix} returns an
#'   object of class 'inflmixNPMLE', a list containing each of the above objects,
#'   for each estimated OIPPMM model with successively more mixture components,
#'   in the search for the non-parametric MLE. An additional object is also provided:
#'   \code{termreasNPMLE} which documents the reason for the termination of the search
#'   for the NPMLE (either NPMLE found, or \code{maxk} reached).
#' @seealso \code{\link{sumOIPPMM}} for a table of summary statistics, and
#' \code{\link{rinflmix}} for the generation of OIPPMM variates.
#' @examples
#' require(maxLik)
#' # Search for the non-parametric MLE
#' sumOIPPMM(inflmix(1:20))
#' # Provide starting values instead of searching for the NPMLE
#' inflmix(1:20, l=c(1, 4), p=c(.4, .4))
#' @import stats
#' @import utils
#' @export
inflmix <- function(y, l=NULL, p=NULL, tol=.00001, maxLikmethod="nr", maxiters=1e4, minlam=0.01, reduntol=0.05, maxk = 4) {

  if(!is.integer(y) && !is.numeric(y)) {stop("y must be of type integer or numeric")}
  if(is.numeric(y) && !floor(y)==y) {stop("y must contain only integers")}
  if(any(y < 1)) {stop("y must be positive")}
  if(is.matrix(y) && ncol(y) > 1) {stop("y must be one-dimensional")}
  if(length(l) != length(p)) {stop("l and p must have the same dimension")}
  if(!is.null(l) && is.null(p) || is.null(l) && !is.null(p)) {stop("l and p must be initialized together, or not at all")}

  y <- as.integer(y)

  oippmmlogl <- function(y, l, p) {
    pmfpp <- function(y,lk) {
      dpois(y, lk) / (1 - dpois(0, lk))
    }
    bigk <- length(l)
    sum(log(rowSums(sapply(1:bigk, function(j) p[j] * pmfpp(y, l[j]))) + (y == 1) * (1 - sum(p))))
  }

  inflmgrid <- function(y, bigk, nlam = 10, npi = 3) {
    lam <- seq(0.1, (max(y) - 2), length.out = nlam)
    lams <- t(combn(lam, bigk))
    pis <- expand.grid(replicate(bigk + 1, 1:npi, simplify=F))
    pis <- pis / rowSums(pis)
    pis <- as.matrix(pis[1:(nrow(pis) - 1), 1:bigk])
    loglmat <- sapply(1:nrow(lams), function(q) sapply(1:nrow(pis), function(r) oippmmlogl(y, lams[q, ], pis[r, ])))
    coords <- which(loglmat == max(loglmat), arr.ind = T)
    list(l = lams[coords[1, 2],], p = pis[coords[1, 1],])
  }

  # PMF of PP distribution
  pmfpp <- function(y,lk) {
    dpois(y, lk) / (1 - dpois(0, lk))
  }

  estimate <- function(l, p) {

    # Not the real log-l. Omits terms not relevant for optimization
    logl <- function(l) {
      sum(sapply(1:bigk, function(j) {sum(w[, j] * (y * log(l[j]) - log(exp(l[j]) - 1)))}))
    }

    # Calculate the weights based on the current values of the lambdas and "pi"s
    getweights <- function(p, l) {
      denom <- rowSums(sapply(1:bigk, function(j) {p[j] * pmfpp(y, l[j])})) + (1 - sum(p)) * (y == 1)
      sapply(1:bigk, function(j) {p[j] * pmfpp(y, l[j]) / denom})
    }

    z <- list()
    iters <- 0L

    repeat {
      iters <- iters + 1L
      if(iters > maxiters) {
        z$termreas <- "Iteration limit reached"
        break
      }

      w <- getweights(p, l)

      # Maximize the log-likelihood (equation (6)), and obtain vector of lambda_hats
      lhat <- maxLik::maxLik(logl, method=maxLikmethod, start=l)$estimate

      # Update the "pi"s
      phat <- colMeans(getweights(p, lhat))

      # Check for convergence
      if(all(abs(c((phat - p) / p, (lhat - l)/ l)) < tol)) {
        z$termreas <- "Convergence reached, within tolerance"
        z$iterations <- iters
        break
      }

      # Update all parameter estimates and continue
      l <- lhat
      p <- phat
    }
    z$lambda <- lhat
    z$pi <- phat
    z$logl <- oippmmlogl(y, lhat, phat)
    z$n <- length(y)
    z$predicted <- z$n * sapply(1:max(y), function(i) {sum(sapply(1:bigk, function(j) {p[j] * pmfpp(i, l[j])})) + (1 - sum(p)) * (i == 1)})
    z$chisq <- sum(((tabulate(y) - z$predicted) ^ 2) / z$predicted)
    z$HTn0 <- sum(sapply(1:bigk, function(j) {(p[j] / sum(p)) * (z$n / (1 - exp(-l[j])) - z$n)}))
    z
  }
  if(is.null(l)) {
    bigk <- 2
    zz <- list()
    class(zz) <- "inflmixNPMLE"
    repeat {
      start <- inflmgrid(y, bigk)
      zz[[bigk - 1]] <- estimate(start$l, start$p)
      names(zz)[bigk - 1] <- paste("K =", bigk)
      if(any(abs(combn(zz[[bigk - 1]]$lambda, 2)[1, ] - combn(zz[[bigk - 1]]$lambda, 2)[2, ]) < reduntol) || any(zz[[bigk - 1]]$lambda < minlam)) {
        zz$termreasNPMLE <- paste("NPMLE found: K =", (bigk - 1))
        return(zz)
      }
      bigk <- bigk + 1
      if(bigk > maxk) {
        zz$termreasNPMLE <- "max K reached"
        return(zz)
      }
    }
  } else {
    bigk <- length(l)
    return(estimate(l, p))
  }
}
