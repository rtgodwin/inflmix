#' Summarize the results of several estimated OIPPMMs.
#'
#' @param zz An object of class "\code{inflmixNPMLE}", a result of a call to
#'   \code{\link{inflmix}} with \code{l=NULL} and \code{p=NULL}.
#' @return Prints a summary matrix of strings.
#' @seealso \code{\link{inflmix}}
#' @examples
#' require(maxLik)
#' # Search for the non-parametric MLE
#' a <- inflmix(1:20)
#' # Summarize information about the various models fit in the search for the NPMLE
#' sumOIPPMM(a)
#' @export
sumOIPPMM <- function(zz) {
  npmle <- length(zz) - 1
  lambda <- signif(do.call(c, lapply(1:npmle, function(j) {zz[[j]]$lambda})), 5)
  lambda <- do.call(c, lapply(1:npmle, function(j) c(lambda[sum(1:j):(sum(1:(j + 1)) - 1)], "")))
  pi <- do.call(c, lapply(1:npmle, function(j) {zz[[j]]$pi}))
  weight <- round(do.call(c, lapply(1:npmle, function(j) {c(zz[[j]]$pi, (1 - sum(pi[sum(1:j):(sum(1:(j + 1)) - 1)])))})), 4)
  logl <- do.call(c, lapply(1:npmle, function(j) {c(signif(zz[[j]]$logl, 5), rep("", (j + 1)))}))
  chisq <- do.call(c, lapply(1:npmle, function(j) {c(signif(zz[[j]]$chisq, 5), rep("", (j + 1)))}))
  HTn0 <- do.call(c, lapply(1:npmle, function(j) {c(signif(zz[[j]]$HTn0, 5), rep("", (j + 1)))}))
  HTN <- do.call(c, lapply(1:npmle, function(j) {c(signif(zz[[j]]$HTn0 + zz[[j]]$n, 5), rep("", (j + 1)))}))
  zmat <- cbind(lambda, weight, logl, chisq, HTn0, HTN)
  rownames(zmat) <- unlist(sapply(2:(npmle + 1), function(j) c(paste("K =", j), rep("", j))))
  rownames(zmat)[sum(1:npmle) - 1] <- "(NPMLE)"
  noquote(zmat)
}
