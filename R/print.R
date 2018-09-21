#' Summarize the results of several estimated OIPPMMs.
#'
#' @param x An object of class "\code{inflmixNPMLE}", a result of a call to
#'   \code{\link{inflmix}} with \code{l=NULL} and \code{p=NULL}, and with
#'   \code{K=NULL}.
#' @param ... Further arguments passed to or from other methods.
#' @return Custom print method for objects of class \code{inflmixNPMLE}.
#' @seealso \code{\link{inflmix}}
#' @examples
#' inflmix(1:20)
#' @export

print.inflmixNPMLE <- function(x, ...) {
  cat("\ntermination reason:", x$termreasNPMLE)
  if(x$termreasNPMLE != "max K reached") {r <- x$KNPMLE + 1}
  else {r <- length(x) - 1}
    lambda <- signif(do.call(c, lapply(1:r, function(j) {x[[j]]$lambda})), 5)
    pi <- signif(do.call(c, lapply(1:r, function(j) {x[[j]]$pi})), 5)
    for(j in 1:r) {
      cat("\n--------------------------------------------\n")
      l <- lambda[(((j - 1) * j / 2) + 1):(j * (j + 1) / 2)]
      p <- pi[(((j - 1) * j / 2) + 1):(j * (j + 1) / 2)]
      o <- rep("", j)
      o[1] <- 1 - sum(p)
      print(data.frame(k = (1:j), lambda = l, pi = p), row.names = F)
      cat("\nomega (1-inflation):", o)
      cat("\niters:", x[[j]]$iterations, "  logl:", signif(x[[j]]$logl,6), "  chisq:", signif(x[[j]]$chisq, 5))
      cat("\nHTn0:", x[[j]]$HTn0, "  HTN:", x[[j]]$HTn0 + x[[j]]$n)
    }
}
