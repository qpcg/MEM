#' GPCM Log-likelihood function
#'
#' Calculates the log-likelihood given a set of thetas.
#' @param thetas 
#' @param constraint c("gpcm","rasch","1PL")
#' @param anchor fixed values of the anchor set
#' @keywords add one here
#'
#' #@export
#' #@examples
#' LogLikGPCM()

LogLikGPCM <-  function (thetas, constraint, anchor = NULL){
  betas <- betas.gpcm(thetas[c(-(length(thetas)-1),-length(thetas))], p-nanchor, ncatg[-1:-nanchor], constraint, anchor.param)
  delta <- thetas[(length(thetas)-1)]
  VarOp <- thetas[length(thetas)]
  log.crf <- crf.GPCM(betas, (Z*VarOp + delta), IRT.param, log = TRUE)
  log.p.xz <- matrix(0, nfreqs, length(Z))
  for (j in 1:p) {
    log.pr <- log.crf[[j]]
    xj <- X[, j]
    na.ind <- is.na(xj)
    log.pr <- log.pr[xj, ]
    if (any(na.ind))
      log.pr[na.ind, ] <- 0
    log.p.xz <- log.p.xz + log.pr
  }
  p.x <- rep(exp(log.p.xz) %*% GHw, obs)
  - sum(log(p.x))
}