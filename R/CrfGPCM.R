#' Cumulative Response Function GPCM
#'
#' This calculates the CRF of the GPCM.
#' @param betas This is a list of betas from the BetasGPCM function
#' @param z I am not sure what z does.
#' @param IRT.param Should IRT parameterization be used? TRUE/FALSE
#' @param log
#' @keywords add one here
#' #@export
#' #@examples
#' CrfGPCM()

CrfGPCM <- function (betas, z, IRT.param = TRUE, log = FALSE, eps = .Machine$double.eps^(1/2)) {
  lapply(linpred.GPCM(betas, z, IRT.param), function (x) {
    num <- exp(apply(x, 2, cumsum))
    if (!is.matrix(num))
      num <- t(num)
    den <- 1 + colSums(num)
    out <- rbind(1/den, num/rep(den, each = nrow(x)))
    if (any(ind <- out == 1))
      out[ind] <- 1 - eps
    if (any(ind <- out == 0))
      out[ind] <- eps
    if (log)
      out <- log(out)
    out
  })
}