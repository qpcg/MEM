#' Linear predictor for the GPCM
#'
#' Takes betas and z to calculate the linear combination of a response category for the GPCM
#' @param betas Item parameters for a particular
#' @param z a particular theta value
#' #@export
#' #@examples
#' LinpredGPCM()

LinpredGPCM <- function (betas, z, IRT.param = TRUE) {
  lapply(betas, function (x) {
    nx <- length(x)
    if (IRT.param)
      t(x[nx] * outer(z, x[-nx], "-"))
    else
      outer(x[-nx], x[nx] * z , "+")
  })
}