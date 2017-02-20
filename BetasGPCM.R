#' BetasGPC
#'
#' This function simply takes a vector of parameters and converts it to a list.
#' @param thetas I don't know what these are.
#' @param nitems The number of items
#' @param ncatg A vector of length nitems containing the number of categories for each item
#' @param constaint c("gpcm", "rasch", "1PL")
#' @param anchor.param a matrix of the parameters for the anchor set
#' @keywords add one here
#' #@export
#' #@examples
#' BetasGPCM()

BetasGPCM <- function (thetas, nitems, ncatg, constraint, anchor.param = NULL,keep.names = FALSE) {
  betas <- if (constraint == "gpcm") {
    ii <- rep(1:nitems, ncatg)
    split(thetas, ii)
  } else if (constraint == "1PL") {
    nt <- length(thetas)
    ii <- rep(1:nitems, ncatg - 1)
    lapply(split(thetas[-nt], ii), function (x) c(x, thetas[nt]))
  } else {
    ii <- rep(1:nitems, ncatg - 1)
    lapply(split(thetas, ii), function (x) c(x, 1))       
  }
  if (!keep.names)
    names(betas) <- NULL
  if(!is.null(anchor.param)) betas <- c(anchor.param,betas)
  betas
}
