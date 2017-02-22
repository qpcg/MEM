#' MEM Estimation for the GPCM
#'
#' This function will allow you to calibrate and link IRT models for the GPCM.
#' @param data A matrix of response
#' @param anchor A vector indicating which items are anchor items
#' @param anchor.param The previous item estimates for the ancor items
#' @param constraint c("gpcm","1PL","rasch")
#' @param hess A logical indicating whether or not to calculate and return the hessian (for standard errors)
#' @param IRT.param Logical indicating whether to use IRT parameterization
#' @param trace Interger specifcying progress of optimization routine
#' @keywords add one here
#' @export
#' #@examples
#' GpcmMEM()

GpcmMEM <- function (data, anchor, anchor.param = NULL,constraint = c("gpcm","1PL","rasch"),hess = F, IRT.param = TRUE, start.val = NULL, na.action = NULL, trace = T,control = list()) {
  
  # Matches the call to allow for fuzzy args specification.
  cl <- match.call()
  
  # Throws an error if not supplied data frame, matrix, or only on column.
  if ((!is.data.frame(data) & !is.matrix(data)) || ncol(data) == 1)
    stop("'data' must be either a numeric matrix or a data.frame, with at least two columns.\n")
  
  # Allows for fuzzy specification of constraint
  constraint <- match.arg(constraint)
  
  # Converts to a data frame.
  X <- if (!is.data.frame(data)) as.data.frame(data) else data
  X <- data.frame(X)
  p <- ncol(X)
  colnamsX <- colnames(X)
  
  # reorders to put anchor first
  if(!is.null(anchor)){
    print(1)
    nanchor <- length(anchor)
    X <- X[,c(anchor, c(1:p)[-anchor])]
  }
  
  # Converts to a factor with levels.
  X[] <- lapply(X, factor)
  
  # Determines the number of categories fo each
  ncatg <- as.vector(sapply(X, function (x) length(levels(x))))
  
  # Convert back to numeric
  X <- sapply(X, unclass)
  
  # Do whatever we are supposed to with the missing data
  #if (!is.null(na.action))
  #  X <- na.action(X)
  
  # Name the columns using original names.
  dimnames(X) <- NULL
  p <- ncol(X)
  
  # Get a list of the response patterns
  pats <- apply(X, 1, paste, collapse = "/")
  
  # Find the frequencies of each of those patterns
  freqs <- table(pats)
  
  # Find the number of unique freqiencies
  nfreqs <- length(freqs)
  
  
  obs <- as.vector(freqs)
  X <- unlist(strsplit(cbind(names(freqs)), "/"))
  X[X == "NA"] <- as.character(NA)
  
  X <- matrix(as.numeric(X), nfreqs, p, TRUE)
  
  # Convert to a pattern score
  XX <- lapply(1:p, function (j) outer(X[, j], seq(1, ncatg[j] - 1), ">") * 1)
  
  # Supply the optimization constraints
  con <- list(iter.qN = 150, GHk = 50, optimizer = "nlminb", optimMethod = "BFGS", numrDeriv = "fd",
              epsHes = 1e-06, parscale = NULL, verbose = getOption("verbose"))
  
  
  #namC <- names(con)
  #con[(namc <- names(control))] <- control
  #if (length(noNms <- namc[!namc %in% namC]) > 0) 
  #  warning("unknown names in control: ", paste(noNms, collapse = ", "))
  
  # Get the Gaussian-Hermite Quadratature points and weights
  GH <- GHpoints(data ~ z1, con$GHk)
  Z <- GH$x[, 2]
  GHw <- GH$w
  
  # Get the starting values for the optimizer
  init.thetas <- start.val.gpcm(NULL, X[,-1:-nanchor], obs, constraint, ncatg[-1:-nanchor], IRT.param =T)
  init.thetas <- c(init.thetas,0,1)
  # Put the environment for with gpcm for scoregpcm and loglikgpcm
  environment(loglikgpcm) <- environment(scoregpcm) <- environment()
  
  # Set the optimization arugments, run and return.
  
  if (is.null(con$parscale) || length(con$parscale) != length(init.thetas))
    con$parscale <- rep(0.1, length(init.thetas))
  
  res.qN <- optim(init.thetas, loglikgpcm,hessian = hess,constraint = constraint,method = con$optimMethod,
                  control = list(maxit = con$iter.qN, parscale = con$parscale, trace = trace))
  
  if(hess){
    Hess <- res.qN$hessian 
  } else{
    Hess <- NA
  }
  
  if(res.qN$convergence != 0) {
    if(!is.null(res.qN$message))
      warning("Not successful convergence: ", res.qN$message, ".\n")
    else
      warning("Not successful convergence.\n")
  }
  if(hess){
    if (all(!is.na(res.qN$hessian) & is.finite(res.qN$hessian))) {
      ev <- eigen(res.qN$hessian, TRUE, TRUE)$values
      if (!all(ev >= -1e-06 * abs(ev[1]))) 
        warning("Hessian matrix at convergence is not positive definite; unstable solution.\n")
    } else 
      warning("Hessian matrix at convergence contains infinite or missing values; unstable solution.\n")
  }
  
  thetas <- betas.gpcm(res.qN$par[c(-length(res.qN$par),-(length(res.qN$par)-1))], p-nanchor, ncatg[-1:-nanchor], constraint)
  names(thetas) <- if (!is.null(colnamsX)) colnamsX[-anchor] else paste("Item", 1:(p-nanchor))
  thetas <- lapply(thetas, function (x) { names(x) <- c(paste("Catgr.", seq(1, length(x) - 1), sep = ""), "Dscrmn"); x })
  #max.sc <- max(abs(scoregpcm(res.qN$par, constraint)), na.rm = TRUE)
  gamma <- res.qN$par[length(res.qN$par)]
  delta <- res.qN$par[(length(res.qN$par)-1)]
  fit <- list(coefficients = thetas, delta = delta, gamma = gamma, log.Lik = -res.qN$value, convergence = res.qN$conv, hessian = res.qN$hessian, 
              counts = res.qN$counts, patterns = list(X = X, obs = obs), GH = list(Z = Z, GHw = GHw), 
              constraint = constraint, IRT.param = IRT.param, X = data, control = con, na.action = na.action,
              call = cl)
  class(fit) <- "gpcm"
  fit
}
