make.Xterm <- function(variable, ncoef, derivative=0, factor=1) {
  if (floor(variable) != variable || variable < 1) {
    stop("Argument VARIABLE is not a positive integer.")
  }
  if (floor(ncoef) != ncoef || ncoef < 1) {
    stop("Argument NCOEF is not a positive integer.")
  }
  if (floor(derivative) != derivative || derivative < 0) {
    stop("Argument DERIVATIVE is not a positive integer.")
  }
  termList <- list(variable=variable, ncoef=ncoef, derivative=derivative, factor=factor)
  return(termList)
}