make.Fterm <- function(ncoef, Ufd, factor=1) {
  if (floor(ncoef) != ncoef || ncoef < 1) {
    stop("Argument NCOEF is not a positive integer.")
  }
  if (!is.fd(Ufd)) {
    stop("Argument UFD is not a functional data object.")
  }
  termList <- list(ncoef=ncoef, Ufd=Ufd, factor=factor)
  return(termList)
}