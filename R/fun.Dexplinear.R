fun.Dexplinear <- function(tvec, bvec, Bbasisobj) {
  bvec   <- as.matrix(bvec)
  nbasis <- length(bvec)
  if (!is.numeric(tvec))  {
    stop("In fun.Dexplinear argument T is not numeric.")
  }
  if (!is.numeric(bvec))  {
    stop("In fun.Dexplinear argument BVEC is not numeric.")
  }
  if (is.fdPar(Bbasisobj)) {
    Bbasisobj <- Bbasisobj$basis
  }
  if (is.fd(Bbasisobj)) {
    Bbasisobj <- Bbasisobj$basis
  }
  if (!is.basis(Bbasisobj)) {
    stop("In fun.Dexplinear argument BASISOBJ is not a basis object.")
  }
  basismat  <- eval.basis(tvec, Bbasisobj)
  bval      <- exp(basismat %*% bvec)
  Dbval     <- basismat*(bval %*% matrix(1,1,nbasis))
  return(Dbval)
}