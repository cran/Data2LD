fun.explinear <- function(tvec, bvec, Bbasisobj) {
  bvec   <- as.matrix(bvec)
  if (!is.numeric(t))  {
    stop("In fun.explinear argument T is not numeric.")
  }
  if (!is.numeric(bvec))  {
    stop("In fun.explinear argument BVEC is not numeric.")
  }
  if (is.fdPar(Bbasisobj)) {
    Bbasisobj <- Bbasisobj$basis
  }
  if (is.fd(Bbasisobj)) {
    Bbasisobj <- Bbasisobj$basis
  }
  if (!is.basis(Bbasisobj)) {
    stop("In fun.explinear BASISOBJ is not a basis object.")
  }
  basismat <- eval.basis(tvec,Bbasisobj)
  bval     <- exp(basismat %*% bvec)
  return(bval)
}
