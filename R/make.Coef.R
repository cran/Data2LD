make.Coef <- function(funobj, parvec, estimate=TRUE) {
  # make.Coef assembles the three arguments into a struct object that is used
  # by function modelcheck to set up a linear dynamic equation object.
  #  Arguments are:
  #  FUNOBJ   ... Either a functional data object that defines the coefficient
  #               and its first derivative, or a user-supplied struct object with 
  #               fields fd and difdip that compute the value of the function and 
  #               its first derivative, respectively, and the field more that
  #               is used for any additional information required by the 
  #               user-supplied coefficient.
  #  PARVEC   ... A vector of values for either the coefficients for the 
  #               functional data object or for the user-supplied struct
  #               object.
  #  ESTIMATE ... If nonzero, estimate the values in PARVEC within function
  #               Data2LD.Opt.  Otherwise, keep them fixed.
  
  #  Last modified 27 December 2018
  
  #  check class of funobj
  
  if (!is.list(funobj) && !is.basis(funobj)  && 
      !is.fd(funobj)   && !is.fdPar(funobj)) {
    stop(paste('Argument FUNOBJ is not a function, basis object, fd object or',
               ' a fdPar object.'))
  } else {
    if (!is.basis(funobj) && !is.fd(funobj) && !is.fdPar(funobj)) {
      if (!is.function(funobj$fun)) {
        stop("The first member of argument FUNOBJ is not a function.")
      }
      if (!is.function(funobj$Dfun)) {
        stop("The second member of argument FUNOBJ is not a function.")
      }
    }
  }
  
  #  check class of parvec
  
  if (!is.numeric(parvec)) {
    stop("Argument PARVEC is not numeric.")
  }
  
  #  check dimension of parvec against number of basis functions
  
  nfunbasis = 0;
  if (is.basis(funobj)) {
    nfunbasis <- funobj$nbasis
  }
  if (is.fd(funobj)) {
    nfunbasis = funobj$basis$nbasis
  }
  if (is.fdPar(funobj)) {
    nfunbasis = funobj$fd$basis$nbasis
  }
  
  parvec <- as.matrix(parvec)
  
  parvecdim <- dim(parvec)
  
  if (parvecdim[2] != 1) {
    stop("Parameter vector is not a column vector.")
  }
  
  if (nfunbasis) { 
    if (parvecdim[1] != nfunbasis) {
      stop(paste("Length of parameter vector not equal to", 
                 "the number of basis functions."))
    }
  }
  
  #  check class of estimate
  
  if (!is.logical(estimate)) {
    stop("Argument ESTIMATE is not logical.")
  }
  
  coefList <- list(fun=funobj, parvec=parvec, estimate=estimate)
  
  return(coefList)
  
}