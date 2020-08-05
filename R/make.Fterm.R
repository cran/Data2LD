make.Fterm <- function(funobj, parvec, estimate, Ufd, factor=1) {
  #  make_Fterm assembles four arguments into a struct object that
  #  defines a single forcing term in a linear differentiable 
  #  equation.
  #  Arguments are:
  #  FUNOBJ   ... Either functional data object that defines the coefficient
  #               and its first derivative, or a user-supplied struct object  
  #               with fields fd and difdip that compute the value of the and 
  #               function its first derivative, respectively, and the field 
  #               more that is used for any additional information required 
  #               by the user-supplied coefficient.
  #  PARVEC   ... A vector of values for either the coefficients for the 
  #               functional data object or for the user-supplied struct
  #               object.
  #  ESTIMATE ... If nonzero, estimate the values in PARVEC within function
  #               Data2LD_Opt.  Otherwise, keep them fixed.
  #               This may also be a logical vector of the same length as
  #               PARVEC specifying which coefficients are fixed and which
  #               are not.  
  #  UFD      ... A functional data object that represents the forcing
  #               variable.
  #  FACTOR   ... The fixed known constant multiplying the coefficient
  #                 function.  Defaults to 1.
 
  #  Last modified 16 April 2020
 
  if (!is.fd(Ufd)) {
    stop("Argument UFD is not a functional data object.")
  }
  
  if (inherits(funobj,c("basisfd","fd","fdPar"))) {
    if (is.basis(funobj))  basisobj <- funobj
    if (is.fd(funobj))     basisobj <- funobj$basis
    if (is.fdPar(funobj))  basisobj <- funobj$fd$basis
    nbasis <- basisobj$nbasis
    if (nbasis != length(parvec)) 
      stop("Number of bases not equal to length of argument PARVEC.")
    nestimate <- length(estimate)
    if (nbasis != nestimate && nestimate > 1)
      stop("Number of bases not equal to length of argumenjt ESTIMATE.")
  }
  
  npar <- length(parvec)
  if (length(estimate) == 1) estimate <- rep(estimate,npar)
  
  index <- NULL
  termList <- list(funobj=funobj, parvec=parvec, index=index, estimate=estimate, 
                   Ufd=Ufd, factor=factor)
  return(termList)
}