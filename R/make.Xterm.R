make.Xterm <- function(funobj, parvec, estimate, variable, derivative=0, factor=1) {
  #  make_Xterm assembles four arguments into a struct object that
  #  defines a single homogeneous term in a linear differentiable 
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
  #  VARIABLE   ... The index of the variable in the system
  #  DERIVATIVE ... The order of the derivative of the variable in the term.
  #                 Defaults to 0.
  #  FACTOR     ... The fixed known constant multiplying the coefficient
  #                 function.  Defaults to 1.
  
  #  Last modified 11 June 2020
  
  if (floor(variable) != variable || variable < 1) {
    stop("Argument VARIABLE is not a positive integer.")
  }
  if (floor(derivative) != derivative || derivative < 0) {
    stop("Argument DERIVATIVE is not a positive integer.")
  }
  
  if (inherits(funobj, c("basisfd","fd","fdPar"))) {
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
  XtermList <- list(funobj=funobj, parvec=parvec, index=index, estimate=estimate, 
                    variable=variable, derivative=derivative, factor=factor)
  return(XtermList)
}