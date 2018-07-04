make.coef <- function(funobj, parvec, estimate=TRUE) {
  #  Make a list object that is a coefficient function.
  #  the following text has been deprecated
  #  , coeftype="beta"
  
  #  Last modified 26 June 2018
  
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
  if (!is.numeric(parvec)) {
    stop("Argument PARVEC is not numeric.")
  }
  if (!is.logical(estimate)) {
    stop("Argument ESTIMATE is not logical.")
  }
  # if (!is.character(coeftype)) {
  #   stop("Argument COEFTYPE is not character.")
  # }
  
  coefList <- list(fun=funobj, parvec=parvec, estimate=estimate)
  
  return(coefList)
}