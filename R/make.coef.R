make.coef <- function(funobj, parvec, estimate=TRUE, coeftype="beta") {
  if (!is.function(funobj) && !is.basis(funobj)  && 
      !is.fd(funobj)       && !is.fdPar(funobj)) {
    stop(paste('Argumnent FUNOBJ is not a function handle, basis, fd object or',
               ' a basis object.'))
  }
  if (!is.numeric(parvec)) {
    stop("Argument PARVEC is not numeric.")
  }
  if (!is.logical(estimate)) {
    stop("Argument ESTIMATE is not logical.")
  }
  if (!is.character(coeftype)) {
    stop("Argument COEFTYPE is not character.")
  }
  
  coefList <- list(fun=funobj, parvec=parvec, estimate=estimate, coeftype=coeftype)
  
  return(coefList)
}