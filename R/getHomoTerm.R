getHomoTerm <- function(XtermList, coefList) {
#  get details for homogeneous term 
#  last modified 5 January 2019
iv        <- XtermList$variable    # index of the variable in this term
ncoef     <- XtermList$ncoef       # index of coefficient object in coefList
factor    <- XtermList$factor      # fixed scale factor
nderiv    <- XtermList$derivative  # order of derivative in this term
coefnList <- coefList[[ncoef]]     # coefficient object itself
Bvec      <- coefnList$parvec      # parameter vector for forcing coefficient
estim     <- coefnList$estimate    # parameter to be estimated? (TRUE/FALSE)
nWbasis   <- length(Bvec)          # number of basis functions
# is a user-defined coefficient (TRUE/FALSE)
funtype   <- !(is.basis(coefnList$fun) ||
               is.fd(   coefnList$fun) ||
               is.fdPar(coefnList$fun))
#  return named list
return(list(coefnList=coefnList, iv=iv, factor=factor, Bvec=Bvec, 
            nWbasis=nWbasis, estim=estim, nderiv=nderiv, funtype=funtype))
}
