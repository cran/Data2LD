getForceTerm <- function(FtermList, coefList) {
  #  get details for a forcing term
  #  last modified 7 January 2019
  ncoef     <- FtermList$ncoef     # index of coefficient object in coefList
  factor    <- FtermList$factor    # fixed scale factor 
  Ufd       <- FtermList$Ufd       # functional data object for forcing term
  coefnList <- coefList[[ncoef]]   # coefficient object itself
  Avec      <- coefnList$parvec    # parameter vector for forcing coefficient
  Aestim    <- coefnList$estimate  # parameter to be estimated? (TRUE/FALSE)
  Ubasis    <- Ufd$basis           # functional basis object for forcing function
  Ucoef     <- Ufd$coef            # coefficient vector for forcing function
  nUbasis   <- Ubasis$nbasis       # number of basis fucntions
  nAbasis   <- length(Avec)        # number of coefficients for B-spline coeff.
  funtype   <- !(is.basis(coefnList$fun) ||
                 is.fd(   coefnList$fun) ||
                 is.fdPar(coefnList$fun))
  #  return named list
  return(list(coefnList=coefnList, factor=factor, Ufd=Ufd, Avec=Avec, Ubasis=Ubasis, 
              Ucoef=Ucoef, nUbasis=nUbasis, Aestim=Aestim, nAbasis=nAbasis, 
              funtype=funtype))
} 