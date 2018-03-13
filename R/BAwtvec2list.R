BAwtvec2list <-function(thetavec, coefList) {   
  #  Places ESTIMATED weight coefficients only in THETAVEC into pre-existing 
  #    list array COEFLIST.
  
  #  Last modified 2 January 2018
  
  if (!is.list(coefList)) {
    stop('COEFCELL is not a list object.')
  } 
  
  ncoef  <- length(coefList)
  ntheta <- length(thetavec)
  coefListnew <- vector("list", ncoef)
  for (icoef in 1:ncoef) {
    coefListi <- coefList[[icoef]]
    if (coefListi$estimate) {
      index <- coefListi$index
      if (max(index) > ntheta) {
        stop("Coefficient index out of range.")
      }
      if (!is.numeric(thetavec[index])) {
        stop("theta value not numeric")
      } else {
        coefListi$parvec <-  thetavec[index]
      }
    }
    coefListnew[[icoef]] <- coefListi
  }
  return(coefListnew)
} 

