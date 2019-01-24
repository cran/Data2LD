modelList2Vec <- function(modelList, coefList) {
  #  Extracts the ESTIMATED weight coefficients only from MODELLIST
  #  and assembles them into a vector.
  
  #  Last modified 28 December 2018
  
  if (!is.list(modelList)) {
    stop('MODELlist is not a list object')
  }
  
  if (!is.list(coefList)) {
    stop('COEFlist is not a list object')
  }

  coefcheckList <- coefCheck(coefList, FALSE)
  coefList      <- coefcheckList$coefList
  ntheta        <- coefcheckList$ntheta
  
  nvar <- length(modelList)
  thetavec <- rep(0,ntheta)
  
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    nXList <- length(modelListi$XList)
    if (nXList > 0) {
      for (iw in 1:nXList) {
        modelListiw <- modelListi$XList[[iw]]
        ncoefiw     <- modelListiw$ncoef
        coefListiw  <- coefList[[ncoefiw]]
        Westimate   <- coefListiw$estimate
        if (Westimate) {
          indexw   <- coefListiw$index
          thetaiw  <- coefListiw$parvec
          thetavec[indexw] <- thetaiw
        }
      }
    }
    nFList <- length(modelListi$FList)
    if (nFList > 0) {
      for (jforce in 1:nFList) {
        modelListj <- modelListi$FList[[jforce]]
        ncoefj     <- modelListj$ncoef
        coefListj  <- coefList[[ncoefj]]
        Aestimate  <- coefListj$estimate
        if (Aestimate) {
          indexj <- coefListj$index
          coefij <- coefListj$parvec
          thetavec[indexj] <- coefij
        }
      }
    }
  }
  
  return(thetavec)
  
}
