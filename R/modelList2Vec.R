modelList2Vec <- function(modelList, nparam) {
  #  Extracts the ESTIMATED weight coefficients only from MODELLIST
  #  and assembles them into a vector.
  
  #  Last modified 16 April 2020
  
  if (!is.list(modelList)) {
    stop('MODELlist is not a list object')
  }
  
  nvar     <- length(modelList)
  thetavec <- rep(0,nparam)
  
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    nXList <- length(modelListi$XList)
    if (nXList > 0) {
      for (iw in 1:nXList) {
        modelListiw <- modelListi$XList[[iw]]
        Westimate   <- modelListiw$estimate
        if (any(Westimate)) {
          parveciw <- modelListiw$parvec
          indexiw  <- modelListiw$index
          thetavec[indexiw] <- parveciw
        }
      }
    }
    nFList <- length(modelListi$FList)
    if (nFList > 0) {
      for (jforce in 1:nFList) {
        modelListij <- modelListi$FList[[jforce]]
        Aestimj     <- modelListij$estim
        if (any(Aestimj)) {
          parvecij <- modelListij$parvec
          indexij  <- modelListij$index
          thetavec[indexij] <- parvecij
        }
      }
    }
  }
  
  return(thetavec)
  
}
