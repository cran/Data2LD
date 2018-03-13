BAwtlist2vec <- function(modelList, coefList) {
  #  Extracts the ESTIMATED weight coefficients only from MODELLIST
  #  and assembles them into a vector.
  
  #  Last revised 2 January 2018
  
  if (!is.list(modelList)) {
    stop('MODELlist is not a list object')
  }
  
  if (!is.list(coefList)) {
    stop('COEFlist is not a list object')
  }

  coefcheckList <- coefCheck(coefList)
  coefList      <- coefcheckList$coefList
  ntheta        <- coefcheckList$ntheta
  
  nvar <- length(modelList)
  thetavec <- rep(0,ntheta)
  
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    for (iw in 1:modelListi$nallXterm) {
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
    if (modelListi$nallFterm > 0) {
      for (jforce in 1:modelListi$nallFterm) {
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
