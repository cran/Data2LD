modelVec2List <- function(modelList, thetavec) {   
  #  Places ESTIMATED weight modelficients only in THETAVEC into pre-existing 
  #    list array modelLIST.
  
  #  Last modified 16 April 2020
  
  if (!is.list(modelList)) {
    stop('modelList is not a list object.')
  } 
  
  nvar   <- length(modelList)
  modelListnew <- vector("list", nvar)
  
  #  homogeneous term portion of parameter vector
  
  m2 <- 0
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    if (modelListi$nallXterm > 0) {
      for (iw in 1:modelListi$nallXterm) {
        modelListiw <- modelListi$XList[[iw]]
        parveciw    <- modelListiw$parvec
        m1 <- m2 + 1
        m2 <- m2 + length(parveciw)
        modelListiw$parvec <- thetavec[m1:m2]
        modelListi$XList[[iw]] <- modelListiw
      }
    }
    modelListnew[[ivar]] <- modelListi
  }
  
  #  forcing term portion of parameter vector
  
  for (ivar in 1:nvar) {
    modelListi <- modelListnew[[ivar]]
    if (modelListi$nallFterm > 0) {
      for (jforce in 1:modelListi$nallFterm) {
        modelListj <- modelListi$FList[[jforce]]
        parvecj      <- modelListj$parvec
        m1 <- m2 + 1
        m2 <- m2 + length(parvecj)
        modelListj$parvec <- thetavec[m1:m2]
        modelListi$FList[[jforce]] <- modelListj        
      }
    }
    modelListnew[[ivar]] <- modelListi
  }
  
  return(modelListnew)
} 