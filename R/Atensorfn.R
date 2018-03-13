Atensorfn <- function(modelList, coefList) {
  #  Set up ATENSORLIST as a list of length NVAR defining products of forcing terms.
  #  NFORCEI is the number of forcing terms for this variable.
  #  Each member of AtensorList[[ivar]] is a list of length NFORCEI
  #  and each member of this is a list of length NFORCEI
  
  #  Last modified 2 January 2018
  
  nvar <- length(modelList)
  AtensorList <- vector("list", nvar)
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    if (modelListi$nallFterm > 0) {
      #  there are NFORCE forcing terms for this variable
      nforce <- modelListi$nallFterm
      AtensorListi <- vector("list",nforce)
      for (jx in 1:nforce) {
        AtensorListi[[jx]] <- vector("list",nforce)
      }
      #  loop through forcing terms
      for (jv in 1:nforce) {
        modelListiv <- modelListi$FList[[jv]]
        ncoefv      <- modelListiv$ncoef
        coefListv   <- coefList[[ncoefv]]
        AfdParv     <- coefListv$fun
        if (is.fdPar(AfdParv) || is.fd(AfdParv) || is.basis(AfdParv)) {
          Ufdv      <- modelListiv$Ufd
          if (is.basis(AfdParv)) {
            Abasisv <- AfdParv
          } else {
            if (is.fd(AfdParv)) {
              Abasisv <- AfdParv$basis
            } else {
              Abasisv <- AfdParv$fd$basis
            }
          }
          Atypev    <- Abasisv$type
          nAbasisv  <- Abasisv$nbasis
          Ubasisv   <- Ufdv$basis
          Utypev    <- Ubasisv$type
          nUbasisv  <- Ubasisv$nbasis
          #  loop through forcing terms again
          for (jx in 1:nforce) {
            modelListix <- modelListi$FList[[jx]]
            ncoefx      <- modelListix$ncoef
            coefListx   <- coefList[[ncoefx]]
            AfdParx     <- coefListx$fun
            if (is.fdPar(AfdParx) || is.fd(AfdParx) || is.basis(AfdParx)) {
              Ufdx      <- modelListix$Ufd
              if (is.basis(AfdParx)) {
                Abasisx <- AfdParx
              } else {
                if (is.fd(AfdParx)) {
                  Abasisx <- AfdParx$basis
                } else {
                  Abasisx <- AfdParx$fd$basis
                }
              }
              Atypex    <- Abasisx$type
              nAbasisx  <- Abasisx$nbasis
              Ubasisx   <- Ufdx$basis
              Utypex    <- Ubasisx$type
              nUbasisx  <- Ubasisx$nbasis
              if (Atypev == "const"   && Atypex == "const" &&
                  Utypev == "bspline" && Utypex == "bspline") {
                #  of both coefficients have constant bases, use inprod.Data2LD
                XWXWmatij <- inprod.Data2LD(Ubasisv, Ubasisx, 0, 0)
                XWXWmatij <- matrix(XWXWmatij, nUbasisv*nUbasisx, 1)
              } else {
                # otherwise use inprod.TPbasis
                XWXWmatij <- inprod.TPbasis(Ubasisv, Abasisv, 
                                            Ubasisx, Abasisx, 
                                            0, 0, 0, 0)
              }
              #  as a single column sparse matrix
              AtensorListi[[jx]][[jv]] <- Matrix(XWXWmatij)
            }
          }
        }
      }
      AtensorList[[ivar]] <- AtensorListi
    } else {
      #  there are no forcing terms for this variable
      AtensorList[[ivar]] <- NULL
    }
    
  }
  
  return(AtensorList)
  
}
