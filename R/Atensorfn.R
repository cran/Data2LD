Atensorfn <- function(modelList) {
  #  Set up ATENSORLIST as a list of length NVAR defining products of forcing terms.
  #  NFORCEI is the number of forcing terms for this variable.
  #  Each member of AtensorList[[ivar]] is a list of length NFORCEI
  #  and each member of this is a list of length NFORCEI

  #  Last modified 15 April 2020

  nvar <- length(modelList)
  AtensorList <- vector("list", nvar)
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    if (!is.null(modelListi$FList)) {
      #  there are NFORCE forcing terms for this variable
      nforce <- modelListi$nallFterm
      AtensorListi <- vector("list",nforce)
      for (jx in 1:nforce) {
        AtensorListi[[jx]] <- vector("list",nforce)
      }
      #  loop through forcing terms
      for (jv in 1:nforce) {
        modelListiv <- modelListi$FList[[jv]]
        AfdParv     <- modelListiv$fun
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
            AfdParx     <- modelListix$fun
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
                XWXWmatij <- inprod(Ubasisv, Ubasisx, 0, 0)
                XWXWmatij <- matrix(XWXWmatij, nUbasisv*nUbasisx, 1)
              } else {
                XWXWmatij <- inprod.TPbasis(Ubasisv, Abasisv,
                                            Ubasisx, Abasisx,
                                            0, 0, 0, 0)
              }
              #  as a single column sparse matrix
              AtensorListi[[jx]][[jv]] <- XWXWmatij
            }
          }
        }
      }
      AtensorList[[ivar]] <- AtensorListi
    } else {
      #  there are no forcing terms for this variable
      AtensorList[[ivar]] <- 0
    }
  }

  return(AtensorList)

}
