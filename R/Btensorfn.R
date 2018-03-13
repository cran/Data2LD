Btensorfn <- function(XbasisList, modelList, coefList) {
  #  Set up BTENSORList defining homogeneous part of LX  
  #  This is a List array of length NVAR
  #    Each List BTENSORList[[i1]] contains a square List array of order 
  #      NALLXTERM+1.
  #      Each List BTENSORList[[i1}[[j1]][[j2]] contains the tensor product
  #      of the coefficient basis functions for basis i in term j1, 
  #         the derivative order j of the X basis functions in term j1,
  #         the derivative order j of the X basis functions in term j2,
  #         the coefficient basis functions for basis i in term j2,
  #         j1, j2 <- 1,...,NALLXTERM+1 where
  #         j1, j2 <- NALLXTERM+1 is for the mth derivative of the X basis
  #         functions.
  
  #  Last modified 2 January 2018
  
  rng     <- XbasisList[[1]]$rangeval
  Wbasism <- create.constant.basis(rng)
  Wfdm    <- fd(1,Wbasism)
  WfdParm <- fdPar(Wfdm, 0, 0, FALSE)
  
  #  set up the structure of BtensorList
  
  nvar <- length(modelList)
  BtensorList <- vector("list", nvar)
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    #  set up a two-dimensiional list object within each member of BtensorList
    nX <- modelListi$nallXterm+1
    BtensorListi <- vector("list",nX)
    for (jx in 1:nX) {
      BtensorListi[[jx]] <- vector("list",nX)
    }
    BtensorList[[ivar]] <- BtensorListi
    #  loop through derivative orders from 1 to nderiv
    nderiv <- nX
    for (iw in 1:nderiv) {        
      if (iw < nderiv) {
        modelListiw <- modelListi$XList[[iw]]
        ncoefw      <- modelListiw$ncoef
        coefListw   <- coefList[[ncoefw]]
        WfdParw     <- coefListw$fun
        Xbasisw     <- XbasisList[[modelListiw$variable]]
        jw          <- modelListiw$derivative
      } else {
          WfdParw <- WfdParm
          Xbasisw <- XbasisList[[ivar]]
          jw <- modelListi$order
      }
      if (is.fdPar(WfdParw) || is.fd(WfdParw) || is.basis(WfdParw)) {
        if (is.basis(WfdParw)) {
          Wbasisw <- WfdParw
        } else {
          if (is.fd(WfdParw)) {
            Wbasisw <- WfdParw$basis
          } else {
            Wbasisw <- WfdParw$fd$basis
          }
        }
        Wtypew   <- Wbasisw$type
        Xtypew   <- Xbasisw$type
        nXbasisw <- Xbasisw$nbasis
        #  loop through derivative orders from 1 to nderiv inside outer loop
        for (ix in 1:nderiv) {
          if (ix < nderiv) {
            modelListix <- modelListi$XList[[ix]]
            ncoefx    <- modelListix$ncoef
            coefListx <- coefList[[ncoefx]]
            WfdParx   <- coefListx$fun
            Xbasisx   <- XbasisList[[modelListix$variable]]
            jx        <- modelListix$derivative
          } else {
            WfdParx <- WfdParm
            Xbasisx <- XbasisList[[ivar]]
            jx      <- modelListi$order
          }
          if (is.fdPar(WfdParx) || is.fd(WfdParx) || is.basis(WfdParx)) {
            if (is.basis(WfdParx)) {
              Wbasisx <- WfdParx
            } else {
              if (is.fd(WfdParx)) {
                Wbasisx <- WfdParx$basis
              } else {
                Wbasisx <- WfdParx$fd$basis
              }
            }
            Wtypex   <- Wbasisx$type
            Xtypex   <- Xbasisx$type
            nXbasisx <- Xbasisx$nbasis
            if (Wtypew == "const"   && Wtypex == "const"   && 
                Xtypew == "bspline" && Xtypex == "bspline" ) {
              # if both coefficients have constant bases, evaluate using inprod.Data2LD
              XWXWmatij <- inprod.Data2LD(Xbasisw, Xbasisx, jw, jx)
              XWXWmatij <- matrix(XWXWmatij, nXbasisw*nXbasisx, 1)
            } else {
              # otherwise evaluate using inprod.TPbasis  
              XWXWmatij <- inprod.TPbasis(Xbasisw, Wbasisw, Xbasisx, Wbasisx, 
                                          jx, 0, jw, 0)
            }
            #  save as a single column sparse matrix
            BtensorList[[ivar]][[ix]][[iw]] <- Matrix(XWXWmatij)
          }
        }
      }
    }
  }
  
  return(BtensorList) 
    
}
