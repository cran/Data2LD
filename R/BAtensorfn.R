BAtensorfn <- function(XbasisList, modelList) {
  #  Set up BATENSORLIST of length NVAR defining products of
  #  derivative terms and forcing terms in LX.
  #  Each list BATENSORLIST[[ivar]] contains a list array of dimensions
  #        NFORCE[j_1] and NDERIVVEC[l] + 1.
  #  Each list BATENSORLIST[[ivar]][[j_1]][[l_2]] contains the inner product
  #  X basis functions for variable i_2,
  #  variable weight basis functions for variable i_2,
  #  basis functions for forcing function j_1 and variable i_1
  #  forcing weight basis for forcing function j_1 and variable i_1.

  #  Last modified 15 January 2019

  rng     <- XbasisList[[1]]$rangeval
  Wbasism <- create.constant.basis(rng)

  #  set up the structure of BAtensorList

  nvar <- length(modelList)
  BAtensorList <- vector("list", nvar)
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    #if (!is.null(modelListi$XList)) {
    nallXterm  <- modelListi$nallXterm
    nallFterm  <- modelListi$nallFterm
    Xbasisi    <- XbasisList[[ivar]]
    Xtypei     <- Xbasisi$type
    nXbasisi   <- Xbasisi$nbasis
    if (!is.null(modelListi$FList)) {
      nX <- nallXterm + 1
      BAtensorListi <- vector("list",nX)
      for (jx in 1:nX) {
        BAtensorListi[[jx]] <- vector("list",nallFterm)
      }
      BAtensorList[[ivar]] <- BAtensorListi
      order <- modelListi$order
      for (jforce in 1:nallFterm) {
        #  store inner products of mth derivative of X bases and forcing functions
        modelListij <- modelListi$FList[[jforce]]
        AfdParj   <- modelListij$fun
        if (is.fdPar(AfdParj) || is.fd(AfdParj) || is.basis(AfdParj)) {
          if (is.basis(AfdParj)) {
            Abasisj <- AfdParj
          } else {
            if (is.fd(AfdParj)) {
              Abasisj <- AfdParj$basis
            } else {
              Abasisj <- AfdParj$fd$basis
            }
          }
          Atypej    <- Abasisj$type
          Ubasisj   <- modelListij$Ufd$basis
          Utypej    <- Ubasisj$type
          nUbasisj  <- Ubasisj$nbasis
          if (Atypej == "const" && Utypej == "bspline" && Xtypei == "bspline") {
            #  both coeffcients are constants, use inprod.Data2LD
            XWXWmatij <- inprod(Xbasisi, Ubasisj, order, 0)
            #  store matrix by rows as a vector
            XWXWmatij <- matrix(t(XWXWmatij), nXbasisi*nUbasisj, 1)
          } else {
            #  otherwise use inprod.TPbasis
            XWXWmatij <- inprod.TPbasis(Xbasisi, Wbasism, Ubasisj, Abasisj,
                                        order, 0, 0, 0)
          }
          #  store vector as a sparse one-column matrix
          BAtensorList[[ivar]][[nX]][[jforce]] <- XWXWmatij
          #  store inner products of right side derivatives of X bases and
          #  forcing functions
          for (iw in 1:nallXterm) {
            modelListiw <- modelListi$XList[[iw]]
            WfdParw     <- modelListiw$fun
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
              derivative <- modelListiw$derivative
              Xbasisw    <- XbasisList[[modelListiw$variable]]
              Wtypew     <- Wbasisw$type
              Xtypew     <- Xbasisw$type
              nXbasisw   <- Xbasisw$nbasis
              if (Wtypew == "const"   && Atypej == "const" &&
                  Xtypew == "bspline" && Utypej == "bspline") {
                #  both coeffcients are constants, use inprod.Data2LD
                XWXWmatwj <- inprod(Xbasisw, Ubasisj, derivative, 0)
                XWXWmatwj <- matrix(t(XWXWmatwj), nXbasisw*nUbasisj, 1)
              } else {
                XWXWmatwj  <- inprod.TPbasis(Xbasisw, Wbasisw,
                                             Ubasisj, Abasisj,
                                             derivative, 0, 0, 0)
              }
              BAtensorList[[ivar]][[iw]][[jforce]] <- XWXWmatwj
            } else {
              BAtensorList[[ivar]][[iw]][[jforce]] <- NULL
            }
          }
        } else {
          BAtensorList[[ivar]][[nX]][[jforce]] <- NULL
        }
      }
    }
    #}
  }

  return(BAtensorList)

}

