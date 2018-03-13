Data2LD.Rmat <- function(XbasisList, modelList, coefList, rhoVec=rep(0.5,nvar), ntheta, 
                         BtensorList) {
  #  D2LD ... stands for "Data to Linear Dynamics"
  #  D2LD_R computes the penalty matrix R associated with the homogeneous
  #  portion of a linear differential operator L as well as its partial
  #  derivative with respect to parameters defining the homogeneous portion.
  #  This version inputs BtensorList as an argument.
  #  For a single variable whose approximated in terms of an exansion in
  #  terms of a vector \phi of basis functions, R is
  #                 R <- \int [L \phi(t)] [L \phi(t)]' dt.
  #  R is of order K, the number of basis functions, symmetric, and of rank
  #  K - m where m is the order of the largest derivative in the operator.
  #  The eigenanalysis of R can be used to construct an alternative basis
  #  expansion defined in terms an increasing order of complexity of shape.
  #
  #  If multiple variables are involved, then R is a composite matrix
  #  constructed from inner products and cross-products of the basis 
  #  function vectors associate with each variable.  It's order will be
  #  \sum K_i.
  #
  #  This version approximates the integrals in the penalty terms by using 
  #  inprod_basis to compute the cross-product matrices for the  
  #  \beta-coefficient basis functions and the corresponding derivative of 
  #  the x-basis functions,and the cross-product matrices for the 
  #  \alpha-coefficients and the corresponding U functions.  
  #  These are computed upon the first call to D2LD4, and then retained 
  #  for subsequent calls by using the persistent command.  See lines about 
  #  560 to 590 for this code.
  #
  #  This version disassociates coefficient functions from equation 
  #  definitions to allow some coefficients to be used repeatedly and for
  #  both homogeneous and forcing terms.  It requires an extra argument
  #  COEFList that contains the coefficients and the position of their
  #  coefficient vectors in vector THETA.
  #
  #  Arguments:
  #
  #  BASISLIST ... A functional data object or a BASIS object.  If so, the 
  #               smoothing parameter LAMBDA is set to 0.
  #
  #  MODELLIST...  A List aray of length NVAR. Each List contains a 
  #                struct object with members:              
  #                XList ... list of length number of homogeneous terms
  #                          Each List contains a struct object with members:
  #                          WfdPar ... a fdPar object for the coefficient
  #                          variable   ... the index of the variable
  #                          derivative ... the order of its derivative
  #                          npar       ... if coefficient estimated, its location
  #                                         in the composite vector 
  #                          estimate   ... 0, held fixed, otherwise, estimated
  #                FList ... List arrau of length number of forcing terms
  #                          Each List contains a struct object with members:
  #                          AfdPar   ... an fdPar object for the coefficient
  #                          Ufd      ... an fd object for the forcing function
  #                          npar     ... if coefficient estimated, its location
  #                                       in the composite vector 
  #                          estimate ... 0, held fixed, otherwise, estimated
  #                order     ... the highest order of derivative
  #                name      ... a  tag for the variable
  #                nallXterm ... the number of homogeneous terms
  #                nallFterm ... the number of forcing functions
  #  COEFLIST  ... A list array of length NCOEF.  Each list contaions a
  #                a list object with members:
  #               parvec   ... a vector of parameters
  #               estimate ... 0, held fixed, otherwise, estimated
  #               coeftype ... homogeneous or forcing
  #               fun      ... functional basis, fd, or fdPar object,
  #                            or a struct object for a general function
  #                            with fields:
  #                 fd       ... function handle for evaluating function
  #                 Dfd      ... function handle for evaluating
  #                              partial derivative with respect to parameter
  #                 more     ... object providing additional information for
  #                             evaluating coefficient function
  #  RHOVEC  ... A vector of length NVAR containing values in [0,1].  
  #                The data sums of squares are weighted by RHO and 
  #               the roughness penalty by 1-RHO.
  #
  #  BTENSORList 
  
  #  Last modified 2 January 2018
  
  #  ------------------------------------------------------------------------
  #                         Set up analysis
  #  ------------------------------------------------------------------------
  
  #  compute number of variables
  
  nvar <- length(modelList)
  
  #  Set up a vector NCOEFVEC containing number of coefficients used
  #  for the expansion of each variable
  
  ncoefvec <- rep(0,nvar)
  for (ivar in 1:nvar) {
    ncoefvec[ivar] <- XbasisList[[ivar]]$nbasis
  }
  ncoefsum <- sum(ncoefvec)
  ncoefcum <- cumsum(c(0,ncoefvec))
  
  #  get the width of the time domain
  
  Xrange <- XbasisList[[1]]$rangeval
  T      <- Xrange[2] - Xrange[1]
  
  conbasis  <- create.constant.basis(Xrange)
  coefListi <- list(parvec=1, estimate=FALSE, fun=fd(1,conbasis))
  
  #  ------------------------------------------------------------------------
  #                         Compute penalty matrix Rmat(theta)
  #  ------------------------------------------------------------------------
  
  Rmat <- matrix(0,ncoefsum,ncoefsum)
  m2 <- 0
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    weighti    <- modelListi$weight
    m1  <- m2 + 1
    m2  <- m2 + ncoefvec[ivar]
    indi <- m1:m2
    #  compute first term in expansion, crossproduct of D^m with itself
    Xbasisi   <- XbasisList[[ivar]]
    nXbasisi  <- ncoefvec[ivar]
    nXtermi   <- modelListi$nallXterm
    order     <- modelListi$order
    if (is.null(BtensorList)) {
      print("BtensorList is null")
      order  <- modelListi$order
      Rmatii <- matrix(inprod.Data2LD(Xbasisi, Xbasisi, order, order),nXbasisi,nXbasisi)
    } else {
      Btensii <- BtensorList[[ivar]][[nXtermi+1]][[nXtermi+1]]
      if (any(is.na(Btensii))) stop("NAs in Btensii")
      Rmatii  <- matrix(RmatFn(nXbasisi, 1, nXbasisi, 1, 1, 1, Btensii),
                       nXbasisi,nXbasisi)
    }
    if (any(is.na(Rmatii))) stop("NAs in Rmatii")
    Rmatii <- rhoVec[ivar]*Rmatii/T
    Rmat[indi,indi] <- Rmat[indi,indi] + weighti*Rmatii
    #  compute second term in expansion, vector of cross-products with D^m
    #  compute third term in expansion, matrix of cross-products
    for (iw in 1:nXtermi) {
      modelListiw <- modelListi$XList[[iw]]
      derivw    <- modelListiw$derivative
      ivw       <- modelListiw$variable
      indw      <- (ncoefcum[ivw]+1):ncoefcum[ivw+1]
      nXbasisw  <- ncoefvec[ivw]
      Xbasisw   <- XbasisList[[ivw]]
      ncoefw    <- modelListiw$ncoef
      coefListw <- coefList[[ncoefw]]
      Bvecw     <- coefListw$parvec
      factorw   <- modelListiw$factor
      nWbasisw  <- length(Bvecw)
      funw      <- coefListw$fun
      funtypew <- !(is.basis(funw) || is.fd(funw) ||is.fdPar(funw))
      for (ix in 1:nXtermi) {
        modelListix <- modelListi$XList[[ix]]
        derivx    <- modelListix$derivative
        ivx       <- modelListix$variable
        indx      <- (ncoefcum[ivx]+1):ncoefcum[ivx+1]
        nXbasisx  <- ncoefvec[ivx]
        Xbasisx   <- XbasisList[[ivx]]
        ncoefx    <- modelListix$ncoef
        coefListx <- coefList[[ncoefx]]
        Bvecx     <- coefListx$parvec
        factorx   <- modelListix$factor
        nWbasisx  <- length(Bvecx)
        funx      <- coefListx$fun
        funtypex  <- !(is.basis(funx) || is.fd(funx) ||is.fdPar(funx))
        #  determine nature of coefficient function
        if (funtypew || funtypex) {
            #  user-coded case
            Rmatwx <- matrix(inprod.basis.Data2LD(Xbasisw, Xbasisx, coefListw, coefListx, 
                                            derivw, derivx),nXbasisw,nXbasisx)
        } else {
            #  fda object case       
            Btenswx <- BtensorList[[ivar]][[iw]][[ix]]
            # Rmatwx  <- matrix(inprodwx(nXbasisw, nWbasisw, nXbasisx, nWbasisx, 
            #                     Bvecw, Bvecx, Btenswx),nXbasisw,nXbasisw)
            Rmatwx <- RmatFn(nXbasisw, nWbasisw, nXbasisx, nWbasisx, 
                             Bvecw, Bvecx, Btenswx)
        }
        Rmatwx <- factorw*factorx*rhoVec[ivar]*Rmatwx/T
        Rmat[indw,indx] <- Rmat[indw,indx] + weighti*Rmatwx
      }
      if (funtypew) {
        Rmatiw <- matrix(inprod.basis.Data2LD(Xbasisi, Xbasisw, coefListi, coefListw, 
                                      order, derivw),nXbasisi,nXbasisw)
      } else {
        Btensiw <- BtensorList[[ivar]][[nXtermi+1]][[iw]]
        # Rmatiw  <- inprodwx(nXbasisi, 1, nXbasisw, nWbasisw, 1, Bvecw, Btensiw)
        Rmatiw <- RmatFn(nXbasisi, 1, nXbasisw, nWbasisw, 1, Bvecw, Btensiw)
      }
      Rmatiw  <- factorw*rhoVec[ivar]*Rmatiw/T
      Rmat[indi,indw] <- Rmat[indi,indw] - weighti*Rmatiw
      Rmat[indw,indi] <- Rmat[indw,indi] - weighti*t(Rmatiw)
    }
  }
  
  #  ------------------------------------------------------------------------
  #  Compute partial derivatives of R with respect to theta 
  #  in parvec(1:nthetaHL)
  #  ------------------------------------------------------------------------
  
  DRarray <- array(0,c(ncoefsum,ncoefsum,ntheta))
  
  m2 <- 0
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    weighti    <- modelListi$weight
    m1   <- m2 + 1
    m2   <- m2 + ncoefvec[ivar]
    indi <- m1:m2
    nXtermi <- modelListi$nallXterm
    nXbasisi <- ncoefvec[ivar]
    Xbasisi  <- XbasisList[[ivar]]
    nderivi  <- modelListi$order
    #  select only active coefficients requiring estimation
    #  select all active coefficients
    #  loop through active variables within equation ivar
    #  whose coefficient require estimation
    for (iw in 1:nXtermi) {
      modelListiw <- modelListi$XList[[iw]]
      nderivw     <- modelListiw$derivative
      ncoefw      <- modelListiw$ncoef
      coefListw   <- coefList[[ncoefw]]
      Westimw     <- coefListw$estimate
      factorw     <- modelListiw$factor
      if (Westimw) {
        #  define coefficient of estimated variable and
        #  it's derivative index
        ivw      <- modelListiw$variable
        jvw      <- modelListiw$derivative
        indw     <- (ncoefcum[ivw]+1):ncoefcum[ivw+1]
        indthw   <- coefListw$index
        nXbasisw <- ncoefvec[ivw] 
        Xbasisw  <- XbasisList[[ivw]] 
        funw     <- coefListw$fun
        funtypew <- !(is.basis(funw) || is.fd(funw) ||is.fdPar(funw))
        #  loop through all active variables within equation ivar
        for (ix in 1:nXtermi) {
          modelListix <- modelListi$XList[[ix]]
          nderivx     <- modelListix$derivative
          #  define coefficient of active variable and
          #  it's derivative index
          ivx       <- modelListix$variable
          jvx       <- modelListix$derivative
          indx      <- (ncoefcum[ivx]+1):ncoefcum[ivx+1]
          nXbasisx  <- ncoefvec[ivx]
          Xbasisx   <- XbasisList[[ivx]]
          ncoefx    <- modelListix$ncoef
          coefListx <- coefList[[ncoefx]]
          Bvecx     <- coefListx$parvec
          funx      <- coefListx$fun
          funtypex  <- !(is.basis(funx) || is.fd(funx) ||is.fdPar(funx))
          factorx   <- modelListix$factor
          nWbasisx  <- length(Bvecx)
          #  get the tensor vector for this pair of coefficients
          #  and derivatives
          if (funtypew || funtypex) {
            #  user-coded case
            DRarraywx <- array(inprod.Dbasis.Data2LD(Xbasisw, Xbasisx, 
                                                    coefListw, coefListx,
                                                    nderivw, nderivx),
                              c(nXbasisw,nXbasisx,length(indthw)))
          } else {
            #  fda object case
            Btenswx   <- BtensorList[[ivar]][[iw]][[ix]]
            # DRarraywx <- array(inprodDix(nXbasisw, nWbasisw, nXbasisx, nWbasisx, 
            #                             Bvecx, Btenswx),
            #                   c(nXbasisw,nXbasisx,length(indthw)))
            DRarraywx <- DRarrayFn(nXbasisw, nWbasisw, nXbasisx, nWbasisx, 
                                  Bvecx, Btenswx)
            # print(c(iw,ix))
            # print(round(DRarraywx,4))
          }
          #  rescale the inner product
          DRarraywx <- factorw*factorx*rhoVec[ivar]*DRarraywx/T
          #  increment the inner product and its transpose for
          #  the appropriate location in DRarray
          if (ivw == ivx && jvw == jvx) {
            temp <- DRarray[indw,indw,indthw,drop=FALSE]
            temp <- temp + 2*weighti*DRarraywx
            DRarray[indw,indw,indthw] <- temp
          } else {
            temp <- DRarray[indw,indx,indthw,drop=FALSE]
            temp <- temp + weighti*DRarraywx
            DRarray[indw,indx,indthw] <- temp
            temp <- DRarray[indx,indw,indthw,drop=FALSE]
            temp <- temp + weighti*aperm(DRarraywx,c(2,1,3))
            DRarray[indx,indw,indthw] <- temp
          }
        }
      }
      #  partial derivatives wrt Wcoef for cross-products with D^m
      #  here x <- ivar, Wbasisx is the constant basis, and
      #  Bvecx <- 1
      #  get the tensor vector for this pair of coefficients
      #  and derivatives
      if (funtypew) {
        DRarraywi <- array(inprod.Dbasis.Data2LD(Xbasisw,   Xbasisi, 
                                                coefListw, coefListi,
                                                nderivw,   nderivi),
                          c(nXbasisw,nXbasisi,length(indthw)))
      } else {
        Btenswi   <- BtensorList[[ivar]][[iw]][[nXtermi+1]]
        DRarraywi <- DRarrayFn(nXbasisw, nWbasisw, nXbasisi, 1, 1, Btenswi)
      }
      #  rescale the inner product
      DRarraywi <- factorw*rhoVec[ivar]*DRarraywi/T
      # print("dim(DRarraywi) = ")
      # print(dim(DRarraywi))
      temp <- DRarray[indw,indi,indthw,drop=FALSE]
      # print("dim(temp) = ")
      # print(dim(temp))
      temp <- temp - weighti*DRarraywi
      # print("DRarray[indw,indi,indthw):")
      DRarray[indw,indi,indthw] <- temp
      # print("DRarray[indi,indw,indthw):")
      DRarray[indi,indw,indthw] <- aperm(temp,c(2,1,3))
    }
  }
  return(list(Rmat=Rmat, DRarray=DRarray))
}

# #  ------------------------------------------------------------------------
# 
# inprodwx <- function(nXbasisw, nWbasisw, nXbasisx, nWbasisx,
#                      Bvecw, Bvecx, Btenswx) {
#   Rmatwx <- matrix(0,nXbasisw,nXbasisx)
#   ncum   <- cumprod(c(nWbasisx, nXbasisx, nWbasisw, nXbasisw))
#   for (i in 1:nXbasisw) {
#     for (k in 1:nXbasisx) {
#       for (j in 1:nWbasisw) {
#         for (l in 1:nWbasisx) {
#           ijkl <- (i-1)*ncum[3] + (j-1)*ncum[2] + (k-1)*ncum[1] + l
#           Rmatwx[i,k] <- Rmatwx[i,k] + Bvecw[j]*Bvecx[l]*Btenswx[ijkl]
#         }
#       }
#     }
#   }
#   return(Rmatwx)
# }
# 
# #  ------------------------------------------------------------------------
# 
# inprodDwx <- function(nXbasisw, nWbasisw, nXbasisx, nWbasisx, Bvecx, Btenswx) {
#   DRarraywx <- array(0,c(nXbasisw,nXbasisx,nWbasisw))
#   ncum <- cumprod(c(nWbasisx, nXbasisx, nWbasisw, nXbasisw))
#   for (i in 1:nXbasisw) {
#     for (j in 1:nWbasisw) {
#       for (k in 1:nXbasisx) {
#         for (l in 1:nWbasisx) {
#           ijkl <- (i-1)*ncum[3] + (j-1)*ncum[2] + (k-1)*ncum[1] + l
#           DRarraywx[i,k,j] <- DRarraywx[i,k,j] + Bvecx[l]*Btenswx[ijkl]
#         }
#       }
#     }
#   }
#   return(DRarraywx)
# }
# 
# #  ------------------------------------------------------------------------
# 
# inprodDwi <- function(nXbasisi, nXbasisw, nWbasisw, Btenswi) {
#   DRarraywi <- array(0,c(nXbasisw,nXbasisi,nWbasisw))
#   ncum <- cumprod(c(1, nXbasisi, nWbasisw, nXbasisw))
#   #  compute inner product with respect to the ix dimension
#   for (i in 1:nXbasisw) {
#     for (j in 1:nWbasisw) {
#       for (k in 1:nXbasisi) {
#         for (l in 1) {
#           ijkl <- (i-1)*ncum[3] + 
#                  (j-1)*ncum[2] + 
#                  (k-1)*ncum[1] + l
#           DRarraywi[i,k,j] <- DRarraywi[i,k,j] + Btenswi[ijkl]
#         }
#       }
#     }
#   }
#   return(DRarraywi)
# }
  