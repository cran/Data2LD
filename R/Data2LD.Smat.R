Data2LD.Smat <- function(XbasisList, modelList, coefList, rhoVec=rep(0.5,nvar), ntheta, 
                         BAtensorList, nrep, nforce) {
  #  D2LD ... stands for "Data to Linear Dynamics"
  #  D2LD_S computes the penalty matrix S associated with the forcing
  #  portion of a linear differential operator L, as well as its partial
  #  derivatives with  respect to the parameter vector.
  #  For a single variable whose approximated in terms of an exansion in
  #  terms of a vector \phi of basis functions, S is
  #                 S <- \int [L \phi(t)] U' dt.
  #  S has dimensions K and NREP, where K is the number of basis
  #  functions in the expansion of the variable, NFORCE is the number of
  #  forcing functions, and NREP is the number of replications.  The
  #  forcing functions are assumed to vary from one replication to another.
  #  This version loads BAtensorList as an argument.
  #
  #  If multiple variables are involved, then S is a composite matrix
  #  constructed from inner products and cross-products of the basis 
  #  function vectors associate with each variable.  It's dimension will be
  #  \sum K_i by NFORCE*NREP.
  #
  #  This version approximates the integrals in the penalty terms by using 
  #  inprod_basis to compute the cross-product matrices for the  
  #  \beta-coefficient basis functions and the corresponding derivative of 
  #  the x-basis functions,and the cross-product matrices for the 
  #  \alpha-coefficients and the corresponding U functions.  
  #  These are computed upon the first call to D2LD, and then retained 
  #  for subsequent calls by using the persistent command.  See lines about 
  #  560 to 590 for this code.
  #
  #  This version disassociates coefficient functions from equation 
  #  definitions to allow some coefficients to be used repeatedly and for
  #  both homogeneous and forcing terms.  It requires an extra argument
  #  COEFLIST that contains the coefficients and the position of their
  #  coefficient vectors in vector THETA.
  #
  #  Arguments:
  #
  #  BASISLIST ... A functional data object or a BASIS object.  If so, the 
  #               smoothing parameter LAMBDA is set to 0.
  #  MODELLIST  ... A list of length NVAR. Each list contains a 
  #                struct object with members:              
  #                XList ... list List of length number of homogeneous terms
  #                          Each list contains a struct object with members:
  #                          WfdPar ... a fdPar object for the coefficient
  #                          variable   ... the index of the variable
  #                          derivative ... the order of its derivative
  #                          npar ... if coefficient estimated, its location
  #                                   in the composite vector 
  #                FList ... list array of length number of forcing terms
  #                          Each list contains a struct object with members:
  #                          AfdPar ... an fdPar object for the coefficient
  #                          Ufd    ... an fd object for the forcing function
  #                          npar ... if coefficient estimated, its location
  #                                   in the composite vector 
  #                order     ... the highest order of derivative
  #                name      ... a  tag for the variable
  #                nallXterm ... the number of homogeneous terms
  #                nallFterm ... the number of forcing functions
  #  COEFLIST  ... A list of length NCOEF.  Each list contaions a
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
  #  RHOVEC      A value in [0,1].  The data are weighted by P and the
  #               roughness penalty by 1-P.
  #  BAtensorList ... A list of four-way tensors required for Smat
  #  NREP       ... The number of replications of the system.  
  #
  #  Last modified 2 January 2018
  
  #  ------------------------------------------------------------------------
  #                         Set up analysis
  #  ------------------------------------------------------------------------
    
  #  compute number of variables
  
  nvar <- length(modelList)
  
  #  Set up a vector NCOEFVEC containing number of coefficients used
  #  for the expansion of each variable
  
  ncoefvec <- matrix(0,nvar,1);
  for (ivar in 1:nvar) {
    ncoefvec[ivar] <- XbasisList[[ivar]]$nbasis
  }
  ncoefcum <- cumsum(c(0,ncoefvec))
  
  #  get the width of the time domain
  
  Xrange <- XbasisList[[1]]$rangeval
  T      <- Xrange[2] - Xrange[1]
  
  #--------------------------------------------------------------------------
  #                 Compute the penalty vector S(theta)
  #--------------------------------------------------------------------------
  
  ncoefsum <- sum(ncoefvec)
  if (sum(nforce)==0) {
    Smat <- NULL
    DSarray <- NULL
    return(list(Smat <- Smat, DSarray <- DSarray))
  }

  Smat <- matrix(0,sum(ncoefvec),nrep)
  m2 <- 0
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    weighti <- modelListi$weight
    m1  <- m2 + 1
    m2  <- m2 + ncoefvec[ivar]
    ind <- m1:m2
    Xbasisi  <- XbasisList[[ivar]]
    order   <- modelListi$order
    if (nforce[ivar] > 0) {
      nXbasisi <- ncoefvec[ivar]
      nXtermi <- modelListi$nallXterm
      nFtermi <- modelListi$nallFterm
      for (jforce in 1:nFtermi) {
        modelListij <- modelListi$FList[[jforce]]
        ncoefj    <- modelListij$ncoef
        coefListj <- coefList[[ncoefj]]
        Avecj     <- coefListj$parvec
        factorj   <- modelListij$factor
        Ufdj      <- modelListij$Ufd
        Ubasisj   <- Ufdj$basis
        Ucoefj    <- Ufdj$coef
        nUbasisj  <- Ubasisj$nbasis
        nAbasisj  <- length(Avecj)
        funj      <- coefListj$fun
        funtypej <- !(is.basis(funj) || is.fd(funj) ||is.fdPar(funj))
        #  Crossproducts of homogeneous terms with forcing terms
        for (iw in 1:nXtermi) {
          modelListiw <- modelListi$XList[[iw]]
          ivw         <- modelListiw$variable
          indw        <- (ncoefcum[ivw]+1):ncoefcum[ivw+1]
          nXbasisw    <- ncoefvec[ivw]
          Xbasisw     <- XbasisList[[ivw]]
          ncoefw      <- modelListiw$ncoef
          coefListw   <- coefList[[ncoefw]]
          WfdParw     <- coefListw$fun
          Bvecw       <- coefListw$parvec
          factorw     <- modelListiw$factor
          nWbasisw    <- length(Bvecw)
          derivw      <- modelListiw$derivative
          funw        <- coefListw$fun
          funtypew <- !(is.basis(funw) || is.fd(funw) ||is.fdPar(funw))
          if (funtypej || funtypew) {
            Smatjw <- matrix(inprod.basis.Data2LD(Xbasisw, Ufdj, coefListw, coefListj,
                                                  derivw, 0),nXbasisw,nrep)
          } else {
            BAtenswj <- BAtensorList[[ivar]][[iw]][[jforce]]
            # Smatjw <- matrix(inprodijw(nXbasisw, nWbasisw, nUbasisj, nAbasisj, nrep, 
            #                           Bvecw, Avecj, Ucoefj, BAtenswj),nXbasisw,nrep)
            Smatjw <- SmatFn(nXbasisw, nWbasisw, nUbasisj, nAbasisj, nrep, 
                            Bvecw, Avecj, Ucoefj, BAtenswj)
          }
          Smatjw <- factorw*factorj*rhoVec[ivar]*Smatjw/T
          temp <- Smat[indw,,drop=FALSE]
          temp <- temp + weighti*Smatjw
          Smat[indw,] <- temp
        }
        #  Crossproducts of D^m with forcing terms
        if (funtypej) {
          Smatji <- matrix(inprod.basis.Data2LD(Xbasisi, Ufdj, 1, coefListj,order, 0),
                          nXbasisi,nrep)
        } else {
          BAtenswj <- BAtensorList[[ivar]][[nXtermi+1]][[jforce]]
          # Smatji <- matrix(inprodijw(nXbasisi, 1, 
          #                           nUbasisj, nAbasisj, nrep, 
          #                           1, Avecj, Ucoefj, BAtenswj),nXbasisi,nrep)
          Smatji <- SmatFn(nXbasisi, 1, nUbasisj, nAbasisj, nrep, 
                          1, Avecj, Ucoefj, BAtenswj)
        }
        Smatji <- factorj*rhoVec[ivar]*Smatji/T
        temp <- Smat[indw,]
        temp <- temp - weighti*Smatji
        Smat[indw,] <- temp
      }
    }
  }
  
  #  ------------------------------------------------------------------------
  #  Compute partial derivatives of Smat if required with respect to theta
  #  in parvec(1:nthetaFL)
  #  ------------------------------------------------------------------------
  DSarray <- rep(0,ncoefsum*nrep*ntheta)
  m2 <- 0
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    weighti <- modelListi$weight
    m1   <- m2 + 1
    m2   <- m2 + ncoefvec[ivar]
    indi <- m1:m2
    if (nforce[ivar] > 0) {
      modelListi <- modelList[[ivar]]
      nXbasisi <- ncoefvec[ivar]
      nXtermi <- modelListi$nallXterm
      nFtermi <- modelListi$nallFterm
      Xbasisi <- XbasisList[[ivar]]
      order   <- modelListi$order
      #  partial derivatives of product of homogeneous terms
      #  and forcing terms with respect to homogeneous coefficients
      #  loop through all active forcing terms
      for (jforce in 1:nFtermi) {
        modelListij <- modelListi$FList[[jforce]]
        ncoefj    <- modelListij$ncoef
        coefListj <- coefList[[ncoefj]]
        Avecj     <- coefListj$parvec
        Aestimj   <- coefListj$estimate
        factorj   <- modelListij$factor
        nAbasisj  <- length(Avecj)
        Ufdj      <- modelListij$Ufd
        Ucoefj    <- Ufdj$coef
        Ubasisj   <- Ufdj$basis
        nUbasisj  <- Ubasisj$nbasis
        funj      <- coefListj$fun
        funtypej <- !(is.basis(funj) || is.fd(funj) ||is.fdPar(funj))
        #  crossproducts of homogeneous terms with forcing terms
        for (iw in 1:nXtermi) {
          modelListiw <- modelListi$XList[[iw]]
          ncoefw      <- modelListiw$ncoef
          coefListw   <- coefList[[ncoefw]]
          Westimw     <- coefListw$estimate
          nderivw     <- modelListiw$derivative
          factorw     <- modelListiw$factor
          funw        <- coefListw$fun
          funtypew <- !(is.basis(funw) || is.fd(funw) ||is.fdPar(funw))
          factorw     <- modelListiw$factor
          if (Westimw) {
            ivw      <- modelListiw$variable
            indw     <- (ncoefcum[ivw]+1):ncoefcum[ivw+1]
            indthw   <- coefListw$index
            nWbasisw <- length(indthw)
            nXbasisw <- ncoefvec[ivw]
            if (funtypej || funtypew) {
              DBSarrayjw <- array(inprod.Dbasis.Data2LD(Xbasisw, Ufdj, 
                                        coefListw, coefListj, nderivw, 0),
                                        c(nXbasisw,nrep,length(indthw)))
            } else {
              BAtenswj  <- BAtensorList[[ivar]][[iw]][[jforce]]
              DBSarrayjw <- DBSarrayFn(nXbasisw, nWbasisw, nUbasisj, nAbasisj, 
                                      nrep, Avecj, Ucoefj, BAtenswj)
            }
            # rescale DSarrayjw
            DBSarrayjw <- weighti*factorj*factorw*rhoVec[ivar]*DBSarrayjw/T
            for (irep in 1:nrep){
              for (l in 1:nWbasisw) {
                offset.out <- nXbasisw*(irep-1 + (l-1)*nrep)
                offset.in  <- nXbasisw*(irep-1 + 
                             (l + indthw[1] - 2)*nrep)
                DSarray[indw+offset.in] <- DSarray[indw+offset.in] +
                  DBSarrayjw[(1:nXbasisw) + offset.out]
              }
            }
            # temp <- DSarray[indw,,indthw,drop=FALSE]
            # temp <- temp + DBSarrayjw
            # DSarray[indw,,indthw] <- temp
          }
        }
        #  partial derivatives wrt forcing term coefficients for
        #  those forcing terms reqXuiring estimation of their
        #  coefficients
        #  loop through all forcing terms with coefficients to be
        #  estimated
        if (Aestimj) {
          indtha   <- coefListj$index
          #  partial derivatives of products of homogeneous terms
          #  with forcing terms wrt forcing coefficients
          for (iw in 1:nXtermi) {
            modelListiw <- modelListi$XList[[iw]]
            ivw       <- modelListiw$variable
            indw      <- (ncoefcum[ivw]+1):ncoefcum[ivw+1]
            nXbasisw  <- ncoefvec[ivw]
            ncoefw    <- modelListiw$ncoef
            coefListw <- coefList[[ncoefw]]
            Bvecw     <- coefListw$parvec
            factorw   <- modelListiw$factor
            nderivw   <- modelListiw$derivative
            nWbasisw  <- length(Bvecw)
            funw      <- coefListw$fun
            funtypew  <- !(is.basis(funw) || is.fd(funw) ||is.fdPar(funw))
            if (funtypew || funtypej) {
              DASarrayjw <- array(inprod.Dbasis.Data2LD(Ufdj, Xbasisw, 
                                          coefListj, coefListw,  nderivw, 0),
                                          c(nXbasisw,nrep,nAbasisj))
            } else {
              BAtenswj   <- BAtensorList[[ivar]][[iw]][[jforce]]
              DASarrayjw <- DASarrayFn(nXbasisw, nWbasisw, nUbasisj, nAbasisj, 
                                      nrep, Bvecw, Ucoefj, BAtenswj)
            }
            #  rescale DSarrayjw
            DASarrayjw <- weighti*factorj*factorw*rhoVec[ivar]*DASarrayjw/T 
            for (irep in 1:nrep){
              for (l in 1:nAbasisj) {
                offset.out <- nXbasisw*(irep-1 + (l-1)*nrep)
                offset.in  <- nXbasisw*(irep-1 + 
                                         (l + indtha[1] - 2)*nrep)
                DSarray[indw + offset.in] <- DSarray[indw + offset.in] +
                                  DASarrayjw[(1:nXbasisw) + offset.out]
              }
            }
          }
          #  partial derivatives of cross-products of D^m
          #  with forcing terms wrt forcing coefficients
          if (funtypej) {
            DASarrayji <- array(inprod.Dbasis.Data2LD(Ufdj, Xbasisi, 
                                      coefListj, 1, 0, order),
                                      c(nXbasisi,nrep,nAbasisj))
          } else {
            BAtensij   <- BAtensorList[[ivar]][[nXtermi+1]][[jforce]]
            DASarrayji <- DASarrayFn(nXbasisi, 1, nUbasisj, nAbasisj, 
                                    nrep, 1, Ucoefj, BAtensij)
          }
          #  rescale DSarrayji
          DASarrayji <- weighti*factorj*rhoVec[ivar]*DASarrayji/T 
          for (irep in 1:nrep){
            for (l in 1:nAbasisj) {
              offset.out <- nXbasisi*(irep-1 + (l-1)*nrep)
              offset.in  <- nXbasisi*(irep-1 + 
                                       (l + indtha[1] - 2)*nrep)
              DSarray[indi + offset.in] <- DSarray[indi + offset.in] -
                DASarrayji[(1:nXbasisi) + offset.out]
            }
          }
        }
      }
    }
  }
  
  return(list(Smat=Smat, DSarray=DSarray))
  
}

# #  ------------------------------------------------------------------------
# 
# inprodDjw <- function(nXbasisw, nWbasisw, nUbasisj, nAbasisj, 
#                      nrep, Bvecw, Ucoefj, BAtenswj) {
#   ncum <- cumprod(c(nAbasisj, nUbasisj, nWbasisw, nXbasisw))
#   DASarray <- array(0,c(nXbasisw,nrep,nAbasisj))
#   for (irep in 1:nrep) {
#     for (i in 1:nXbasisw) {
#       for (j in 1:nWbasisw) {
#         for (k in 1:nUbasisj) {
#           for (l in 1:nWbasisw) {
#             ijkl <- (i-1)*ncum[3] + (j-1)*ncum[2] + (k-1)*ncum[1] + l
#             DASarray[i,irep,l] <- DASarray[i,irep,l] + 
#               Bvecw[j]*BAtenswj[ijkl]*Ucoefj[k,irep]
#           }
#         }
#       }
#     }
#   }
#   
#   return(DASarray) 
#   
# }
# 
# #  ------------------------------------------------------------------------
# 
# inprodDwj <- function(nXbasisw, nWbasisw, nUbasisj, nAbasisj, 
#                       nrep, Avecj, Ucoefj, BAtenswj) {
#   ncum <- cumprod(c(nAbasisj, nUbasisj, nWbasisw, nXbasisw))
#   DBSarray <- array(0,c(nXbasisw,nrep,nWbasisw))
#   for (irep in 1:nrep) {
#     for (i in 1:nXbasisw) {
#       for (j in 1:nWbasisw) {
#         for (k in 1:nUbasisj) {
#           for (l in 1:nAbasisj) {
#             ijkl <- (i-1)*ncum[3] + (j-1)*ncum[2] + (k-1)*ncum[1] + l
#             DBSarray[i,irep,j] <- DBSarray[i,irep,j] + 
#               Avecj[l]*BAtenswj[ijkl]*Ucoefj[k,irep]
#           }
#         }
#       }
#     }
#   }
#   
#   return(DBSarray)
#   
# }
# 
# #  ------------------------------------------------------------------------
#   
# inprodijw <- function(nXbasisw, nWbasisw, nUbasisj, nAbasisj, 
#                      nrep, Bvecw, Avecj, Ucoefj, BAtenswj) {
#   ncum <- cumprod(c(nAbasisj, nUbasisj, nWbasisw, nXbasisw))
#   Smatjw <- matrix(0,nXbasisw,nrep)
#   for (irep in 1:nrep) {
#     for (i in 1:nXbasisw) {
#       for (j in 1:nWbasisw) {
#         for (k in 1:nUbasisj) {
#           for (l in 1:nAbasisj) {
#             ijkl <- (i-1)*ncum[3] + (j-1)*ncum[2] + (k-1)*ncum[1] + l
#             Smatjw[i,irep] <- Smatjw[i,irep] + 
#                    Bvecw[j]*Avecj[l]*BAtenswj[ijkl]*Ucoefj[k,irep]
#           }
#         }
#       }
#     }
#   }
#   return(Smatjw)
# }
  