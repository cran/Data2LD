Data2LD.ISE <- function(XbasisList, modelList, coefList, coef, 
                        Rmat, Smat,  nrep, nforce, rhoVec=rep(0.5,nvar), AtensorList) {
  #  D2LD ... stands for "Data to Linear Dynamics"
  #  D2LD_ISE computes the value of the penalty term, the integrated
  #  squared difference between the right and left sides of a linear
  #  differential equation of the form
  #      D^m x_i <- sum_k^d sum_j^{m_k} beta_{kj} D^{j-1} x_k + 
  #                sum_{\ell}^{M_i} \alpha_{\ell,i} u_{\ell,i}, 
  #      i=1,...,d
  #  where
  #  and where each coefficient is expanded in terms of its own number of
  #  B-spline basis functions:
  #      \beta_{ij}(t)      <- \bbold_{ij}'     \phibold_{ij}(t),
  #      \alpha_{\ell,i}(t) <- \abold_{\ell,i}' \psibold_{\ell,i}(t)
  #  This version inputs AtensorList as an argument.
  #
  #  This version approximates the integrals in the penalty terms by using 
  #  inprod_basis to compute the cross-product matrices for the  
  #  \beta-coefficient basis functions and the corresponding derivative of 
  #  the x-basis functions,and the cross-product matrices for the 
  #  \alpha-coefficients and the corresponding U functions.  
  #  These are computed upon the first call to D2LD, and then retained 
  #  for subsequent calls by using the R-Cache command.  
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
  #
  #  MODELLIST ...  A list of length NVAR. Each list contains a 
  #                struct object with members:              
  #                XList ... list of length number of homogeneous terms
  #                          Each list contains a struct object with members:
  #                          WfdPar ... a fdPar object for the coefficient
  #                          variable   ... the index of the variable
  #                          derivative ... the order of its derivative
  #                          npar ... if coefficient estimated, its location
  #                                   in the composite vector 
  #                FList ... list of length number of forcing terms
  #                          Each list contains a struct object with members:
  #                          AfdPar ... an fdPar object for the coefficient
  #                          Ufd    ... an fd object for the forcing function
  #                          npar ... if coefficient estimated, its location
  #                                   in the composite vector 
  #                order     ... the highest order of derivative
  #                name      ... a  tag for the variable
  #                nallXterm ... the number of homogeneous terms
  #                nallFterm ... the number of forcing functions
  #
  #  COEFLIST  ... A list of length NCOEF.  Each list contains a
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
  #  RHOVEC    ... A vector of length NVAR containing values in [0,1].  
  
  #                The data sums of squares are weighted by P and 
  #                the roughness penalty by 1-rho.
  #
  #  COEF      ... The coefficient matrix
  #
  #  RMAT      ... The penalty matrix for the homogeneous part of L
  #
  #  SMAT      ... The penalty matrix for the forcing part of L
  #
  #  NREP      ... The number of replications
  #
  #  NFORCE    ... The vector containing the number of forcing functions
  #                per variable.
  #
  #  ATENSORLIST ...
  
  #  Last modified 2 January 2018
  
  #  Set up a vector NCOEFVEC containing number of coefficients used
  #  for the expansion of each variable
  
  nvar <- length(modelList)
  
  ncoefvec <- rep(0,nvar)
  for (ivar in 1:nvar) {
    ncoefvec[ivar] <- XbasisList[[ivar]]$nbasis
  }
  
  #  get the width of the time domain
  
  Xrange <- XbasisList[[1]]$rangeval
  T      <- Xrange[2] - Xrange[1]
  
  ISE1 <- rep(0,nvar)
  ISE2 <- rep(0,nvar)
  ISE3 <- rep(0,nvar)
  for (irep in 1:nrep) {
    ISE1i <- 0
    ISE2i <- 0
    ISE3i <- 0
    coefi <- coef[,irep]
    for (ivar in 1:nvar) {
      ISE1i <- ISE1i + t(coefi) %*% Rmat %*% coefi
      modelListi <- modelList[[ivar]]
      weighti <- modelListi$weight
      nforcei  <- nforce[ivar]
      if (nforcei > 0) {
        ISE2i <- ISE2i + 2 %*% t(coefi) %*% Smat[,irep]
        ISE3i <- 0
        for (jforce in 1:nforcei) {
          modelListij <- modelListi$FList[[jforce]]
          ncoefj     <- modelListij$ncoef
          coefListj  <- coefList[[ncoefj]]
          Avecj      <- coefListj$parvec
          nAbasisj   <- length(Avecj)
          factorj    <- modelListij$factor
          Ufdj       <- modelListij$Ufd
          Ubasisj    <- Ufdj$basis
          Ucoefj     <- Ufdj$coef
          nUbasisj   <- Ubasisj$nbasis
          funj       <- coefListj$fun
          funtypej   <- !(is.basis(funj) || is.fd(funj) ||is.fdPar(funj))
          for (kforce in 1:nforcei) {
            modelListik <- modelListi$FList[[kforce]]
            ncoefk    <- modelListik$ncoef
            coefListk <- coefList[[ncoefk]]
            Aveck     <- coefListk$parvec
            factork   <- modelListik$factor
            nAbasisk  <- length(Aveck)
            Ufdk      <- modelListik$Ufd
            Ubasisk   <- Ufdk$basis
            Ucoefk    <- Ufdk$coef
            nUbasisk  <- Ubasisk$nbasis
            # Atensijk  <- AtensorList[[ivar]]
            Atensijk  <- AtensorList[[ivar]][[jforce]][[kforce]]
            funk      <- coefListk$fun
            funtypek <- !(is.basis(funk) || is.fd(funk) ||is.fdPar(funk))
            if (funtypej || funtypek) {
              ISE3i <- inprod.basis.Data2LD(Ufdk,  Ufdj, 
                                           coefListk, coefListj,  0,   0)
              
            } else {
              ISE3i <- 0
              Atensijk <- AtensorList[[ivar]][[jforce]][[kforce]]
              ncum <- cumprod(c(nAbasisk, nUbasisk, nAbasisj, nUbasisj))
              for (i in 1:nUbasisj) {
                for (j in 1:nAbasisj) {
                  for (k in 1:nUbasisk) {
                    for (l in 1:nAbasisk) {
                      ijkl <- (i-1)*ncum[3] + (j-1)*ncum[2] + (k-1)*ncum[1] + l
                      ISE3i <- ISE3i +  
                        Ucoefj[i,irep]*Avecj[j]*Ucoefk[k,irep]*Aveck[l]*Atensijk[ijkl]
                    }
                  }
                }
              }
            }
            ISE3i <- factorj*factork*rhoVec[ivar]*ISE3i/T
          }
        }
      } else {
        ISE2i <- 0
        ISE3i <- 0
      }
      ISE1[ivar] <- ISE1[ivar] + weighti*ISE1i
      ISE2[ivar] <- ISE2[ivar] + weighti*ISE2i
      ISE3[ivar] <- ISE3[ivar] + weighti*ISE3i
    }
  }
  ISE <- (ISE1 + ISE2 + ISE3)/nrep
  
  return(ISE)
  
}
