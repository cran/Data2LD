Data2LD <- function(yList, XbasisList, modelList, coefList,
                    rhoVec=0.5*rep(1,nvar), summary=TRUE) {
#  Data2LD ... stands for "Data to Linear Dynamics"
#  It approximates the data in argument YLIST by one or smooth 
#  functions x_i, i=1,...,d.  This approximation is defined by a set
#  of linear differential or algebraic equations defined by a set of
#  parameters some of which require estimation from the data. 
#
#  The approximation minimizes the sum of squared residuals, expressed as 
#  follows in Latex notation:
#
#    H(\theta) <- \sum_i^d \sum_j^n \sum_\ell^N [y_{ij \ell} - x_i(t_j)]^2
#
#  where:
#  i    <- 1,...,d indexes equations in a system differential and algebraic
#                 equations.
#  j    <- 1,...,n indexes times of observation of a variable
#  \ell <- 1,...,N indexes replications of observations
#
#  But there is additional flexibility not captured in this expression:
#  1.  Only a subset of the variables may be observed, so that not all 
#      values of index i are actually used.
#  2.  The number and location of times of observation t_j  can vary
#      from one observed variable to another.
#  using a roughness penaltylinear differential operator that depends
#  on unknown parameters in list COEFLIST, which is described below.
#
#  The fitting functions x_i(t) in turn, here assumed to be defined over
#  the interval [0,T], are defined by basis function expansions:
#
#         x_i(t_j) <- \sum_k^K_i c_{ik}(\theta|\rho) \phi_{ik}(t_j)
#
#  where \phi_{ik}(t) is the kth function in a system of K_i basis 
#  functions used to approximate the ith variable.
#  The number of K_i basis functions and the type of basis function system
#  can vary from one variable to another.  This information is contained
#  in argument BASISLIST described below.
#
#  The coefficients c_{ik}(\theta|\rho) defining the smoothing functions  
#  are functions of the unknown parameters in vector \theta that define
#  the differential and algebraic equations and that require estimation
#  from the data.  The smoothing parameter $\rho is a value in the interval
#  [0,1).
#  The coefficient functions c_{ik}(\theta|\rho) minimize the inner 
#  least squares criterion, expressed here for simplicity for a single  
#  variable or equation:
#
#    J(c|\theta) <- (1-\rho) \sum_j^n [y_j - c_k^K \phi_k(t_j)]^2/n +
#                  \rho \int_0^T {L(\theta)x(t)]^2 dt/T
#
#  Linear differential operator L(\theta) is a linear differential or
#  algebraic equation rewritten as an operator by subtracting the right
#  side of the equation from the left, so that a solution x of the 
#  equation satisfies Lx <- 0.
#  Each L in a system of equation depends on one or more parameters that 
#  need to be estimated and are contained in parameter vector \theta.
#
#  The linear differential equation is of the form
#      D^m x_i <- sum_k^d sum_j^{m_k} beta_{kj}(t) D^{j-1} x_k + 
#                sum_f^{F_i} \alpha_{fi}(t) u_{f,i}, 
#      i=1,...,d,  f=1,...,F_i
#  where
#  and where each coefficient is expanded in terms of its own number of
#  B-spline basis functions:
#      \beta_{ij}(t)  <- \bbold_{ij}" \phibold_{ij}(t),
#      \alpha_{fi}(t) <- \abold_{fi}" \psibold_{fi}(t)
#
#  As smoothing parameter \rho increases toward its upper limit of 1,
#  the second roughness penalty term in J is more and more emphasized
#  and the fit to the data less and less emphasized, so that x_i is
#  required to approach a solution to its respective equation.
#
#  The highest order of derivative can vary from one equation to another,
#  and in particular can be larger than 1.
#  Each right side may contain contributions from all variables in the
#  equation.  Morever, these contributions may be from all permissible
#  derivatives of each equation, wher "permissible" means up to one less 
#  than the highest order in the variable.
#  In addition, each equation may be forced by a number of forcing terms
#  that can vary from equation to equation.
#
#  This version approximates the integrals in the penalty terms by using 
#  inprod_basis to compute the cross-product matrices for the  
#  \beta-coefficient basis functions and the corresponding derivative of 
#  the x-basis functions,and the cross-product matrices for the 
#  \alpha-coefficients and the corresponding U functions.  
#  These are computed upon the first call to D2LD, and then retained 
#  for subsequent calls by using the persistent command.  This technique
#  for speeding up computaton is called memoization.  
#  See lines 224 and 352 to 388 for this code.
#
#  The structure of the model is defined in list MODELLIST, which is
#  described below.
#
#  This version disassociates coefficient functions from equation 
#  definitions to allow some coefficients to be used repeatedly and for
#  both homogeneous and forcing terms.  It requires an extra argument
#  COEFLIST that contains the coefficients and the position of their
#  coefficient vectors in vector THETA.
#
#  ------------------------------------------------------------------------
#
#  Arguments:
#
#  YLIST     ... A list of length NVAR.  Each list contains in turn
#                a struct object with fields:
#                  "argvals" is a vector of length n_i of observation times
#                  "y" contains a matrix with n_i rows and NREP columns.
#                The number of columns must be the same for all variables,
#                except that, if a list is empty, that variable is taken to 
#                be not observed.
#
#  BASISLIST ... A list array of length NVAR.  Each member contains in turn
#                a functional data object or a BASIS object.
#
#  MODELLIST ... A list of length NVAR. Each list contains a
#                list object with members:              
#                XList ... list of length number of homogeneous terms
#                          Each list contains a struct object with members:
#                          fun        ... a fdPar object for the coefficient
#                          variable   ... the index of the variable
#                          derivative ... the order of its derivative
#                          ncoef      ... if coefficient estimated, its location
#                                         in the composite vector 
#                          factor     ... a scalar multiplier (def. 1)
#                          estimate   ... 0, held fixed, otherwise, estimated
#                FList ... list of length number of forcing terms
#                          Each list contains a struct object with members:
#                          AfdPar ... an fdPar object for the coefficient
#                          Ufd    ... an fd object for the forcing function
#                          ncoef  ... if coefficient estimated, its location
#                                     in the composite vector 
#                          factor ... a scalar multiplier (def. 1)
#                          estimate   ... 0, held fixed, otherwise, estimated
#                order ... the highest order of derivative
#                name  ... a  tag for the variable
#                nallXterm ... the number of homogeneous terms
#                nallFterm ... the number of forcing functions
#
#  COEFLIST  ... A list of length NCOEF.  Each member contains a
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
#
#  RHOVEC    ... A vector of length NVAR containing values in [0,1].
#                The data sums of squares are weighted by P and 
#                the roughness penalty by 1-P.
#
#  LOAD_TENSOR ... If nonzero, attempt to load the lists
#                BtensorList, BAtensorList and AtensorList.  These must
#                have set up before any call to D2LD and saves as 
#                .mat files with names BtensorList.mat, BAtensorList.mat
#                and AtensorList.mat.  
#                For information on how these are set up, see the functions 
#                Btensorfn, BAtensorfn and Atensorfn.
#
#  ------------------------------------------------------------------------
#
#  Output objects (d <- number of equations, 
#                  NTHETA is the total number of estimated parameters):
#
#  MSE     ... The weighted mean squared errors computed over the variables
#              with data.    
#  DpMSE   ... The gradient of the objective function MSE with respect to 
  #            the estimated parameters.
#  D2ppMSE ... A square symmetric matrx of order NTHETA that contains
#              the second partial derivatives of the objective function.
#  XFDPARLIST ... A list of length d containing functional parameter
#              objects of class fdPar for the estimated functions x_i(t).
#  DF      ... An equivalent degrees of freedom value 
#                   df <- trace(2*YM - YM*YM) where YM is the matrix
#              fitMap described below.
#  GCV     ... The generalized cross-validation measure.  The value of
#              \rho corresponding to the minimim of GCV across values of
#              smoothing parameter \rho is often chose for an automatic 
#              data-driven level of smoothing.
#  ISE     ... The sum across variables of the integrated squared value of
#              the differential operator.  This value multiplied by \rho
#              and divided by T, the width of the domain, is the second
#              term the objective function.
#  RMAT    ... A square symmetric matrix of order equal to the sum of the
#              coefficients in the basis function expansions of the 
#              variables.  This matrix defines the size of the second term
#              in the objective function.
#  SMAT    ... Either a vector of length equal to the order of RMAT if
#              there is only one replication, or a matrix with number of
#              columns equal to the number of replications NREP.
#  fitMap  ... A matrix with number of rows equal to the total number of
#              coefficients in the basis expansions of variables and 
#              number of columns equal the total number of observations.
#              This matrix is the linear map from the data to the 
#              combined coefficients.

#  Last modified 2 January 2018

#  ------------------------------------------------------------------------
#                    Check modelList 
#  ------------------------------------------------------------------------

modelList <- modelCheck(modelList, coefList)

nvar <- length(modelList)

#  ------------------------------------------------------------------------
#                    Check coefList 
#  ------------------------------------------------------------------------
  
coefCheckList <- coefCheck(coefList)
coefList      <- coefCheckList$coefList 
ntheta        <- coefCheckList$ntheta

#  ------------------------------------------------------------------------
#  Store the number of forcing functions for each variable and load the
#  starting position in the composite vector of estimated coefficients.
#  ------------------------------------------------------------------------

nhomog   <- rep(0,nvar)
nforce   <- rep(0,nvar)
nthetaH  <- 0
nthetaF  <- 0
for (ivar in 1:nvar) {
  modelListi <- modelList[[ivar]]
  nXterm <- modelListi$nallXterm
  #  process homogeneous terms
  if (nXterm > 0) {
    nhomog[ivar] <- nXterm
    for (iterm in 1:nXterm) {
      XListi    <- modelListi$XList[[iterm]]
      ncoefi    <- XListi$ncoef
      coefListi <- coefList[[ncoefi]]
      if (coefListi$estimate == TRUE) {
        nthetaH <- nthetaH + length(coefListi$parvec)
      }
    }
  }
  nFterm <- modelListi$nallFterm
  #  process forcing terms
  if (nFterm > 0) {
    nforce[ivar] <- nFterm
    for (iterm in 1:nFterm) {
      FListi <- modelListi$FList[[iterm]]
      ncoefi <- FListi$ncoef
      coefListi <- coefList[[ncoefi]]
      if (coefListi$estimate) {
        nthetaF <- nthetaF + length(coefListi$parvec)
      }
    }
  }
}

#  ------------------------------------------------------------------------
#                        Check YLIST
#  ------------------------------------------------------------------------

ycheckList <- yListCheck(yList, nvar)
nrep    <- ycheckList$nrep
nvec    <- ycheckList$nvec 
dataWrd <- ycheckList$dataWrd

#  ------------------------------------------------------------------------
#                        Check structure of BASISLIST
#  ------------------------------------------------------------------------

#  Retrieve basis object for each variable from BASISLIST and install it
#  in basis List.
#  And set up a vector NCOEFVEC containing number of coefficients used
#  for the expansion of each variable

if (!is.list(XbasisList)) {
  stop("BASISLIST is not a list.")
}

if (length(XbasisList) != nvar) {
  stop("BASISLIST is not of length NVAR.")
}

errwrd <- FALSE
ncoefvec <- matrix(0,nvar,1)
for (ivar in 1:nvar) {
    basisi <- XbasisList[[ivar]]
    if (!is.basis(basisi)) {
      warning(paste("BASIS is not a BASIS object for variable ",ivar,"."))
      errwrd <- TRUE
    } else {
        ncoefvec[ivar] <- basisi$nbasis
        XbasisList[[ivar]] <- basisi
    }
}

if (errwrd) {
    stop("One or more terminal err encountered in BASISLIST.")
}

#  get the width of the time domain

Xrange <- XbasisList[[1]]$rangeval
T      <- Xrange[2] - Xrange[1]

#  check rhoVec

if (length(rhoVec) != nvar) {
    stop("RHOVEC not of length NVAR.")
}

for (ivar in 1:nvar) {
    if (rhoVec[ivar] > 1 || rhoVec[ivar] < 0) {
        stop(paste("P is not in [0,1] for variable ",ivar,"."))
    }
}

#  ------------------------------------------------------------------------
#  Compute crossproduct matrices for tensor products of 
#  D^j X-basis functions and W-basis functions, j=0,,nderivvec 
#  if not already set up
#  ------------------------------------------------------------------------

Btensorfn  <- addMemoization( Btensorfn)
BAtensorfn <- addMemoization(BAtensorfn)
Atensorfn  <- addMemoization( Atensorfn)

BtensorList  <-   Btensorfn(XbasisList, modelList, coefList)
BAtensorList <-  BAtensorfn(XbasisList, modelList, coefList)
AtensorList  <-   Atensorfn(            modelList, coefList)

#  ------------------------------------------------------------------------
#  set up list for matrices of basis function values
#  ------------------------------------------------------------------------

basismatList <- list(nvar,1)
for (ivar in 1:nvar) {
    if (!is.null(yList[[ivar]])) {
        yListi  <- yList[[ivar]]
        basisi  <- XbasisList[[ivar]]
        argvals <- as.vector(yListi$argvals)
        basismati <- eval.basis(argvals, basisi)
        basismatList[[ivar]] <- basismati
    }
}

#  ------------------------------------------------------------------------
#                  Compute coefficient matrix Bmat
#  ------------------------------------------------------------------------

nsum     <- sum(nvec)
ncoefsum <- sum(ncoefvec)
Bmat     <- matrix(0,ncoefsum,ncoefsum)
basismat <- Matrix(0,nsum,ncoefsum)
ymat     <- matrix(0,nsum,nrep)
m2 <- 0
n2 <- 0
for (ivar in 1:nvar) {
    m1  <- m2 + 1
    m2  <- m2 + ncoefvec[ivar]
    ind <- m1:m2
    if (dataWrd[ivar]) {
        modelListi <- modelList[[ivar]]
        weighti    <- modelListi$weight
        n1   <- n2 + 1
        n2   <- n2 + nvec[ivar]
        indn <- n1:n2
        yListi   <- yList[[ivar]]
        ymat[indn,] <- as.matrix(yListi$y)
        basismati   <- basismatList[[ivar]]
        Bmat[ind,ind] <- weighti*(1-rhoVec[ivar])*crossprod(basismati)/nvec[ivar]
        basismat[indn,ind] <- basismati
    }
}

#  ------------------------------------------------------------------------
#      Compute roughness penalty matrices Rmat(theta) and Smat(theta)
#      and their derivatives with respect to estimated parameters
#  ------------------------------------------------------------------------

#  Matrices R and DR

#  There are parameters to estimate defining the homogeneous terms

Data2LDRList <- Data2LD.Rmat(XbasisList, modelList, coefList, rhoVec, ntheta,
                            BtensorList)
Rmat <- Data2LDRList$Rmat 

if (nthetaH > 0) {
  DRarray <- Data2LDRList$DRarray
} else {
  #  No estimated parameters are involved in homogeneous terms
  DRarray <- NULL
}

#  Matrices S and DS for variables having forcing functions

Data2LDSList <- Data2LD.Smat(XbasisList, modelList, coefList, rhoVec, 
                            ntheta, BAtensorList, nrep, nforce)
Smat <- Data2LDSList$Smat

if (nthetaF > 0) {
  DSarray <- Data2LDSList$DSarray
} else {
  #  No estimated parameters are involved in forcing terms
  DSarray <- NULL
}

Cmat <- Bmat + Rmat

#  ------------------------------------------------------------------------
#                     Set up right side of equation
#  ------------------------------------------------------------------------

Dmat <- matrix(0,ncoefsum,nrep)
m2 <- 0
for (ivar in 1:nvar) {
    m1  <- m2 + 1
    m2  <- m2 + ncoefvec[ivar]
    ind <- m1:m2
    if (dataWrd[ivar]) {
        modelListi <- modelList[[ivar]]
        weighti    <- modelListi$weight
        yListi     <- yList[[ivar]]
        basismati  <- basismatList[[ivar]]
        yi         <- yListi$y
        ni         <- nvec[ivar]
        Dmat[ind,] <- weighti*(1-rhoVec[ivar])*t(basismati) %*% yi/ni
    }
}

#  --------------------------------------------------------------------
#               Compute coefficient matrix
#  --------------------------------------------------------------------

if (is.null(Smat)) {
  coef <- solve(Cmat,Dmat)
} else {
  coef <- solve(Cmat,(Dmat-Smat))
}

#  ------------------------------------------------------------------------
#  Compute the vector of unpenalized error sum of squares, MSE, 
#  the sum of which is the outer objective function H(\theta|\rho).
#  Each MSE_i is normalized by dividing by NREP and by the n_i"s.
#  ------------------------------------------------------------------------

xmat <- basismat %*% coef

MSE    <- 0
SSEtot <- 0
m2     <- 0
for (ivar in 1:nvar) {
  if (dataWrd[ivar]) {
    modelListi <- modelList[[ivar]]
    weighti    <- modelListi$weight
    m1 <- m2 + 1
    m2 <- m2 + nvec[ivar]
    xmati  <- xmat[m1:m2,]
    ymati  <- yList[[ivar]]$y
    rmati  <- ymati - xmati
    SSEi   <- sum(rmati^2)
    SSEtot <- SSEtot + weighti*SSEi
    MSE    <- MSE + weighti*SSEi/nrep/nvec[ivar]
  }
}

#  compute residual variance

Rvar <- SSEtot/nsum

#  ------------------------------------------------------------------------
#                     Compute summary values
#  ------------------------------------------------------------------------

if (summary) {
  y2cFac <- t(basismat)
  m2 <- 0
  n2 <- 0
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    weighti    <- modelListi$weight
    m1  <- m2 + 1
    m2  <- m2 + ncoefvec[ivar]
    ind <- m1:m2
    if (dataWrd[ivar]) {
      n1   <- n2 + 1
      n2   <- n2 + nvec[ivar]
      indn <- n1:n2
      y2cFac[ind,indn] <- weighti*(1-rhoVec[ivar])*y2cFac[ind,indn]/nvec[ivar]
    }
  }
  RgtFactor <- solve(Cmat,y2cFac)
  fitMap <- Matrix(basismat %*% RgtFactor)
  
  #  Use fitMap to compute a equivalent degrees of freedom measure
  
  df <- sum(diag(fitMap))
  
  #  set up the functional data object for the output smooth of the data 
  
  XfdParList <- vector("list",nvar)
  m2 <- 0
  for (ivar in 1:nvar) {
    m1  <- m2 + 1
    m2  <- m2 + ncoefvec[ivar]
    ind <- m1:m2
    Xbasisi <- XbasisList[[ivar]]
    Xfdobji <- fd(coef[ind,],Xbasisi)
    XfdParList[[ivar]] <- fdPar(Xfdobji)
  }
  
  #  compute GCV
  if (df < nsum) {
    gcv <- Rvar/((nsum - df)/nsum)^2
  } else {
    gcv <- NA
  }
  
  #  ------------------------------------------------------------------------
  #  Compute unpenalized error integrated squares, ISE, the penalty term
  #  loop through number of replications NREP
  #  ------------------------------------------------------------------------
  
  ISE <- Data2LD.ISE(XbasisList, modelList, coefList, coef, Rmat, Smat, 
                    nrep, nforce, rhoVec, AtensorList)
} else {
  df  <- NULL
  gcv <- NULL
  ISE <- NULL
  XfdParList <- NULL
  fitMap <- NULL
}

#  ------------------------------------------------------------------------
#       Compute total derivative of MSE wrt theta if required
#  ------------------------------------------------------------------------

#  ---------  Compute the total derivative  ----------------

#  Compute the partial derivatives of the coefficients with respect to the
#  estimated parameters,  dc/dtheta

Dcoef <- array(0, c(ncoefsum,ntheta,nrep))
for (itheta in 1:ntheta) {
  if (nthetaH > 0) {
    DRmati <- DRarray[,,itheta]
    if (nthetaF > 0) {
      for (irep in 1:nrep) {
        DRi <- -DRmati %*% as.matrix(coef[,irep])
        DSi <- matrix(0,ncoefsum,1)
        m2 <- 0
        for (ivar in 1:nvar) {
          m1 <- m2 + 1
          m2 <- m2 + ncoefvec[ivar]
          indi <- m1:m2
          offset <- ncoefvec[ivar]*(irep - 1 + (itheta-1)*nrep)
          DSi[indi,1] <- -DSarray[indi + offset]
        }
        Dcoef[,itheta,irep] <- solve(Cmat,(DRi + DSi))
      }
    } else {
      for (irep in 1:nrep) {
        DRi <- -DRmati %*% coef[,irep]
        Dcoef[,itheta,irep] <- solve(Cmat,DRi)
      }
    }
  } else {
    if (nthetaF > 0) {
      for (irep in 1:nrep) {
        DSi <- matrix(0,ncoefsum,1)
        m2 <- 0
        for (ivar in 1:nvar) {
          m1 <- m2 + 1
          m2 <- m2 + ncoefvec[ivar]
          indi <- m1:m2
          offset <- ncoefvec[ivar]*(irep - 1 + (itheta-1)*nrep)
          DSi[indi,1] <- -DSarray[indi + offset]
        }
        Dcoef[,itheta,irep] <- solve(Cmat,DSi)
      }
    }
  }
}

#  ------------------------------------------------------------------------
#            Compute the total theta-gradient of H
#  ------------------------------------------------------------------------

xmat    <- basismat %*% coef
DpMSE   <- matrix(0,ntheta, 1)
D2ppMSE <- matrix(0,ntheta, ntheta)
D2pyMSE <- matrix(0,ntheta, nsum)

m2 <- 0
for (ivar in 1:nvar) {
    if (dataWrd[ivar]) {
        m1 <- m2 + 1
        m2 <- m2 + nvec[ivar]
        modelListi <- modelList[[ivar]]
        weighti    <- modelListi$weight
        basismati  <- basismat[m1:m2,]
        xmati      <- xmat[m1:m2,]
        yVeci      <- yList[[ivar]]$y
        rmati      <- as.matrix(yVeci - xmati)
        for (irep in 1:nrep) {
          Dcoefi <- Dcoef[,,irep]
          BasDcoefi <- basismati %*% Dcoefi
          DpMSE   <- DpMSE   - 2*weighti*
                 (t(BasDcoefi) %*% as.matrix(rmati[,irep]))/nrep/nvec[ivar]
          D2ppMSE <- D2ppMSE + 2*weighti*
                 (t(BasDcoefi) %*% BasDcoefi)/nrep/nvec[ivar]
          temp <- 2*weighti*t(BasDcoefi)/nrep/nvec[ivar]
          D2pyMSE[,m1:m2] <- D2pyMSE[,m1:m2] - as.matrix(temp)
        }
    }
}

DpMSE   <- as.matrix(DpMSE)
D2ppMSE <- as.matrix(D2ppMSE)

if (summary) {
  DpDy      <- -solve(D2ppMSE,D2pyMSE)
  sigmasq   <- SSEtot/(nsum - df)
  Var.theta <- as.matrix(sigmasq*(DpDy %*% t(DpDy)))
}

if (summary) {
  return(list(MSE=MSE, DpMSE=DpMSE, D2ppMSE=D2ppMSE, XfdParList=XfdParList, 
              df=df, gcv=gcv, ISE=ISE, Var.theta=Var.theta, 
              Rmat=Rmat, Smat=Smat))
  
} else {
  return(list(MSE=MSE, DpMSE=DpMSE, D2ppMSE=D2ppMSE))
}

}

