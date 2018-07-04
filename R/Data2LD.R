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
#  These are computed upon the first call to Data2LD, and then retained
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
#                have set up before any call to Data2LD and saves as
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

#  Last modified 27 June 2018

#  ------------------------------------------------------------------------
#                    Check modelList
#  ------------------------------------------------------------------------

modelList <- modelCheck(modelList, coefList, report=FALSE)

nvar <- length(modelList)

#  ------------------------------------------------------------------------
#                    Check coefList
#  ------------------------------------------------------------------------

coefCheckList <- coefCheck(coefList, report=FALSE)
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
nrep       <- ycheckList$nrep
nvec       <- ycheckList$nvec
dataWrd    <- ycheckList$dataWrd

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

Ceigvals <- eigen(Cmat)$values

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
if (summary) {
  D2pyMSE <- matrix(0,ntheta, nsum)
}

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
          D2ppMSE <- D2ppMSE + 
            2*weighti*crossprod(BasDcoefi)/nrep/nvec[ivar]
          if (summary) {
            temp <- 2*weighti*t(BasDcoefi)/nrep/nvec[ivar]
            D2pyMSE[,m1:m2] <- D2pyMSE[,m1:m2] - as.matrix(temp)
          }
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

#   ---------------------------------------------------------------------------

Data2LD.Rmat <- function(XbasisList, modelList, coefList, rhoVec=rep(0.5,nvar),
                         ntheta,BtensorList) {
  #  Data2LD_Rmat computes the penalty matrix R associated with the homogeneous
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
  #  These are computed upon the first call to Data2LD4, and then retained
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

  #  Last modified 1 May 2018

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
      order  <- modelListi$order
      Rmatii <- matrix(inprod(Xbasisi, Xbasisi, order, order),nXbasisi,nXbasisi)
    } else {
      Btensii <- BtensorList[[ivar]][[nXtermi+1]][[nXtermi+1]]
      if (any(is.na(Btensii))) stop("NAs in Btensii")
      Rmatii  <- RmatFn(nXbasisi, 1, nXbasisi, 1, 1, 1, Btensii)
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
          Rmatwx  <- RmatFn(nXbasisw, nWbasisw, nXbasisx, nWbasisx,
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
        Rmatiw  <- RmatFn(nXbasisi, 1, nXbasisw, nWbasisw, 1, Bvecw, Btensiw)
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
        Bvecw    <- coefListw$parvec
        factorw  <- modelListiw$factor
        nWbasisw <- length(Bvecw)
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
            DRarraywx <- DRarrayFn(nXbasisw, nWbasisw, nXbasisx, nWbasisx,
                                   Bvecx, Btenswx)
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
          DRarraywi <- DRarrayFn(nXbasisw, nWbasisw, nXbasisi, 1, 1.0, Btenswi)
        }
        #  rescale the inner product
        DRarraywi <- factorw*rhoVec[ivar]*DRarraywi/T
        temp <- DRarray[indw,indi,indthw,drop=FALSE]
        temp <- temp - weighti*DRarraywi
        DRarray[indw,indi,indthw] <- temp
        DRarray[indi,indw,indthw] <- aperm(temp,c(2,1,3))
      }
    }
  }
  return(list(Rmat=Rmat, DRarray=DRarray))
}

#  ----------------------------------------------------------------------------

Data2LD.Smat <- function(XbasisList, modelList, coefList, rhoVec=rep(0.5,nvar),
                         ntheta, BAtensorList, nrep, nforce) {
  #  Data2LD ... stands for "Data to Linear Dynamics"
  #  Data2LD_S computes the penalty matrix S associated with the forcing
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
  #  These are computed upon the first call to Data2LD, and then retained
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
  #  Last modified 23 May 2018

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
            Smatjw   <- SmatFn(nXbasisw, nWbasisw, nUbasisj, nAbasisj, nrep,
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
  #  in parvec(1:nthetaL)
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
              BAtenswj   <- BAtensorList[[ivar]][[iw]][[jforce]]
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

#  ----------------------------------------------------------------------------

Data2LD.ISE <- function(XbasisList, modelList, coefList, coef,
                        Rmat, Smat,  nrep, nforce, rhoVec=rep(0.5,nvar),
                        AtensorList) {
  #  Data2LD ... stands for "Data to Linear Dynamics"
  #  Data2LD_ISE computes the value of the penalty term, the integrated
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
  #  These are computed upon the first call to Data2LD, and then retained
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

  #  Last modified 2 May 2018

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

#  ----------------------------------------------------------------------------

inprod.basis.Data2LD <- function(fdobj1, fdobj2, coefList1, coefList2,
                                 Lfdobj1=int2Lfd(0), Lfdobj2=int2Lfd(0),
                                 EPS=1e-6, JMAX=15, JMIN=5)
{
  #  INPROD.BASIS.DATA2LD  Computes matrix of inner products of bases by numerical
  #    integration using Romberg integration with the trapezoidal rule.
  #  This version multiplies bases by respective coefficient functions
  #    betafn1 and betafn2, and is intended for use within function
  #    Data2LD.  It produces a matrix.

  #  Arguments:
  #  FDOBJ1 and FDOBJ2:    These are functional data objects.
  #  COEFLIST1 and COEFLIST2 ...  List objects for BASIS1 and BASIS2
  #  containing members:
  #               parvec   ... a vector of parameters
  #               estimate ... 0, held fixed, otherwise, estimated
  #               coeftype ... homogeneous or forcing
  #               fun      ... functional basis, fd, or fdPar object,
  #                            or a struct object for a general function
  #                            with fields:
  #                 fd      ... function handle for evaluating function
  #                 Dfd     ... function handle for evaluating
  #                              derivative with respect to parameter
  #                 more    ... object providing additional information for
  #                             evaluating coefficient function
  #  However, these may also be real constants that will be used as fixed coefficients.  An
  #  example is the use of 1 as the coefficient for the left side of the equation.
  #  Lfdobj1 and Lfdobj2:  order of derivatives for inner product of
  #               fdobj1 and fdobj2, respectively, or functional data
  #               objects defining linear differential operators
  #  EPS    convergence criterion for relative stop
  #  JMAX   maximum number of allowable iterations
  #  JMIN   minimum number of allowable iterations

  #  Return:
  #  A matrix of NREP1 by NREP2 of inner products for each possible pair
  #  of functions.

  #  Last modified 12 January 2018

  #  Determine where fdobj1 and fdobj2 are basis or fd objects, define
  #  BASIS1 and BASIS2, and check for common range

  errwrd <- FALSE
  if (is.basis(fdobj1)) {
    basis1  <- fdobj1
    nbasis1 <- basis1$nbasis - length(basis1$dropind)
    ndim1   <- nbasis1
  } else {
    if (is.fd(fdobj1)) {
      basis1  <- fdobj1$basis
      nbasis1 <- basis1$nbasis - length(basis1$dropind)
      ndim1   <- 1
    } else {
      errwrd <- TRUE
      warning("First argument is not a basis or an fd object.")
    }
  }

  if (is.basis(fdobj2)) {
    basis2  <- fdobj2
    nbasis2 <- basis2$nbasis - length(basis2$dropind)
    ndim2   <- nbasis2
  } else {
    if (is.fd(fdobj2)) {
      basis2  <- fdobj2$basis
      nbasis2 <- basis2$nbasis - length(basis2$dropind)
      ndim2    <- 1
    } else {
      errwrd <- TRUE
      warning("Second argument is not a basis or an fd object.")
    }
  }

  if (errwrd) stop("Terminal error encountered.")

  #  get coefficient vectors

  if (is.numeric(coefList1)) {
    bvec1 <- coefList1
    conbasis <- create.constant.basis(basis1$rangeval)
    List1 <- list(fun=conbasis, parvec=bvec1, estimate=FALSE)
    coefList1 <- List1
  } else {
    bvec1 <- coefList1$parvec
  }

  if (is.numeric(coefList2)) {
    bvec2 <- coefList2
    conbasis <- create.constant.basis(basis2$rangeval)
    List2 <- list(fun=conbasis, parvec=bvec2, estimate=FALSE)
    coefList2 <- List2
  } else {
    bvec2 <- coefList2$parvec
  }

  #  Set up beta functions

  if (!is.basis(coefList1$fun) && !is.fd(coefList1$fun) && !is.fdPar(coefList1$fun)) {
    type1   <- TRUE
    betafd1 <- coefList1$fun$fd
    more1   <- coefList1$fun$more
  } else {
    type1   <- FALSE
    fdobj   <- coefList1$fun
    if (is.basis(fdobj)) {
      betafd1 <- fd(coefList1$parvec,fdobj)
    }
    if (is.fd(fdobj)) {
      betafd1 <- fd(coefList1$parvec, fdobj$basis)
    }
    if (is.fdPar(fdobj)) {
      betafd1 <- fd(coefList1$parvec, fdobj$fd$basis)
    }
  }

  if (!is.basis(coefList2$fun) && !is.fd(coefList2$fun) && !is.fdPar(coefList2$fun)) {
    type2   <- TRUE
    betafd2 <- coefList2$fun$fd
    more2   <- coefList2$fun$more
  } else {
    type2   <- FALSE
    fdobj   <- coefList2$fun
    if (is.basis(fdobj)) {
      betafd2 <- fd(coefList2$parvec,fdobj)
    }
    if (is.fd(fdobj)) {
      betafd2 <- fd(coefList2$parvec, fdobj$basis)
    }
    if (is.fdPar(fdobj)) {
      betafd2 <- fd(coefList2$parvec, fdobj$fd$basis)
    }
  }

  #  check for any knot multiplicities in either argument

  knotmult <- numeric(0)
  if (type1 == "bspline") knotmult <- knotmultchk(basis1, knotmult)
  if (type2 == "bspline") knotmult <- knotmultchk(basis2, knotmult)

  #  Modify RNGVEC defining subinvervals if there are any
  #  knot multiplicities.

  rng <- basis1$rangeval
  if (length(knotmult) > 0) {
    knotmult <- sort(unique(knotmult))
    knotmult <- knotmult[knotmult > rng[1] && knotmult < rng[2]]
    rngvec   <- c(rng[1], knotmult, rng[2])
  } else {
    rngvec <- rng
  }

  #  -----------------------------------------------------------------
  #                   loop through sub-intervals
  #  -----------------------------------------------------------------

  nrng <- length(rngvec)
  for (irng  in  2:nrng) {
    rngi <- c(rngvec[irng-1],rngvec[irng])
    #  change range so as to avoid being exactly on
    #  multiple knot values
    if (irng > 2   ) rngi[1] <- rngi[1] + 1e-10
    if (irng < nrng) rngi[2] <- rngi[2] - 1e-10

    #  set up first iteration

    iter  <- 1
    width <- rngi[2] - rngi[1]
    JMAXP <- JMAX + 1
    h <- rep(1,JMAXP)
    h[2] <- 0.25
    s <- array(0,c(JMAXP,ndim1,ndim2))
    #sdim <- length(dim(s))
    x <- rngi
    nx <- 2
    #  first argument
    if (is.basis(fdobj1)) {
      if (type1) betamat1  <- matrix(betafd1(x, bvec1, more1),nx,nbasis1)
      else       betamat1  <- matrix(eval.fd(x, betafd1),     nx,nbasis1)
      basismat1 <- eval.basis(x, basis1, Lfdobj1) * betamat1
    } else {
      if (type1) {
        betamat1  <- matrix(betafd1(x, bvec1, more1),nx,nbasis1)
      } else {
        betamat1  <- matrix(eval.fd(x, betafd1),nx,nbasis1)
      }
      basismat1 <- eval.fd(x, basis1, Lfdobj1) * betamat1
    }
    #  second argument
    if (is.basis(fdobj2)) {
      if (type2) {
        betamat2 <- matrix(betafd2(x, bvec2, more2),nx,nbasis2)
      } else {
        betamat2 <- matrix(eval.fd(x, betafd2), nx, nbasis2)
      }
      basismat2 <- eval.basis(x, basis2, Lfdobj2) * betamat2
    } else {
      if (type2) {
        betamat2 <- matrix(betafd2(x, bvec2, more2),nx,nbasis2)
      } else {
        betamat2 <- matrix(eval.fd(x, betafd2),nx,ndim2)
      }
      basismat2 <- eval.fd(x, fdobj2, Lfdobj2) * betamat2
    }
    iter <- 1
    tnm  <- 0.5
    # initialize array s
    chs <- width*crossprod(basismat1,basismat2)/2
    s[1,,] <- chs

    #  now iterate to convergence

    for (iter in 2:JMAX) {
      tnm <- tnm*2
      if (iter == 2) {
        x <- mean(rngi)
        nx <- 1
      } else {
        del <- width/tnm
        x   <- seq(rngi[1]+del/2, rngi[2]-del/2, del)
        nx  <- length(x)
      }
      #  first argument
      if (is.basis(basis1)) {
        if (type1) {
          betamat1  <- matrix(betafd1(x, bvec1, more1),nx,nbasis1)
        } else {
          betamat1  <- matrix(eval.fd(x, betafd1),nx,nbasis1)
        }
        basismat1 <- eval.basis(x, basis1, Lfdobj1) * betamat1
      } else {
        if (type1) {
          betamat1  <- matrix(betafd1(x, bvec1, more1),nx,nbasis1)
        } else {
          betamat1  <- matrix(eval.fd(x, betafd1),nx,nbasis1)
        }
        basismat1 <- eval.fd(x, basis1, Lfdobj1) * betamat1
      }
      #  second argument
      if (is.basis(fdobj2)) {
        if (type2) {
          betamat2 <- matrix(betafd2(x, bvec2, more2),nx,nbasis2)
        } else {
          betamat2 <- matrix(eval.fd(x, betafd2),nx,nbasis2)
        }
        basismat2 <- eval.basis(x, basis2, Lfdobj2) * betamat2
      } else {
        if (type2) {
          betamat2 <- matrix(betafd2(x, bvec2, more2),nx,nbasis2)
        } else {
          betamat2 <- matrix(eval.fd(x, betafd2),nx,ndim2)
        }
        basismat2 <- eval.fd(x, fdobj2, Lfdobj2) * betamat2
      }
      # update array s
      chs <- width*crossprod(basismat1,basismat2)/tnm
      chsold <- matrix(s[iter-1,,],dim(chs))
      s[iter,,] <- (chsold + chs)/2
      # predict next values on fifth or higher iteration
      if (iter >= 5) {
        ind <- (iter-4):iter
        ya <- array(s[ind,,],c(length(ind),ndim1,ndim2))
        xa <- h[ind]
        absxa <- abs(xa)
        absxamin <- min(absxa)
        ns <- min((1:length(absxa))[absxa == absxamin])
        cs <- ya
        ds <- ya
        y  <- matrix(ya[ns,,],ndim1,ndim2)
        ns <- ns - 1
        for (m in 1:4) {
          for (i in 1:(5-m)) {
            ho      <- xa[i]
            hp      <- xa[i+m]
            w       <- (cs[i+1,,] - ds[i,,])/(ho - hp)
            ds[i,,] <- hp*w
            cs[i,,] <- ho*w
          }
          if (2*ns < 5-m) {
            dy <- matrix(cs[ns+1,,],ndim1,ndim2)
          } else {
            dy <- matrix(ds[ns  ,,],ndim1,ndim2)
            ns <- ns - 1
          }
          y <- y + dy
        }
        ss     <- y
        errval <- max(abs(dy))
        ssqval <- max(abs(ss))
        # test for convergence
        if (all(ssqval > 0)) {
          crit <- errval/ssqval
        } else {
          crit <- errval
        }
        if (crit < EPS && iter >= JMIN) break
      }
      s[iter+1,,] <- s[iter,,]
      h[iter+1]   <- 0.25*h[iter]
      if (iter == JMAX) warning("Failure to converge.")
    }
  }

  return(ss)

}

#  -------------------------------------------------------------------------------

inprod.Dbasis.Data2LD <- function(fdobj1, fdobj2, coefList1, coefList2,
                                  Lfdobj1=int2Lfd(0), Lfdobj2=int2Lfd(0),
                                  EPS=1e-5, JMAX=15, JMIN=5)
{

  #  INPROD.DBASIS.DATA2LD  Computes the three-way tensor of inner products
  #  where the first two dimensions are the number of basis functions
  #  of fd functions in basis1 and basis2 respectively;
  #  and the third dimension is the partial derivatives of the first
  #  coefficient function with respect to its defining parameter vector.
  #  The  integration is approximated using Romberg integration with the
  #  trapezoidal rule.

  #  Arguments:
  #  BASIS1 and BASIS2:    These are funmctional basis objects.
  #  COEFLIST1 and COEFLIST2 ...  List objects for BASIS1 and BASIS2
  #  containing members:
  #               parvec   ... a vector of parameters
  #               estimate ... 0, held fixed, otherwise, estimated
  #               coeftype ... homogeneous or forcing
  #               fun      ... functional basis, fd, or fdPar object,
  #                            or a struct object for a general function
  #                            with fields:
  #                 fd      ... function handle for evaluating function
  #                 Dfd     ... function handle for evaluating
  #                              derivative with respect to parameter
  #                 more    ... object providing additional information for
  #                             evaluating coefficient function
  #  However, these may also be real constants that will be used as fixed coefficients.  An
  #  example is the use of 1 as the coefficient for the left side of the equation.
  #  Lfdobj1 and Lfdobj2:  order of derivatives for inner product for
  #               basis1 and basis2, respectively, or functional data
  #               objects defining linear differential operators
  #  EPS    convergence criterion for relative stop
  #  JMAX   maximum number of allowable iterations
  #  JMIN   minimum number of allowable iterations

  #  Return:
  #  A matrix of NREP1 by NREP2 of inner products for each possible pair
  #  of functions.

  #  Last modified 12 January 2018

  errwrd <- FALSE
  if (is.basis(fdobj1)) {
    basis1  <- fdobj1
    nbasis1 <- basis1$nbasis - length(basis1$dropind)
    ndim1   <- nbasis1
  } else {
    if (is.fd(fdobj1)) {
      basis1 <- fdobj1
      ndim1  <- 1
    } else {
      errwrd <- TRUE
      warning("First argument is not a basis object or an fd object.")
    }
  }

  if (is.basis(fdobj2)) {
    basis2  <- fdobj2
    nbasis2 <- basis2$nbasis - length(basis2$dropind)
    ndim2   <- nbasis2
  } else {
    if (is.fd(fdobj2)) {
      basis2  <- fdobj2
      ndim2   <- 1
    } else {
      errwrd <- TRUE
      warning("Second argument is not a basis object or an fd object.")
    }
  }

  if (errwrd) stop("Terminal error encountered.")

  bvec1 <- coefList1$parvec
  bvec2 <- coefList2$parvec

  #  Set up beta functions

  if (!is.basis(coefList1$fun) && !is.fd(coefList1$fun) && !is.fdPar(coefList1$fun)) {
    type1   <- TRUE
    betaDfd1 <- coefList1$fun$Dfd
    more1   <- coefList1$fun$more
  } else {
    type1   <- FALSE
    fdobj   <- coefList1$fun
    if (is.basis(fdobj)) {
      betabasis1 <- fdobj
    }
    if (is.fd(fdobj)) {
      betabasis1 <- fdobj$basis
    }
    if (is.fdPar(fdobj)) {
      betabasis1 <- fdobj$fd$basis
    }
  }

  if (!is.basis(coefList2$fun) && !is.fd(coefList2$fun) && !is.fdPar(coefList2$fun)) {
    type2   <- TRUE
    betafd2 <- coefList2$fun$fd
    more2   <- coefList2$fun$more
  } else {
    type2   <- FALSE
    fdobj   <- coefList2$fun
    if (is.fd(fdobj)) {
      betafd2 <- fd(coefList2$parvec, fdobj$basis)
    }
    if (is.basis(fdobj)) {
      betafd2 <- fd(coefList2$parvec, fdobj)
    }
    if (is.fdPar(fdobj)) {
      betafd2 <- fd(coefList2$parvec, fdobj$fd$basis)
    }
  }

  npar <- length(bvec1)

  #  check for any knot multiplicities in either argument

  knotmult <- numeric(0)
  if (type1 == "bspline") knotmult <- knotmultchk(basis1, knotmult)
  if (type2 == "bspline") knotmult <- knotmultchk(basis2, knotmult)

  #  Modify RNGVEC defining subinvervals if there are any
  #  knot multiplicities.

  if (is.basis(basis1)) {
    rng <- basis1$rangeval
  } else {
    rng <- basis1$basis$rangeval
  }
  if (length(knotmult) > 0) {
    knotmult <- sort(unique(knotmult))
    knotmult <- knotmult[knotmult > rng[1] && knotmult < rng[2]]
    rngvec   <- c(rng[1], knotmult, rng[2])
  } else {
    rngvec <- rng
  }

  #  -----------------------------------------------------------------
  #                   loop through sub-intervals
  #  -----------------------------------------------------------------

  nrng <- length(rngvec)
  for (irng  in  2:nrng) {
    rngi <- c(rngvec[irng-1],rngvec[irng])
    #  change range so as to avoid being exactly on
    #  multiple knot values
    if (irng > 2   ) rngi[1] <- rngi[1] + 1e-10
    if (irng < nrng) rngi[2] <- rngi[2] - 1e-10

    #  set up first iteration

    iter  <- 1
    width <- rngi[2] - rngi[1]
    JMAXP <- JMAX + 1
    h <- rep(1,JMAXP)
    h[2] <- 0.25
    tnm  <- 0.5
    iter <- 1
    s <- array(0,c(JMAX,ndim1,ndim2,npar))
    sdim <- length(dim(s))
    #  the first iteration uses just the endpoints
    x <- rngi

    #  For first argument:
    #  basis matrix evaluated with Lfdobj1
    if (is.basis(basis1)) {
      basismat1 <- eval.basis(x, basis1, Lfdobj1)
    } else {
      basismat1 <- eval.fd(x, basis1, Lfdobj1)
    }
    #  matrix of partial derivative values of first coefficient
    if (type1) {
      betamat1  <- betaDfd1(x, bvec1, more1)
    } else {
      betamat1  <- eval.basis(x, betabasis1, 1)
    }

    #  For second argument:
    #  basis matrix evaluated with Lfdobj1
    if (is.basis(basis2)) {
      basismat2 <- eval.basis(x, basis2, Lfdobj2)
    } else {
      basismat2 <- eval.fd(x, fdobj2, Lfdobj2)
    }
    #  vector of values of second coefficient
    if (type2) {
      betavec2 <- betafd2(x, bvec2, more2)
    } else {
      betavec2 <- eval.fd(x, betafd2)
    }

    temp2 <- basismat2*rep(betavec2,ndim2)
    for (k in 1:npar) {
      temp1 <- basismat1*rep(betamat1[,k],ndim1)
      chs <- width*crossprod(temp1,temp2)/2
      s[1,,,k] <- chs
    }

    #  now iterate to convergence

    for (iter in 2:JMAX) {
      tnm <- tnm*2
      if (iter == 2) {
        x <- mean(rngi)
      } else {
        del <- width/tnm
        x   <- seq(rngi[1]+del/2, rngi[2]-del/2, del)
      }

      #  For first argument:
      #  basis matrix evaluated with Lfdobj1
      if (is.basis(basis1)) {
        basismat1 <- as.matrix(eval.basis(x, basis1, Lfdobj1))
      } else {
        basismat1 <- as.matrix(eval.fd(x, fdobj1, Lfdobj1))
      }
      #  matrix of partial derivative values of first coefficient
      if (type1) {
        betamat1  <- betaDfd1(x, bvec1, more1)
      } else {
        betamat1  <- eval.basis(x, betabasis1, 1)
      }

      #  For second argument:
      #  basis matrix evaluated with Lfdobj1
      if (is.basis(basis2)) {
        basismat2 <- as.matrix(eval.basis(x, basis2, Lfdobj2))
      } else {
        basismat2 <- as.matrix(eval.fd(x, fdobj2, Lfdobj2))
      }
      #  vector of values of second coefficient
      if (type2) {
        betavec2 <- betafd2(x, bvec2, more2)
      } else {
        betavec2 <- eval.fd(x, betafd2)
      }

      temp2 <- basismat2*rep(betavec2,ndim2)
      for (k in 1:npar) {
        temp1 <- as.matrix(basismat1*betamat1[,k])
        chs <- width*crossprod(temp1,temp2)/tnm
        chsold <- s[iter-1,,,k]
        if (is.null(dim(chs))) chs <- matrix(chs,dim(chsold))
        s[iter,,,k] <- (chsold + chs)/2
      }

      if (iter >= 5) {
        ind <- (iter-4):iter
        ya <- array(s[ind,,,],c(length(ind),ndim1,ndim2,npar))
        xa <- h[ind]
        absxa <- abs(xa)
        absxamin <- min(absxa)
        ns <- min((1:length(absxa))[absxa == absxamin])
        cs <- ya
        ds <- ya
        y  <- array(ya[ns,,,],c(ndim1,ndim2,npar))
        ns <- ns - 1
        for (m in 1:4) {
          for (i in 1:(5-m)) {
            ho      <- xa[i]
            hp      <- xa[i+m]
            w       <- (cs[i+1,,,] - ds[i,,,])/(ho - hp)
            ds[i,,,] <- hp*w
            cs[i,,,] <- ho*w
          }
          if (2*ns < 5-m) {
            dy <- array(cs[ns+1,,,],c(ndim1,ndim2,npar))
          } else {
            dy <- array(ds[ns  ,,,],c(ndim1,ndim2,npar))
            ns <- ns - 1
          }
          y <- y + dy
        }
        ss     <- y
        errval <- max(abs(dy))
        ssqval <- max(abs(ss))
        if (all(ssqval > 0)) {
          crit <- errval/ssqval
        } else {
          crit <- errval
        }
        if (crit < EPS && iter >= JMIN) break
      }
      s[iter+1,,,] <- s[iter,,,]
      h[iter+1]    <- 0.25*h[iter]
      if (iter == JMAX) warning("Failure to converge.")
    }
  }

  return(ss)

}

#  -------------------------------------------------------------------------------

knotmultchk <- function(basisobj, knotmult) {
  type <- basisobj$type
  if (type == "bspline") {
    # Look for knot multiplicities in first basis
    params  <- basisobj$params
    nparams <- length(params)
    if (nparams > 1) {
      for (i in 2:nparams) {
        if (params[i] == params[i-1]) {
          knotmult <- c(knotmult, params[i])
        }
      }
    }
  }
  return(knotmult)
}


#  -------------------------------------------------------------------------------

DASarrayFn <- function(nXbasisw, nWbasisw, nUbasisj, nAbasisj, 
                       nrep,     Bvecw,    Ucoefj,   BAtens){
  # print("DASarrayFn")
  DASarray = .Call("DASarrayFnC", as.integer(nXbasisw), as.integer(nWbasisw), 
                   as.integer(nUbasisj), as.integer(nAbasisj), 
                   as.integer(nrep),     
                   as.double(Bvecw),    as.double(Ucoefj),   
                   as.double(BAtens))
  # print(length(DASarray))
  DASarray = array(DASarray,c(nXbasisw,nrep,nAbasisj))
  return(DASarray)
  
}

#  -------------------------------------------------------------------------------

DBSarrayFn <- function(nXbasisw, nWbasisw, nUbasisj, nAbasisj,  
                       nrep,     Avecj,    Ucoefj,   BAtens){
  # print("DBSarrayFn")
  DBSarray = .Call("DBSarrayFnC", as.integer(nXbasisw), as.integer(nWbasisw), 
                   as.integer(nUbasisj), as.integer(nAbasisj),  
                   as.integer(nrep),     
                   as.double(Avecj),    as.double(Ucoefj),   
                   as.double(BAtens))
  # print(length(DBSarray))
  DBSarray - array(DBSarray,c(nXbasisw,nrep,nAbasisj))
  return(DBSarray)
}

#  -------------------------------------------------------------------------------

DRarrayFn <- function(nXbasisw, nWbasisw, nXbasisx, nWbasisx, 
                      Bvecx,    Btens) {
  DRarray = .Call("DRarrayFnC", as.integer(nXbasisw), as.integer(nWbasisw), 
                  as.integer(nXbasisx), as.integer(nWbasisx), 
                  as.double(Bvecx),     as.double(Btens))
  DRarray <- array(DRarray, c(nXbasisw, nXbasisx, nWbasisw))
  return(DRarray)
}

#  -------------------------------------------------------------------------------

RmatFn <- function(nXbasisw, nWbasisw, nXbasisx, nWbasisx,
                   Bvecw,    Bvecx,    Btens) {
  # print("Rmat")
  Rmat <- .Call("RmatFnC", as.integer(nXbasisw), as.integer(nWbasisw), 
                as.integer(nXbasisx), as.integer(nWbasisx),
                as.double(Bvecw),     as.double(Bvecx),    
                as.double(Btens))
  Rmat <- matrix(Rmat,nXbasisw,nXbasisx)
  return(Rmat)
}

SmatFn <- function(nXbasisw, nWbasisw, nUbasisj, 
                   nAbasisj, nrep,     
                   Bvecw, Avecj, Ucoefj, BAtens){
  # print("Smat")
  Smat = .Call("SmatFnC", as.integer(nXbasisw), as.integer(nWbasisw), 
               as.integer(nUbasisj), as.integer(nAbasisj), 
               as.integer(nrep),
               as.double(Bvecw),     as.double(Avecj),     
               as.double(Ucoefj),    as.double(BAtens))
  Smat = matrix(Smat,nXbasisw,nrep)
  return(Smat)
}

#  -------------------------------------------------------------------------------

yListCheck <- function(yList, nvar) {
  
  #  Last modified 7 June 2018
  
  if (!is.list(yList)) {
    stop("YLIST is not a list vector.")
  }
  if (!is.vector(yList)) {
    if (nvar == 1) {
      ynames = names(yList)
      if (ynames[[1]] == "argvals" && ynames[[2]] == "y") {
        yList.tmp = yList
        yList = vector("list",1)
        yList[[1]] = yList.tmp
      } else {
        stop("yList is not a vector and does not have correct names.")
      }
    } else {
      stop("yList is not a vector.")
    }
  }
  errwrd <- FALSE
  
  nvec    <- rep(0,nvar)
  dataWrd <- rep(FALSE,nvar)
  ydim    <- rep(0,nvar)
  nobs    <- 0
  #  find number of replications for first non-empty cell
  for (ivar in 1:nvar) {
    if (is.list(yList[[ivar]])) {
      yListi <- yList[[ivar]]
      nrep <- dim(as.matrix(yListi$y))[2]
      break
    }
  }
  #  loop through variables
  for (ivar in 1:nvar) {
    if (is.list(yList[[ivar]]) && !is.null(yList[[ivar]]$y)) {
      dataWrd[ivar] <- TRUE
      yListi <- yList[[ivar]]
      if (is.null(yListi$argvals)) {
        warning(paste("ARGVALS is not a member for (YLIST[[", ivar,"]]."))
        errwrd <- TRUE
      }
      ni <- length(yListi$argvals)
      nvec[ivar] <- ni
      ydimi <- dim(as.matrix(yListi$y))
      if (length(ydimi) > 2) {
        warning(paste("More than two dimensions for (y in YLIST[[",
                      ivar,"]]."))
        errwrd <- TRUE
      } else {
        ydim[ivar] <- ydimi[1]
      }
      #  set up and check NREP
      nrepi <- ydimi[2]
      if (nrepi != nrep) {
        warning("Second dimensions of YList.y are not equal.")
        errwrd <- TRUE
      }
      nobs <- nobs + 1
      if (ni != ydimi[1]) {
        print(c(ni,ydimi[1]))
        warning(paste("Length of ARGVALS and first dimension of Y",
                      "are not equal."))
        errwrd <- TRUE
      }
    } else {
      dataWrd[ivar] <- FALSE
    }
  }
  
  if (nobs == 0) {
    warning("No variables have observations.")
    errwrd <- TRUE
  }
  
  if (errwrd) {
    stop("One or more terminal stop encountered in YLIST.")
  }
  
  return(list(nrep=nrep, nvec=nvec, dataWrd=dataWrd))
  
}

#  -------------------------------------------------------------------------------

inprod.Data2LD <- function(fdobj1, fdobj2=NULL, Lfdobj1=int2Lfd(0), Lfdobj2=int2Lfd(0),
                           rng = range1, wtfd = 0, returnMatrix=FALSE)
{
  
  #  computes matrix of inner products of functions by numerical
  #    integration using Romberg integration
  
  #  Arguments:
  #  FDOBJ1 and FDOBJ2    These may be either functional data or basis
  #               function objects.  In the latter case, a functional
  #               data object is created from a basis function object
  #               by using the identity matrix as the coefficient matrix.
  #               Both functional data objects must be univariate.
  #               If inner products for multivariate objects are needed,
  #               use a loop and call inprod(FDOBJ1[i],FDOBJ2[i]).
  #     If FDOBJ2 is not provided or is NULL, it defaults to a function
  #     having a constant basis and coefficient 1 for all replications.
  #     This permits the evaluation of simple integrals of functional data
  #     objects.
  #  LFDOBJ1 and LFDOBJ2  order of derivatives for inner product for
  #               FDOBJ1 and FDOBJ2, respectively, or functional data
  #               objects defining linear differential operators
  #  RNG    Limits of integration
  #  WTFD   A functional data object defining a weight
  #  JMAX   maximum number of allowable iterations
  #  EPS    convergence criterion for relative stop
  #  RETURNMATRIX ... If False, a matrix in sparse storage model can be returned
  #               from a call to function BsplineS.  See this function for
  #               enabling this option.
  
  #  Return:
  #  A matrix of NREP1 by NREP2 of inner products for each possible pair
  #  of functions.
  
  #  Last modified 21 November 2017 by Jim Ramsay
  
  #  Check FDOBJ1 and get no. replications and basis object
  
  result1   <- fdchk(fdobj1)
  nrep1     <- result1[[1]]
  fdobj1    <- result1[[2]]
  coef1     <- fdobj1$coefs
  basisobj1 <- fdobj1$basis
  type1     <- basisobj1$type
  range1    <- basisobj1$rangeval
  
  #  Default FDOBJ2 to a constant function, using a basis that matches
  #  that of FDOBJ1 if possible.
  
  if (is.null(fdobj2)) {
    tempfd    <- fdobj1
    tempbasis <- tempfd$basis
    temptype  <- tempbasis$type
    temprng   <- tempbasis$rangeval
    if (temptype == "bspline") {
      basis2 <- create.bspline.basis(temprng, 1, 1)
    } else {
      if (temptype == "fourier") basis2 <- create.fourier.basis(temprng, 1)
      else                       basis2 <- create.constant.basis(temprng)
    }
    fdobj2 <- fd(1,basis2)
  }
  
  #  Check FDOBJ2 and get no. replications and basis object
  
  result2   <- fdchk(fdobj2)
  nrep2     <- result2[[1]]
  fdobj2    <- result2[[2]]
  coef2     <- fdobj2$coefs
  basisobj2 <- fdobj2$basis
  type2     <- basisobj2$type
  range2    <- basisobj2$rangeval
  
  # check ranges
  
  if (rng[1] < range1[1] || rng[2] > range1[2]) stop(
    "Limits of integration are inadmissible.")
  
  #  Call B-spline version if
  #  [1] both functional data objects are univariate
  #  [2] both bases are B-splines
  #  (3) the two bases are identical
  #  (4) both differential operators are integers
  #  (5) there is no weight function
  #  (6) RNG is equal to the range of the two bases.
  
  if (is.fd(fdobj1)                    && 
      is.fd(fdobj2)                    &&
      type1 == "bspline"               && 
      type2 == "bspline"               &&
      is.eqbasis(basisobj1, basisobj2) &&
      is.integer(Lfdobj1)              && 
      is.integer(Lfdobj2)              &&
      length(basisobj1$dropind) == 0   &&
      length(basisobj1$dropind) == 0   &&
      wtfd == 0                        && all(rng == range1)) {
    
    inprodmat <- inprod.bspline(fdobj1, fdobj2,
                                Lfdobj1$nderiv, Lfdobj2$nderiv)
    return(inprodmat)
  }
  
  #  check LFDOBJ1 and LFDOBJ2
  
  Lfdobj1 <- int2Lfd(Lfdobj1)
  Lfdobj2 <- int2Lfd(Lfdobj2)
  
  #  Else proceed with the use of the Romberg integration.
  
  #  ------------------------------------------------------------
  #  Now determine the number of subintervals within which the
  #  numerical integration takes.  This is important if either
  #  basis is a B-spline basis and has multiple knots at a
  #  break point.
  #  ------------------------------------------------------------
  
  #  set iter
  
  iter <- 0
  
  # The default case, no multiplicities.
  
  rngvec <- rng
  
  #  check for any knot multiplicities in either argument
  
  knotmult <- numeric(0)
  if (type1 == "bspline") knotmult <- knotmultchk(basisobj1, knotmult)
  if (type2 == "bspline") knotmult <- knotmultchk(basisobj2, knotmult)
  
  #  Modify RNGVEC defining subinvervals if there are any
  #  knot multiplicities.
  
  if (length(knotmult) > 0) {
    knotmult <- sort(unique(knotmult))
    knotmult <- knotmult[knotmult > rng[1] && knotmult < rng[2]]
    rngvec   <- c(rng[1], knotmult, rng[2])
  }
  
  #  check for either coefficient array being zero
  
  if ((all(c(coef1) == 0) || all(c(coef2) == 0)))
    return(matrix(0,nrep1,nrep2))
  
  #  -----------------------------------------------------------------
  #                   loop through sub-intervals
  #  -----------------------------------------------------------------
  
  #  Set constants controlling convergence tests
  
  JMAX <- 15
  JMIN <-  5
  EPS  <- 1e-4
  
  inprodmat <- matrix(0,nrep1,nrep2)
  
  nrng <- length(rngvec)
  for (irng  in  2:nrng) {
    rngi <- c(rngvec[irng-1],rngvec[irng])
    #  change range so as to avoid being exactly on
    #  multiple knot values
    if (irng > 2   ) rngi[1] <- rngi[1] + 1e-10
    if (irng < nrng) rngi[2] <- rngi[2] - 1e-10
    
    #  set up first iteration
    
    iter  <- 1
    width <- rngi[2] - rngi[1]
    JMAXP <- JMAX + 1
    h <- rep(1,JMAXP)
    h[2] <- 0.25
    s <- array(0,c(JMAXP,nrep1,nrep2))
    sdim <- length(dim(s))
    #  the first iteration uses just the endpoints
    fx1 <- eval.fd(rngi, fdobj1, Lfdobj1, returnMatrix)
    fx2 <- eval.fd(rngi, fdobj2, Lfdobj2, returnMatrix)
    #  multiply by values of weight function if necessary
    if (!is.numeric(wtfd)) {
      wtd <- eval.fd(rngi, wtfd, 0, returnMatrix)
      fx2 <- matrix(wtd,dim(wtd)[1],dim(fx2)[2]) * fx2
    }
    s[1,,] <- width*matrix(crossprod(fx1,fx2),nrep1,nrep2)/2
    tnm  <- 0.5
    
    #  now iterate to convergence
    
    for (iter in 2:JMAX) {
      tnm <- tnm*2
      if (iter == 2) {
        x <- mean(rngi)
      } else {
        del <- width/tnm
        x   <- seq(rngi[1]+del/2, rngi[2]-del/2, del)
      }
      fx1 <- eval.fd(x, fdobj1, Lfdobj1, returnMatrix)
      fx2 <- eval.fd(x, fdobj2, Lfdobj2, returnMatrix)
      if (!is.numeric(wtfd)) {
        wtd <- eval.fd(wtfd, x, 0, returnMatrix)
        fx2 <- matrix(wtd,dim(wtd)[1],dim(fx2)[2]) * fx2
      }
      chs <- width*matrix(crossprod(fx1,fx2),nrep1,nrep2)/tnm
      s[iter,,] <- (s[iter-1,,] + chs)/2
      if (iter >= 5) {
        ind <- (iter-4):iter
        ya <- s[ind,,]
        ya <- array(ya,c(5,nrep1,nrep2))
        xa <- h[ind]
        absxa <- abs(xa)
        absxamin <- min(absxa)
        ns <- min((1:length(absxa))[absxa == absxamin])
        cs <- ya
        ds <- ya
        y  <- ya[ns,,]
        ns <- ns - 1
        for (m in 1:4) {
          for (i in 1:(5-m)) {
            ho      <- xa[i]
            hp      <- xa[i+m]
            w       <- (cs[i+1,,] - ds[i,,])/(ho - hp)
            ds[i,,] <- hp*w
            cs[i,,] <- ho*w
          }
          if (2*ns < 5-m) {
            dy <- cs[ns+1,,]
          } else {
            dy <- ds[ns,,]
            ns <- ns - 1
          }
          y <- y + dy
        }
        ss     <- y
        errval <- max(abs(dy))
        ssqval <- max(abs(ss))
        if (all(ssqval > 0)) {
          crit <- errval/ssqval
        } else {
          crit <- errval
        }
        if (crit < EPS && iter >= JMIN) break
      }
      s[iter+1,,] <- s[iter,,]
      h[iter+1]   <- 0.25*h[iter]
      if (iter == JMAX) warning("Failure to converge.")
    }
    inprodmat <- inprodmat + ss
    
  }
  
  if((!returnMatrix) && (length(dim(inprodmat)) == 2)) {
    #  coerce inprodmat to be nonsparse
    return(as.matrix(inprodmat))
  } else {
    #  allow inprodmat to be sparse if it already is
    return(inprodmat)
  }
  
}

#  -------------------------------------------------------------------------------

fdchk <- function(fdobj) {
  
  #  check the class of FDOBJ and extract coefficient matrix
  
  if (inherits(fdobj, "fd")) {
    coef  <- fdobj$coefs
  } else {
    if (inherits(fdobj, "basisfd")) {
      coef  <- diag(rep(1,fdobj$nbasis - length(fdobj$dropind)))
      fdobj <- fd(coef, fdobj)
    } else { 
      stop("FDOBJ is not an FD object.")
    }
  }
  
  #  extract the number of replications and basis object
  
  coefd <- dim(as.matrix(coef))
  if (length(coefd) > 2) stop("Functional data object must be univariate")
  nrep     <- coefd[2]
  basisobj <- fdobj$basis
  
  return(list(nrep, fdobj))
  
}


#  -------------------------------------------------------------------------------

is.eqbasis <- function(basisobj1, basisobj2) {
  
  #  tests to see of two basis objects are identical
  
  eqwrd <- TRUE
  
  #  test type
  
  if (basisobj1$type != basisobj2$type) {
    eqwrd <- FALSE
    return(eqwrd)
  }
  
  #  test range
  
  if (any(basisobj1$rangeval != basisobj2$rangeval)) {
    eqwrd <- FALSE
    return(eqwrd)
  }
  
  #  test nbasis
  
  if (basisobj1$nbasis != basisobj2$nbasis) {
    eqwrd <- FALSE
    return(eqwrd)
  }
  
  #  test params
  
  if (any(basisobj1$params != basisobj2$params)) {
    eqwrd <- FALSE
    return(eqwrd)
  }
  
  #  test dropind
  
  if (any(basisobj1$dropind != basisobj2$dropind)) {
    eqwrd <- FALSE
    return(eqwrd)
  }
  
  return(eqwrd)
  
}


