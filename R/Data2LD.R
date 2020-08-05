Data2LD <- function(yList, XbasisList, modelList, rhoVec=0.5*rep(1,nvar), 
                    summary=TRUE) {
#  Data2LD  stands for "Data to Linear Dynamics"
#  It approximates the data in argument YLIST by one or smooth
#  functions x_i, i=1,,d.  This approximation is defined by a set
#  of linear differential or algebraic equations defined by a set of
#  parameters some of which require estimation from the data.
#
#  The approximation minimizes the sum of squared residuals, expressed as
#  follows in Latex notation:
#
#    H(\theta) <- \sum_i^d \sum_j^n \sum_\ell^N [y_{ij \ell} - x_i(t_j)]^2
#
#  where:
#  i    <- 1,,d indexes equations in a system differential and algebraic
#                 equations.
#  j    <- 1,,n indexes times of observation of a variable
#  \ell <- 1,,N indexes replications of observations
#
#  But there is additional flexibility not captured in this expression:
#  1.  Only a subset of the variables may be observed, so that not all
#      values of index i are actually used.
#  2.  The number and location of times of observation t_j  can vary
#      from one observed variable to another.
#  using a roughness penaltylinear differential operator that depends
#  on unknown parameters in list MODELLIST, which is described below.
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
#      i=1,,d,  f=1,,F_i
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
#  ------------------------------------------------------------------------
#
#  Arguments:
#
#  YLIST      A list of length NVAR.  Each list contains in turn
#                a list object with fields:
#                  "argvals" is a vector of length n_i of observation times
#                  "y" contains the n_i observations.
#                The number of columns must be the same for all variables,
#                except that, if a list is empty, that variable is taken to
#                be not observed.
#
#  BASISLIST  A list array of length NVAR.  Each member contains in turn
#                a functional data object or a BASIS object.
#
#  MODELLIST  A list of length NVAR. Each list contains a
#                list object with members:
#                XList  list of length number of homogeneous terms
#                          Each list contains a list object with members:
#                          fun         a fdPar object for the coefficient
#                          variable    the index of the variable
#                          derivative  the order of its derivative
#                          ncoef       if coefficient estimated, its location
#                                         in the composite vector
#                          factor      a scalar multiplier (def. 1)
#                          estimate    0, held fixed, otherwise, estimated
#                FList  list of length number of forcing terms
#                          Each list contains a list object with members:
#                          AfdPar  an fdPar object for the coefficient
#                          Ufd     an fd object for the forcing function
#                          ncoef   if coefficient estimated, its location
#                                     in the composite vector
#                          factor  a scalar multiplier (def. 1)
#                          estimate    0, held fixed, otherwise, estimated
#                order  the highest order of derivative
#                name   a  tag for the variable
#                nallXterm  the number of homogeneous terms
#                nallFterm  the number of forcing functions
#
#  RHOVEC     A vector of length NVAR containing values in [0,1].
#                The data sums of squares are weighted by P and
#                the roughness penalty by 1-P.
#
#  ------------------------------------------------------------------------
#
#  Output objects (d <- number of equations,
#                  nparams is the total number of estimated parameters):
#
#  MSE      The weighted mean squared errors computed over the variables
#              with data.
#  DpMSE    The gradient of the objective function MSE with respect to
  #            the estimated parameters.
#  D2ppMSE  A square symmetric matrx of order nparams that contains
#              the second partial derivatives of the objective function.
#  XFDPARLIST  A list of length d containing functional parameter
#              objects of class fdPar for the estimated functions x_i(t).
#  DF       An equivalent degrees of freedom value
#                   df <- trace(2*YM - YM*YM) where YM is the matrix
#              fitMap described below.
#  GCV      The generalized cross-validation measure.  The value of
#              \rho corresponding to the minimim of GCV across values of
#              smoothing parameter \rho is often chose for an automatic
#              data-driven level of smoothing.
#  ISE      The sum across variables of the integrated squared value of
#              the differential operator.  This value multiplied by \rho
#              and divided by T, the width of the domain, is the second
#              term the objective function.
#  RMAT     A square symmetric matrix of order equal to the sum of the
#              coefficients in the basis function expansions of the
#              variables.  This matrix defines the size of the second term
#              in the objective function.
#  SMAT     A column vector of length equal to the order of RMAT.
#  FITMAP   A matrix with number of rows equal to the total number of
#              coefficients in the basis expansions of variables and
#              number of columns equal the total number of observations.
#              This matrix is the linear map from the data to the
#              combined coefficients.

#  Last modified 3 June 2020

nvar <- length(modelList)

#  ------------------------------------------------------------------------
#  Store the number of forcing functions for each variable and load the
#  starting position in the composite vector of estimated coefficients.
#  ------------------------------------------------------------------------

nhomog  <- rep(0,nvar)
nforce  <- rep(0,nvar)
nthetaH <- 0
nthetaF <- 0
nparams <- 0
for (ivar in 1:nvar) {
  modelListi <- modelList[[ivar]]
  nXterm <- modelListi$nallXterm
  #  process homogeneous terms
  if (nXterm > 0) {
    nhomog[ivar] <- nXterm
    for (iterm in 1:nXterm) {
      XListi    <- modelListi$XList[[iterm]]
      nparamsi <- length(XListi$parvec)
      nthetaH  <- nthetaH + nparamsi
      nparams  <- nparams + nparamsi
    }
  }
  nFterm <- modelListi$nallFterm
  #  process forcing terms
  if (nFterm > 0) {
    nforce[ivar] <- nFterm
    for (iterm in 1:nFterm) {
      FListi <- modelListi$FList[[iterm]]
      nparamsi <- length(FListi$parvec)
      nthetaF  <- nthetaF + nparamsi
      nparams  <- nparams + nparamsi
    }
  }
}

#  ------------------------------------------------------------------------
#                        Check YLIST
#  ------------------------------------------------------------------------

ycheckList <- yListCheck(yList, nvar)
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

#  check that cells contain basis objects

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
#  set up list for matrices of basis function values
#  ------------------------------------------------------------------------

basismatList <- list(nvar,1)
for (ivar in 1:nvar) {
    if (dataWrd[[ivar]]) {
        yListi  <- yList[[ivar]]
        basisi  <- XbasisList[[ivar]]
        argvals <- as.vector(yListi$argvals)
        basismati <- eval.basis(argvals, basisi)
        basismatList[[ivar]] <- basismati
    } else {
      basismatList[[ivar]] <- matrix(0, max(nvec), 1)
    }
}

#  ------------------------------------------------------------------------
#                  Compute coefficient matrix Bmat
#  ------------------------------------------------------------------------

nsum     <- sum(nvec)
ncoefsum <- sum(ncoefvec)
Bmat     <- matrix(0,ncoefsum,ncoefsum)
basismat <- matrix(0,nsum,ncoefsum)
ymat     <- matrix(0,nsum,1)
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
        Bmat[ind,ind] <- 
          weighti*(1-rhoVec[ivar])*crossprod(basismati)/nvec[ivar]
        basismat[indn,ind] <- basismati
    }
}

#  ------------------------------------------------------------------------
#      Compute roughness penalty matrices Rmat(theta) and Smat(theta)
#      and their derivatives with respect to estimated parameters
#  ------------------------------------------------------------------------

#  Matrices R and DR

#  There are parameters to estimate defining the homogeneous terms

Data2LDRList <- Data2LD.Rmat(XbasisList, modelList, rhoVec, nthetaH)
Rmat <- Data2LDRList$Rmat
if (nthetaH > 0) {
  DRarray <- Data2LDRList$DRarray
} else {
  #  No estimated parameters are involved in homogeneous terms
  DRarray <- NULL
}

#  Matrices S and DS for variables having forcing functions

Data2LDSList <- Data2LD.Smat(XbasisList, modelList, rhoVec,
                             nthetaH, nthetaF, nforce, nparams)
Smat <- Data2LDSList$Smat

if (nthetaF > 0) {
  DSmat <- Data2LDSList$DSmat
} else {
  #  No estimated parameters are involved in forcing terms
  DSmat <- NULL
}

Cmat <- Bmat + Rmat

#  Ceigvals <- eigen(Cmat)$values

#  ------------------------------------------------------------------------
#                     Set up right side of equation
#  ------------------------------------------------------------------------

Dmat <- matrix(0,ncoefsum,1)
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
coef = as.matrix(coef)

#  ------------------------------------------------------------------------
#  Compute the vector of unpenalized error sum of squares, MSE,
#  the sum of which is the outer objective function H(\theta|\rho).
#  Each MSE_i is normalized by dividing by the n_i"s.
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
    MSE    <- MSE + weighti*SSEi/nvec[ivar]
  }
}

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
  fitMap    <- basismat %*% RgtFactor

  #  Use fitMap to compute a equivalent degrees of freedom measure

  df <- sum(diag(2*fitMap)) - sum(fitMap^2)
   
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
  
  #  compute residual variance
  
  Rvar <- SSEtot/nsum
  
  #  compute GCV
  if (df < nsum) {
    gcv <- Rvar/((nsum - df)/nsum)^2
  } else {
    gcv <- NA
  }

  #  ------------------------------------------------------------------------
  #  Compute unpenalized error integrated squaresP
  #  ------------------------------------------------------------------------

  ISE <- Data2LD.ISE(XbasisList, modelList, coef, Rmat, Smat, nforce, rhoVec)
} else {
  df  <- NULL
  gcv <- NULL
  ISE <- NULL
  XfdParList <- NULL
  fitMap <- NULL
}

#  ------------------------------------------------------------------------
#       Compute total derivative of MSE wrt theta
#  ------------------------------------------------------------------------

#  Compute the partial derivatives of the coefficients with respect to the
#  estimated parameters,  dc/dtheta

Dcoef <- array(0, c(ncoefsum,nparams))

for (iparam in 1:nparams) {
  if (iparam <= nthetaH) {
    DRmati <- DRarray[,,iparam]
    if (is.null(DSmat)) {
      DRi <- -DRmati %*% coef
      Dcoef[,iparam] <- matrix(solve(Cmat,DRi),ncoefvec[ivar],1)
    } else {
      DRi <- -DRmati %*% coef
      DSi <- -DSmat[,iparam]
      Dcoef[,iparam] <- matrix(solve(Cmat,(DRi + DSi)),ncoefsum,1)
    }
  } else {
    DSi <- -DSmat[,iparam]
    Dcoef[,iparam] <- matrix(solve(Cmat, DSi),ncoefsum,1)
  }
}

#  ------------------------------------------------------------------------
#            Compute the total theta-gradient of H
#  ------------------------------------------------------------------------

xmat    <- basismat %*% coef
DpMSE   <- matrix(0,nparams, 1)
D2ppMSE <- matrix(0,nparams, nparams)
if (summary) {
  D2pyMSE <- matrix(0,nparams, nsum)
}

xmat    <- basismat %*% coef
DpMSE   <- matrix(0,nparams,1)
D2ppMSE <- matrix(0,nparams,nparams)
m2 <- 0
for (ivar in 1:nvar) {
  modelListi <- modelList[[ivar]]
  nXterm <- modelListi$nallXterm
  if (nXterm > 0) {
    for (iterm in 1:nXterm) {
      XListi   <- modelListi$XList[[iterm]]
      m1 <- m2 + 1
      m2 <- m2 + length(XListi$parvec)
      estimate <- as.numeric(XListi$estimate)
      for (m in m1:m2) {
        if (estimate[m - m1 + 1] == 0) {
          D2ppMSE[m,m] <- 1
        }
      }
    }
  }
  nFterm <- modelListi$nallFterm
  if (nFterm > 0) {
    for (iterm in 1:nFterm) {
      FListi   <- modelListi$FList[[iterm]]
      m1 <- m2 + 1
      m2 <- m2 + length(FListi$parvec)
      estimate <- as.numeric(FListi$estimate)
      for (m in m1:m2) {
        if (estimate[m - m1 + 1] == 0) {
          D2ppMSE[m,m] <- 1
        }
      }
    }
  }
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
    resmati    <- as.matrix(yVeci - xmati)
    BasDcoefi  <- basismati %*% Dcoef
    DpMSE      <- DpMSE   - 2*weighti*
      crossprod(BasDcoefi,resmati)/nvec[ivar]
    D2ppMSE    <- D2ppMSE + 
      2*weighti*crossprod(BasDcoefi)/nvec[ivar]
    if (summary) {
         D2pyMSE[,m1:m2] <- D2pyMSE[,m1:m2] - 
           2*weighti*t(BasDcoefi)/nvec[ivar]
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
              Rmat=Rmat, Smat=Smat, DRarray=DRarray, DSmat=DSmat,
              fitMap=fitMap))

} else {
  return(list(MSE=MSE, DpMSE=DpMSE, D2ppMSE=D2ppMSE, XfdParList=XfdParList,
              df=df, gcv=gcv, ISE=ISE))
}

}

#   ---------------------------------------------------------------------------

Data2LD.Rmat <- function(XbasisList, modelList,  rhoVec=rep(0.5,nvar), 
                         nthetaH) {
  #  Data2LD_Rmat computes the penalty matrix R associated with the homogeneous
  #  portion of a linear differential operator L as well as its partial
  #  derivative with respect to parameters defining the homogeneous portion.
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
  #  BASISLIST  A functional data object or a BASIS object.  If so, the
  #               smoothing parameter LAMBDA is set to 0.
  #
  #  MODELLIST  A List aray of length NVAR. Each List contains a
  #                list object with members:
  #                XList  list of length number of homogeneous terms
  #                          Each List contains a list object with members:
  #                          WfdPar  a fdPar object for the coefficient
  #                          variable    the index of the variable
  #                          derivative  the order of its derivative
  #                          npar        if coefficient estimated, its location
  #                                         in the composite vector
  #                          estimate    0, held fixed, otherwise, estimated
  #                FList  List arrau of length number of forcing terms
  #                          Each List contains a list object with members:
  #                          AfdPar    an fdPar object for the coefficient
  #                          Ufd       an fd object for the forcing function
  #                          npar      if coefficient estimated, its location
  #                                       in the composite vector
  #                          estimate  0, held fixed, otherwise, estimated
  #                order      the highest order of derivative
  #                name       a  tag for the variable
  #                nallXterm  the number of homogeneous terms
  #                nallFterm  the number of forcing functions
  #  RHOVEC   A vector of length NVAR containing values in [0,1].
  #                The data sums of squares are weighted by RHO and
  #               the roughness penalty by 1-RHO.
  
  #  Last modified 3 June 2020
  
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
  XListi <- list(parvec=1, estimate=FALSE, fun=fd(1,conbasis))
  
  #  ------------------------------------------------------------------------
  #                         Compute penalty matrix Rmat(theta)
  #  ------------------------------------------------------------------------
  
  #  The order of Rmat is the sum of the number of numbers of spline
  #  coefficients for each variable or trajectory
  
  Rmat <- matrix(0,ncoefsum,ncoefsum)
  
  #  loop through variables or equations
  
  m2 <- 0
  for (ivar in 1:nvar) {
    #  select list object for this variable that defines this variable's
    #  structure.
    modelListi <- modelList[[ivar]]
    #  weight for this variable for computing fitting criterion
    weighti    <- modelListi$weight
    #  indices in Rmat for this order nXbasisi submatrix
    nXbasisi  <- ncoefvec[ivar]
    m1  <- m2 + 1
    m2  <- m2 + ncoefvec[ivar]
    indi <- m1:m2
    #  order of the equation
    order     <- modelListi$order
    #  index of homogeneous terms
    nXtermi   <- modelListi$nallXterm
    #  the functional basis object
    Xbasisi   <- XbasisList[[ivar]]
    #  First we compute the term in Rmat for the left side of the equation
    #  otherwise reformat the vector form of the submatrix
    Btensii <- modelListi$Btens[[nXtermi+1]][[nXtermi+1]]
    if (any(is.na(Btensii))) stop("NAs in Btensii")
    Rmatii  <- .Call("RmatFnCpp", as.integer(nXbasisi), as.integer(1), 
                     as.integer(nXbasisi), as.integer(1),      
                     as.double(1.0),  as.double(1.0),  
                     as.double(Btensii))    
    Rmatii <- matrix(Rmatii, nXbasisi, nXbasisi)
    if (any(is.na(Rmatii))) stop("NAs in Rmatii")
    #  multiply by weight, value of smoothing parameter rho, and divide by
    #  duration of interval
    Rmatii <- weighti*rhoVec[ivar]*Rmatii/T
    #  initialize the left side submatrix withinin supermatrix Rmat
    Rmat[indi,indi] <- Rmat[indi,indi] + Rmatii
    #  now we compute the contribution to R mat for each pair of homogeneous
    #  terms in the equation
    if (nXtermi > 0) {
    for (iw in 1:nXtermi) {
      #  select homogeneous term iw
      modelListiw <- modelListi$XList[[iw]]
      TListw   <- getHomoTerm(modelListiw)
      indw     <- (ncoefcum[TListw$iv]+1):ncoefcum[TListw$iv+1]
      nXbasisw <- ncoefvec[TListw$iv]
      Xbasisw  <- XbasisList[[TListw$iv]]
      if (nXtermi > 0) {
      for (ix in 1:nXtermi) {
        modelListix <- modelListi$XList[[ix]]
        #  get the details for homogeneous coefficient ix
        TListx    <- getHomoTerm(modelListix)
        indx      <- (ncoefcum[TListx$iv]+1):ncoefcum[TListx$iv+1]
        nXbasisx  <- ncoefvec[TListx$iv]
        Xbasisx   <- XbasisList[[TListx$iv]]
        #  set up the inner products for the two basis systems
        if (TListw$funtype || TListx$funtype) {
          #  the user-supplied coefficient function case ...
          #  must be computed anew for each pair
          Rmatwx <- matrix(inprod.basis.Data2LD(
            Xbasisw, Xbasisx, modelListiw, modelListix,
            TListw$nderiv, TListx$nderiv), nXbasisw, nXbasisx)
        } else {
          #  the much faster case where both homogeneous coefficients
          #  are B-spline functional data objects
          Btenswx <- modelListi$Btens[[iw]][[ix]]
          #  reformat the vector of values
          Rmatwx  <- .Call("RmatFnCpp", as.integer(nXbasisw),     
                                        as.integer(TListw$nWbasis), 
                                        as.integer(nXbasisx),     
                                        as.integer(TListx$nWbasis),
                                        as.double(TListw$Bvec),  
                                        as.double(TListx$Bvec), 
                                        as.double(Btenswx))
          Rmatwx <- matrix(Rmatwx, nXbasisw, nXbasisx)
        }
        #  apply the scale factors
        Rmatwx <- weighti*TListw$factor*TListx$factor*rhoVec[ivar]*Rmatwx/T
        #  increment the submatrix
        Rmat[indw,indx] <- Rmat[indw,indx] + Rmatwx
      }
      }
      #  now we need to compute the submatrices for the product of
      #  the left side of the equation and each homogeneous term in
      #  the equation
      if (TListw$funtype) {
        #  user-supplied coefficient case
        Rmatiw <- matrix(inprod.basis.Data2LD(
          Xbasisi, Xbasisw, XListi, modelListiw,
          order, TListw$nderiv),
          nXbasisi,nXbasisw)
      } else {
        #  reformat the vector of previous computed values
        Btensiw <- modelListi$Btens[[nXtermi+1]][[iw]]
        # Rmatiw  <- RmatFn(nXbasisi, as.integer(1), nXbasisw, TListw$nWbasis, as.integer(1),  
        #                   TListw$Bvec, Btensiw)
        Rmatiw  <- .Call("RmatFnCpp", as.integer(nXbasisi),     
                                      as.integer(1), 
                                      as.integer(nXbasisw),     
                                      as.integer(TListw$nWbasis),
                                      as.double(1.0),  
                                      as.double(TListw$Bvec), 
                                      as.double(Btensiw))
        Rmatiw <- matrix(Rmatiw, nXbasisi, nXbasisw)
      }
      #  apply the scale factors
      Rmatiw  <- weighti*TListw$factor*rhoVec[ivar]*Rmatiw/T
      #  subtract the increment for both off-diagonal submatrices
      Rmat[indi,indw] <- Rmat[indi,indw] -   Rmatiw
      Rmat[indw,indi] <- Rmat[indw,indi] - t(Rmatiw)
    }
    }
  }
  
  #  ------------------------------------------------------------------------
  #  Compute partial derivatives with respect to the parameters in vector
  #  theta that are to be estimated if they are required
  #  This is a three-dimensional array DRarray with the final dimension being
  #  The verbose commentary above will not be continued.
  #  ------------------------------------------------------------------------
  
  DRarray <- array(0,c(ncoefsum,ncoefsum,nthetaH))
  
  m2 <- 0
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    weighti    <- modelListi$weight
    m1   <- m2 + 1
    m2   <- m2 + ncoefvec[ivar]
    indi <- m1:m2
    nXtermi  <- modelListi$nallXterm
    nXbasisi <- ncoefvec[ivar]
    Xbasisi  <- XbasisList[[ivar]]
    nderivi  <- modelListi$order
    #  select only active coefficients requiring estimation
    #  select all active coefficients
    #  loop through active variables within equation ivar
    #  whose coefficient require estimation
    if (nXtermi > 0) {
    for (iw in 1:nXtermi) {
      modelListiw <- modelListi$XList[[iw]]
      TListw  <- getHomoTerm(modelListiw)
      Bestimw <- TListw$estim
      if (any(Bestimw)) {
        #  define coefficient of estimated variable and
        #  it's derivative index
        indw     <- (ncoefcum[TListw$iv]+1):ncoefcum[TListw$iv+1]
        indthw   <- modelListiw$index
        nXbasisw <- ncoefvec[TListw$iv]
        Xbasisw  <- XbasisList[[TListw$iv]]
        #  loop through all active variables within equation ivar
        if (nXtermi > 0) {
        for (ix in 1:nXtermi) {
          modelListix <- modelListi$XList[[ix]]
          TListx  <- getHomoTerm(modelListix)
          #  define coefficient of active variable and
          #  it's derivative index
          indx      <- (ncoefcum[TListx$iv]+1):ncoefcum[TListx$iv+1]
          nXbasisx  <- ncoefvec[TListx$iv]
          Xbasisx   <- XbasisList[[TListx$iv]]
          #  get the tensor vector for this pair of coefficients
          #  and derivatives
          if (TListw$funtype || TListx$funtype) {
            #  user-coded case
            DRarraywx <- array(inprod.Dbasis.Data2LD(
              Xbasisw, Xbasisx,
              modelListiw, modelListix,
              TListw$nderiv, TListx$nderiv),
              c(nXbasisw,nXbasisx,length(indthw)))
          } else {
            #  fda object case
            Btenswx <- modelListi$Btens[[iw]][[ix]]
            # DRarrayFnCpp returns a vector
            DRarraywx <- .Call("DRarrayFnCpp", as.integer(nXbasisw), 
                                               as.integer(TListw$nWbasis),  
                                               as.integer(nXbasisx), 
                                               as.integer(TListx$nWbasis), 
                                               as.double(TListx$Bvec), 
                                               as.double(Btenswx))
            #  reformat vector into an array
            DRarraywx <- array(DRarraywx, c(nXbasisw, nXbasisx, TListw$nWbasis))
          }
          #  rescale the inner product
          DRarraywx <- weighti*TListw$factor*TListx$factor*rhoVec[ivar]*
                       DRarraywx/T
          #  increment DRarray
          if (iw == ix) {
            DRarray[indw,indw,indthw[Bestimw]] <-
            DRarray[indw,indw,indthw[Bestimw],drop=FALSE] + 
              2*DRarraywx[,,Bestimw,drop=FALSE]
          } else {
            DRarray[indw,indx,indthw[Bestimw]] <-
            DRarray[indw,indx,indthw[Bestimw],drop=FALSE] + 
              DRarraywx[,,Bestimw,drop=FALSE]
            DRarraywxt <- aperm(DRarraywx,c(2,1,3))
            DRarray[indx,indw,indthw[Bestimw]] <-
            DRarray[indx,indw,indthw[Bestimw],drop=FALSE] + 
              DRarraywxt[,,Bestimw,drop=FALSE]
          }
        }
        }
        #  partial derivatives wrt Wcoef for cross-products with D^m
        #  here x <- ivar, Wbasisx is the constant basis, and
        #  Bvecx <- 1
        #  get the tensor vector for this pair of coefficients
        #  and derivatives
        if (TListw$funtype) {
          DRarraywi <- array(inprod.Dbasis.Data2LD(
            Xbasisw,   Xbasisi,
            modelListiw, XListi,
            TListw$nderiv,    nderivi),
            c(nXbasisw,nXbasisi,length(indthw)))
        } else {
          Btenswi   <- modelListi$Btens[[iw]][[nXtermi+1]]
          DRarraywi <- .Call("DRarrayFnCpp", as.integer(nXbasisw), 
                                             as.integer(TListw$nWbasis),  
                                             as.integer(nXbasisi), 
                                             as.integer(1), 
                                             as.double(1.0), 
                                             as.double(Btenswi))
          DRarraywi <- array(DRarraywi, c(nXbasisw, nXbasisi, TListw$nWbasis))
        }
        #  rescale the inner product
        DRarraywi <- weighti*TListw$factor*rhoVec[ivar]*DRarraywi/T
        #  update DRarray
        DRarray[indw,indi,indthw[Bestimw]] <- 
        DRarray[indw,indi,indthw[Bestimw],drop=FALSE] - DRarraywi[,,Bestimw,drop=FALSE]
        DRarraywit <- aperm(DRarraywi,c(2,1,3))
        DRarray[indi,indw,indthw[Bestimw]] <-
        DRarray[indi,indw,indthw[Bestimw],drop=FALSE] - DRarraywit[,,Bestimw,drop=FALSE]
      }
    }
    }
  }
  return(list(Rmat=Rmat, DRarray=DRarray))
}

#  ----------------------------------------------------------------------------

Data2LD.Smat <- function(XbasisList, modelList, rhoVec=rep(0.5,nvar),
                         nthetaH, nthetaF, nforce, nparams) {
  #  Data2LD  stands for "Data to Linear Dynamics"
  #  Data2LD_S computes the penalty matrix S associated with the forcing
  #  portion of a linear differential operator L, as well as its partial
  #  derivatives with  respect to the parameter vector.
  #  For a single variable whose approximated in terms of an exansion in
  #  terms of a vector \phi of basis functions, S is
  #                 S <- \int [L \phi(t)] U' dt.
  #  S has dimensions K and 1, where K is the number of basis
  #  functions in the expansion of the variable, NFORCE is the number of
  #  forcing functions.  
  #  If multiple variables are involved, then S is a composite matrix
  #  constructed from inner products and cross-products of the basis
  #  function vectors associate with each variable.  It's dimension will be
  #  \sum K_i by NFORCE.
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
  #  Arguments:
  #
  #  BASISLIST  A functional data object or a BASIS object.  If so, the
  #               smoothing parameter LAMBDA is set to 0.
  #  MODELLIST   A list of length NVAR. Each list contains a
  #                list object with members:
  #                XList  list List of length number of homogeneous terms
  #                          Each list contains a list object with members:
  #                          WfdPar  a fdPar object for the coefficient
  #                          variable    the index of the variable
  #                          derivative  the order of its derivative
  #                          npar  if coefficient estimated, its location
  #                                   in the composite vector
  #                FList  list array of length number of forcing terms
  #                          Each list contains a list object with members:
  #                          AfdPar  an fdPar object for the coefficient
  #                          Ufd     an fd object for the forcing function
  #                          npar  if coefficient estimated, its location
  #                                   in the composite vector
  #                order      the highest order of derivative
  #                name       a  tag for the variable
  #                nallXterm  the number of homogeneous terms
  #                nallFterm  the number of forcing functions
  #  RHOVEC      A value in [0,1].  The data are weighted by P and the
  #               roughness penalty by 1-P.
  #
  #  Last modified 3 June 2020
  
  #  ------------------------------------------------------------------------
  #                         Set up analysis
  #  ------------------------------------------------------------------------
  
  #  compute number of variables
  
  nvar <- length(modelList)
  
  #  Set up a vector NCOEFVEC containing number of coefficients used
  #  for the expansion of each variable
  
  ncoefvec <- matrix(0,nvar,1)
  for (ivar in 1:nvar) {
    ncoefvec[ivar] <- XbasisList[[ivar]]$nbasis
  }
  ncoefcum <- cumsum(c(0,ncoefvec))
  
  #  get the width of the time domain
  
  Xrange <- XbasisList[[1]]$rangeval
  T      <- Xrange[2] - Xrange[1]
  
  XListi <- vector("list",1)
  XListi$parvec   <- 1
  XListi$estimate <- 0
  XListi$funobj   <- fd(1,create.constant.basis(Xrange))
  
  #--------------------------------------------------------------------------
  #         Compute the penalty vector  or matrix Smat(theta)
  #--------------------------------------------------------------------------
  
  #  If there are no forcing functions, return Smat and DSmat as empty
  
  ncoefsum <- sum(ncoefvec)
  if (sum(nforce)==0) {
    Smat <- NULL
    DSmat <- NULL
    return(list(Smat <- Smat, DSmat <- DSmat))
  }
  
  #  Loop through variables
  
  Smat <- matrix(0,sum(ncoefvec),1)
  m2 <- 0
  for (ivar in 1:nvar) {
    m1  <- m2 + 1
    m2  <- m2 + ncoefvec[ivar]
    ind <- m1:m2
    if (nforce[ivar] > 0) {
      modelListi <- modelList[[ivar]]
      weighti  <- modelListi$weight
      Xbasisi  <- XbasisList[[ivar]]
      order    <- modelListi$order
      nXbasisi <- ncoefvec[ivar]
      nXtermi  <- modelListi$nallXterm
      nFtermi  <- modelListi$nallFterm
      for (jforce in 1:nFtermi) {
        # For this forcing term select parameter vectors for 
        # coefficient, factor, known forcing function and whether
        # coefficient has a basis expansion (funtypej = FALSE) or is
        # user-defined (funtypej = TRUE)
        modelListij <- modelListi$FList[[jforce]]
        TListj      <- getForceTerm(modelListij)
        #  Crossproducts of homogeneous terms with forcing terms
        if (nXtermi > 0) {
          for (iw in 1:nXtermi) {
            # For this homogeneous term obtain variable basis,
            # coefficient parameter vector and type of coefficient
            modelListiw <- modelListi$XList[[iw]]
            TListw      <- getHomoTerm(modelListiw)
            indw        <- (ncoefcum[TListw$iv]+1):ncoefcum[TListw$iv+1]
            nXbasisw    <- ncoefvec[TListw$iv]
            Xbasisw     <- XbasisList[[TListw$iv]]
            WfdParw     <- TListw$coefnList$fun
            if (TListj$funtype || TListw$funtype) {
              # Either the homogeneous coefficient or the 
              # forcing coefficient or both are user-defined,
              # use numerical integration to get inner products
              Smatjw <- matrix(inprod.basis.Data2LD(Xbasisw, TListj$Ufd, 
                                                    modelListiw, modelListij,
                                                    TListw$nderiv, 0), nXbasisw)
            } else {
              # Both coefficients are basis function expansions,
              # use previously computed inner product values
              BAtenswj <- modelListi$BAtens[[iw]][[jforce]]
              Smatjw <- .Call("SmatFnCpp", as.integer(nXbasisw), 
                              as.integer(TListw$nWbasis), 
                              as.integer(TListj$nUbasis), 
                              as.integer(TListj$nAbasis), 
                              as.double(TListw$Bvec),     
                              as.double(TListj$Avec),     
                              as.double(TListj$Ucoef),    
                              as.double(BAtenswj))
              Smatjw <- matrix(Smatjw, nXbasisw,1)
            }
            #  rescale Smatjw
            Smatjw <- weighti*TListw$factor*TListj$factor*rhoVec[ivar]*Smatjw/T
            #  update Smat
            Smat[indw,1] <- Smat[indw,1,drop=FALSE] + Smatjw[,,drop=FALSE]
          }
        }
        #  Crossproducts of D^m with forcing terms
        if (TListj$funtype) {
          Smatji <- matrix(inprod.basis.Data2LD(Xbasisi, TListj$Ufd, 
                                                XListi,  modelListij, 
                                                order, 0),nXbasisi,1)
        } else {
          BAtensji <- modelListi$BAtens[[nXtermi+1]][[jforce]]
          Smatji <- .Call("SmatFnCpp", as.integer(nXbasisi), 
                          as.integer(1), 
                          as.integer(TListj$nUbasis), 
                          as.integer(TListj$nAbasis), 
                          as.double(1.0),     
                          as.double(TListj$Avec),     
                          as.double(TListj$Ucoef),    
                          as.double(BAtensji))
          Smatji <- matrix(Smatji, nXbasisi, 1)
        }
        # rescale Smatji
        Smatji <- weighti*TListj$factor*rhoVec[ivar]*Smatji/T
        #  update Smat
        Smat[indw,] <- Smat[indw,,drop=FALSE] - Smatji[,,drop=FALSE]
      }
    }
  }
  
  #  ------------------------------------------------------------------------
  #  Compute partial derivatives of Smat if required with respect to theta
  #  in parvec(1:nparamsL)
  #  ------------------------------------------------------------------------
  
  # DSmat is a matrix with second dimension containing theta
  DSmat <- matrix(0, ncoefsum, nparams)
  
  #  ---------------  Computation of DASmat -------------------
  #  partial derivatives of product of homogeneous terms
  #  and forcing terms with respect to all non-fixed forcing coefficients
  m2 <- 0
  for (ivar in 1:nvar) {
    m1 <- m2 + 1
    m2 <- m2 + ncoefvec[ivar]
    modelListi <- modelList[[ivar]]
    weighti <- modelListi$weight
    indi <- m1:m2
    nXbasisi <- ncoefvec[ivar]
    nXtermi <- modelListi$nallXterm
    nFtermi <- modelListi$nallFterm
    Xbasisi <- XbasisList[[ivar]]
    order <- modelListi$order
    if (nforce[ivar] > 0) {
      for (jforce in 1:nFtermi) {
        modelListij <- modelListi$FList[[jforce]]
        TListj  <- getForceTerm(modelListij)
        indthj  <- modelListij$index
        Aestimj <- TListj$estim
        if (any(Aestimj)) {
          #  This coefficient is estimated, get its details
          #  The index set for the parameters for this forcing
          #  coefficient 
          #  partial derivatives wrt forcing term coefficients for
          #  those forcing terms requiring estimation of their
          #  coefficients
          if (nXtermi > 0) {
            for (iw in 1:nXtermi) {
              modelListiw <- modelListi$XList[[iw]]
              TListw <- getHomoTerm(modelListiw)
              iv = TListw$iv
              indw   <- (ncoefcum[iv] + 1):ncoefcum[iv + 1]
              indthw <- modelListiw$index
              nXbasisw <- ncoefvec[iv]
              #  compute the matrix of products of
              #  partial derivatives of these forcing
              #  coefficients and homogeneous basis functions
              if (TListw$funtype || TListj$funtype) {
                # one or more user-defined coefficients
                DASmatjw <- inprod.Dbasis.Data2LD(TListj$Ufd, Xbasisw, 
                                                    modelListij, modelListiw,  
                                                    0, TListw$nderiv)
                DASmatjw <- aperm(DASmatjw,c(2,3,1))
              } else {
                BAtenswj   <- modelListi$BAtens[[iw]][[jforce]]
                #  DASmatFnCpp returns a vector
                DASmatjw <- .Call("DASarrayFnCpp", 
                                    as.integer(nXbasisw),  
                                    as.integer(TListw$nWbasis), 
                                    as.integer(TListj$nUbasis), 
                                    as.integer(TListj$nAbasis), 
                                    as.double(TListw$Bvec), 
                                    as.double(TListj$Ucoef), 
                                    as.double(BAtenswj))
                #  reformat the vector into a matrix
                DASmatjw <- matrix(DASmatjw, nXbasisw, length(indthj))
              }
              #  rescale DSmatjw
              DASmatjw <- weighti * TListj$factor * TListw$factor * 
                rhoVec[ivar] * DASmatjw/T
              DSmat[indw, indthj[Aestimj]] <- 
              DSmat[indw, indthj[Aestimj],drop=FALSE] + 
                DASmatjw[,Aestimj,drop=FALSE]
            }
          }
        }
        #  partial derivatives of cross-products of D^m
        #  with forcing terms wrt forcing coefficients
        if (TListj$funtype) {
          DASmatji <- inprod.Dbasis.Data2LD(TListj$Ufd, Xbasisi, 
                                              modelListij, XListi, 
                                              0, order)
        } else {
          BAtensij <- modelListi$BAtens[[nXtermi + 1]][[jforce]]
          DASmatji <- .Call("DASarrayFnCpp", 
                              as.integer(nXbasisi), 
                              as.integer(1),  
                              as.integer(TListj$nUbasis), 
                              as.integer(TListj$nAbasis), 
                              as.double(1),             
                              as.double(TListj$Ucoef), 
                              as.double(BAtensij))
        }
        DASmatji <- matrix(DASmatji, nXbasisi, length(indthj))
        #  rescale DSmatji
        DASmatji <- weighti * TListj$factor * rhoVec[ivar] * DASmatji/T
        DSmat[indi, indthj[Aestimj]] <- 
        DSmat[indi, indthj[Aestimj], drop=FALSE] - 
            DASmatji[,Aestimj,drop=FALSE]
      }
    }
  }
  
  #  ---------------  Computation of DBSmat -------------------
  #  partial derivatives of product of homogeneous terms and
  #  forcing terms with respect to all non-fixed homogeneous terms.
  
  for (ivar in 1:nvar) {
    modelListi = modelList[[ivar]]
    weighti = modelListi$weight
    nXtermi = modelListi$nallXterm
    nFtermi = modelListi$nallFterm
    #  Crossproducts of homogeneous terms with forcing terms,
    #  derivative with respect to homogeneous coefficient
    if (nXtermi > 0) {
      for (iw in 1:nXtermi) {
        modelListiw <- modelListi$XList[[iw]]
        TListw <- getHomoTerm(modelListiw)
        Bestimw <- TListw$estim
        if (any(Bestimw)) {
          ivw <- TListw$iv
          indw <- (ncoefcum[ivw] + 1):ncoefcum[ivw + 1]
          indthw <- modelListiw$index
          nXbasisw <- ncoefvec[ivw]
          Xbasisw <- XbasisList[[ivw]]
          if (nFtermi > 0) {
            for (jforce in 1:nFtermi) {
              modelListij <- modelListi$FList[[jforce]]
              TListj <- getForceTerm(modelListij)
              #  compute the matrix of products of 
              #  partial derivatives of these forcing 
              #  coefficients and homogeneous basis functions
              if (TListj$funtype || TListw$funtype) {
                DBSmatjw <- inprod.Dbasis.Data2LD(
                  Xbasisw, TListj$Ufd, modelListiw, 
                  modelListij, TListw$nderiv, 0)
              } else {
                BAtenswj <- modelListi$BAtens[[iw]][[jforce]]
                DBSmatjw <- .Call("DBSarrayFnCpp", as.integer(nXbasisw),           
                                    as.integer(TListw$nWbasis), 
                                    as.integer(TListj$nUbasis), 
                                    as.integer(TListj$nAbasis), 
                                    as.double(TListj$Avec), 
                                    as.double(TListj$Ucoef), 
                                    as.double(BAtenswj))
              }
              DBSmatjw <- matrix(DBSmatjw, nXbasisw, length(indthw))
              # rescale DSmatjw
              DBSmatjw <- weighti * TListj$factor * TListw$factor *  
                rhoVec[ivar] * DBSmatjw/T
              DSmat[indw, indthw[Bestimw]] <- 
              DSmat[indw, indthw[Bestimw],drop=FALSE] + 
                DBSmatjw[,,drop=FALSE]
            }
          }
        }
      }
    }
  }
  
  return(list(Smat=Smat, DSmat=DSmat))
  
}

#  ----------------------------------------------------------------------------

Data2LD.ISE <- function(XbasisList, modelList, coef, 
                        Rmat, Smat,  nforce, rhoVec=rep(0.5,nvar)) {
  #  Data2LD  stands for "Data to Linear Dynamics"
  #  Data2LD_ISE computes the value of the penalty term, the integrated
  #  squared difference between the right and left sides of a linear
  #  differential equation of the form
  #      D^m x_i <- sum_k^d sum_j^{m_k} beta_{kj} D^{j-1} x_k +
  #                sum_{\ell}^{M_i} \alpha_{\ell,i} u_{\ell,i},
  #      i=1,,d
  #  where
  #  and where each coefficient is expanded in terms of its own number of
  #  B-spline basis functions:
  #      \beta_{ij}(t)      <- \bbold_{ij}'     \phibold_{ij}(t),
  #      \alpha_{\ell,i}(t) <- \abold_{\ell,i}' \psibold_{\ell,i}(t)
  #
  #  This version approximates the integrals in the penalty terms by using
  #  inprod_basis to compute the cross-product matrices for the
  #  \beta-coefficient basis functions and the corresponding derivative of
  #  the x-basis functions,and the cross-product matrices for the
  #  \alpha-coefficients and the corresponding U functions.
  #  These are computed upon the first call to Data2LD, and then retained
  #  for subsequent calls by using the R-Cache command.
  #
  #  Arguments:
  #
  #  BASISLIST  A functional data object or a BASIS object.  If so, the
  #               smoothing parameter LAMBDA is set to 0.
  #
  #  MODELLIST   A list of length NVAR. Each list contains a
  #                list object with members:
  #                XList  list of length number of homogeneous terms
  #                          Each list contains a list object with members:
  #                          WfdPar  a fdPar object for the coefficient
  #                          variable    the index of the variable
  #                          derivative  the order of its derivative
  #                          npar  if coefficient estimated, its location
  #                                   in the composite vector
  #                FList  list of length number of forcing terms
  #                          Each list contains a list object with members:
  #                          AfdPar  an fdPar object for the coefficient
  #                          Ufd     an fd object for the forcing function
  #                          npar  if coefficient estimated, its location
  #                                   in the composite vector
  #                order      the highest order of derivative
  #                name       a  tag for the variable
  #                nallXterm  the number of homogeneous terms
  #                nallFterm  the number of forcing functions
  #
  #  RHOVEC     A vector of length NVAR containing values in [0,1].
  
  #                The data sums of squares are weighted by P and
  #                the roughness penalty by 1-rho.
  #
  #  COEF       The coefficient matrix
  #
  #  RMAT       The penalty matrix for the homogeneous part of L
  #
  #  SMAT       The penalty matrix for the forcing part of L
  #
  #  NFORCE     The vector containing the number of forcing functions
  #                per variable.
  
  #  Last modified 3 June 2020
  
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
  ISE1i <- 0
  ISE2i <- 0
  ISE3i <- 0
  for (ivar in 1:nvar) {
    ISE1i <- ISE1i + t(coef) %*% Rmat %*% coef
    modelListi <- modelList[[ivar]]
    weighti <- modelListi$weight
    nforcei  <- nforce[ivar]
    if (nforcei > 0) {
      ISE2i <- ISE2i + 2 %*% t(coef) %*% Smat
      ISE3i <- 0
      for (jforce in 1:nforcei) {
        modelListij <- modelListi$FList[[jforce]]
        TListj  <- getForceTerm(modelListij)
        for (kforce in 1:nforcei) {
          modelListik <- modelListi$FList[[kforce]]
          TListk    <- getForceTerm(modelListik)
          if (TListj$funtype || TListk$funtype) {
            ISE3i <- 
              inprod.basis.Data2LD(TListk$Ufd,  TListj$Ufd,
                                   modelListik, modelListij,  
                                   0,   0)
          } else {
            ISE3i <- 0
            ncum <- cumprod(c(TListk$nAbasis, TListk$nUbasis, 
                              TListj$nAbasis, TListj$nUbasis))
            Atensijk  <- modelListi$Atens[[jforce]][[kforce]]
            for (i in 1:TListj$nUbasis) {
              for (j in 1:TListj$nAbasis) {
                for (k in 1:TListk$nUbasis) {
                  for (l in 1:TListk$nAbasis) {
                    ijkl <- (i-1)*ncum[3] + (j-1)*ncum[2] + (k-1)*ncum[1] + l
                    ISE3i <- ISE3i +
                      TListj$Ucoef[i,1]*TListj$Avec[j]*
                      TListk$Ucoef[k,1]*TListk$Avec[l]*Atensijk[ijkl]
                  }
                }
              }
            }
          }
          ISE3i <- TListj$factor*TListk$factor*rhoVec[ivar]*ISE3i/T
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
  ISE <- (sum(ISE1) + sum(ISE2) + sum(ISE3))
  
  return(ISE)
  
}

#  ----------------------------------------------------------------------------

inprod.basis.Data2LD <- function(fdobj1, fdobj2, modelList1, modelList2,
                                 Lfdobj1=int2Lfd(0), Lfdobj2=int2Lfd(0),
                                 EPS=1e-6, JMAX=15, JMIN=5) {
  #  INPROD.BASIS.DATA2LD  Computes matrix of inner products of bases by numerical
  #    integration using Romberg integration with the trapezoidal rule.
  #  This version multiplies bases by respective coefficient functions
  #    betafn1 and betafn2, and is intended for use within function
  #    Data2LD.  It produces a matrix.

  #  Arguments:
  #  FDOBJ1 and FDOBJ2:    These are functional data objects.
  #  MODELLIST1 and MODELLIST2   List objects for BASIS1 and BASIS2
  #  containing members:
  #               funobj   functional basis, fd, or fdPar object,
  #                            or a list object for a general function
  #                            with fields:
  #                 fd       function handle for evaluating function
  #                 Dfd      function handle for evaluating
  #                              derivative with respect to parameter
  #                 more     object providing additional information for
  #                             evaluating coefficient function
  #               parvec    a vector of parameters
  #               index     position within parameter vector of parameters
  #               estimate  0, held fixed, otherwise, estimated
  #               for homogeneous terms:
  #                 variable   index of variable in the system
  #                 derivative order of derivative for the left side
  #                 factor     constant multiplier of the term
  #               for forcing terms:
  #                 Ufd        single function data object for the forcing function
  #                 factor     constant multiplier of the term
  #  However, these may also be real constants that will be used as fixed coefficients.  An
  #  example is the use of 1 as the coefficient for the left side of the equation.
  #  Lfdobj1 and Lfdobj2:  order of derivatives for inner product of
  #               fdobj1 and fdobj2, respectively, or functional data
  #               objects defining linear differential operators
  #  EPS    convergence criterion for relative stop
  #  JMAX   maximum number of allowable iterations
  #  JMIN   minimum number of allowable iterations

  #  Return:
  #  A matrix of inner products for each possible pair of functions.

  #  Last modified 3 June 2020

  #  Determine where fdobj1 and fdobj2 are basis or fd objects, define
  #  BASIS1 and BASIS2, and check for common range

  errwrd <- FALSE
  dimResults1 <- find.dim(fdobj1)
  ndim1  <- dimResults1$ndim
  basis1 <- dimResults1$basis
  errwrd <- dimResults1$errwrd
  dimResults2 <- find.dim(fdobj2)
  ndim2  <- dimResults2$ndim
  basis2 <- dimResults2$basis
  errwrd <- dimResults2$errwrd
  if (errwrd) stop("Terminal error encountered.")

  #  get coefficient vectors

  bvec1 <- modelList1$parvec
  bvec2 <- modelList2$parvec

  #  Set up beta functions

  funobj1 <- modelList1$funobj
  funobj2 <- modelList2$funobj
  
  if (!inherits(funobj1, c("basisfd", "fd", "fdPar"))) {
    userbeta1 <- TRUE
  } else {
    userbeta1  <- FALSE
  }
  
  if (!inherits(funobj2, c("basisfd", "fd", "fdPar"))) {
    userbeta2 <- TRUE
  } else {
    userbeta2  <- FALSE
  }
  
  #  check for any knot multiplicities in either argument

  knotmult <- numeric(0)
  if (basis1$type == "bspline") knotmult <- knotmultchk(basis1, knotmult)
  if (basis2$type == "bspline") knotmult <- knotmultchk(basis2, knotmult)

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

    JMAXP <- JMAX + 1
    s <- array(0,c(JMAXP,ndim1,ndim2))
    
    iter  <- 1
    width <- rngi[2] - rngi[1]
    h <- rep(1,JMAXP)
    h[2] <- 0.25
    x <- rngi
    nx <- 2
    #  first argument
    termmat1 <- make.termmat(x, fdobj1, modelList1, ndim1, userbeta1)   
    #  second argument
    termmat2 <- make.termmat(x, fdobj2, modelList2, ndim2, userbeta2)   
    tnm  <- 0.5
    # initialize array s
    chs <- width*crossprod(termmat1,termmat2)/2
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
      termmat1 <- make.termmat(x, fdobj1, modelList1, ndim1, userbeta1)   
      #  second argument
      termmat2 <- make.termmat(x, fdobj2, modelList2, ndim2, userbeta2)   
      # update array s
      chs <- width*crossprod(termmat1,termmat2)/tnm
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

inprod.Dbasis.Data2LD <- function(fdobj1, fdobj2, modelList1, modelList2,
                                  Lfdobj1=int2Lfd(0), Lfdobj2=int2Lfd(0),
                                  EPS=1e-5, JMAX=16, JMIN=5)
{

  #  INPROD.DBASIS.DATA2LD  Computes the three-way tensor of inner products
  #  where the first two dimensions are the number of basis functions
  #  of fd functions in basis1 and basis2 respectively
  #  and the third dimension is the partial derivatives of the first
  #  coefficient function with respect to its defining parameter vector.
  #  The  integration is approximated using Romberg integration with the
  #  trapezoidal rule.

  #  Arguments:
  #  FDOBJ1 and FDOBJ2:    These are functional data objects.
  #  MODELLIST1 and MODELLIST2   List objects for BASIS1 and BASIS2
  #  containing members:
  #               funobj   functional basis, fd, or fdPar object,
  #                            or a list object for a general function
  #                            with fields:
  #                 fd       function handle for evaluating function
  #                 Dfd      function handle for evaluating
  #                              derivative with respect to parameter
  #                 more     object providing additional information for
  #                             evaluating coefficient function
  #               parvec    a vector of parameters
  #               index     position within parameter vector of parameters
  #               estimate  0, held fixed, otherwise, estimated
  #               for homogeneous terms:
  #                 variable   index of variable in the system
  #                 derivative order of derivative for the left side
  #                 factor     constant multiplier of the term
  #               for forcing terms:
  #                 Ufd        single function data object for the forcing function
  #                 factor     constant multiplier of the term
  #  However, these may also be real constants that will be used as fixed coefficients.  An
  #  example is the use of 1 as the coefficient for the left side of the equation.
  #  Lfdobj1 and Lfdobj2:  order of derivatives for inner product of
  #               fdobj1 and fdobj2, respectively, or functional data
  #               objects defining linear differential operators
  #  EPS    convergence criterion for relative stop
  #  JMAX   maximum number of allowable iterations
  #  JMIN   minimum number of allowable iterations
  
  #  Return:
  #  A matrix of inner products for each possible pair of functions.
  
  #  Last modified 3 June 2020
  
  errwrd <- FALSE
  dimResults1 <- find.dim(fdobj1)
  ndim1  <- dimResults1$ndim
  basis1 <- dimResults1$basis
  errwrd <- dimResults1$errwrd
  dimResults2 <- find.dim(fdobj2)
  ndim2  <- dimResults2$ndim
  basis2 <- dimResults2$basis
  errwrd <- dimResults2$errwrd
  if (errwrd) stop("Terminal error encountered.")

  #  extract coefficient function coefficients
  
  bvec1 <- modelList1$parvec
  bvec2 <- modelList2$parvec

  #  Set up beta functions
  
  funobj1 <- modelList1$funobj
  funobj2 <- modelList2$funobj
  
  if (!inherits(funobj1, c("basisfd", "fd", "fdPar"))) {
    userbeta1 <- TRUE
  } else {
    userbeta1  <- FALSE
  }
  
  if (!inherits(funobj2, c("basisfd", "fd", "fdPar"))) {
    userbeta2 <- TRUE
  } else {
    userbeta2  <- FALSE
  }
  
  npar <- length(bvec1)

  #  check for any knot multiplicities in either argument

  knotmult <- numeric(0)
  if (basis1$type == "bspline") knotmult <- knotmultchk(basis1, knotmult)
  if (basis2$type == "bspline") knotmult <- knotmultchk(basis2, knotmult)

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

    JMAXP <- JMAX + 1
    s <- array(0,c(JMAX,ndim1,ndim2,npar))
    
    iter  <- 1
    width <- rngi[2] - rngi[1]
    h <- rep(1,JMAXP)
    h[2] <- 0.25
    tnm  <- 0.5
    iter <- 1
    sdim <- length(dim(s))
    #  the first iteration uses just the endpoints
    x <- rngi
    #  first argument
    Dtermarray1 <- make.Dtermarray(x, fdobj1, modelList1, ndim1, userbeta1)   
    #  second argument
    termmat2 <- make.termmat(x, fdobj2, modelList2, ndim2, userbeta2)   
    for (ipar in 1:npar) {
      chs <- width*crossprod(Dtermarray1[,,ipar],termmat2)/2
      s[1,,,ipar] <- chs
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
      nx <- length(x)
      #  first argument
      Dtermarray1 <- make.Dtermarray(x, fdobj1, modelList1, ndim1, userbeta1)   
      #  second argument
      termmat2 <- make.termmat(x, fdobj2, modelList2, ndim2, userbeta2)   
      for (ipar in 1:npar) {
        temp <- as.matrix(Dtermarray1[,,ipar])
        if (nx == 1) temp <- t(temp)
        chs <- width*crossprod(temp,termmat2)/tnm
        chsold <- s[iter-1,,,ipar]
        if (is.null(dim(chs))) chs <- matrix(chs,dim(chsold))
        s[iter,,,ipar] <- (chsold + chs)/2
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
        if (crit < EPS && iter >= JMIN) {
          return(ss)
        }
      }
      if (iter == JMAX) {
        warning("Failure to converge.")
        return(ss)
      }
      s[iter+1,,,] <- s[iter,,,]
      h[iter+1]    <- 0.25*h[iter]
    }
  }
}

#  -------------------------------------------------------------------------------

getHomoTerm <- function(XtermList) {
  #  get details for homogeneous term 
  
  index     <- XtermList$index       # positions in parameter vector 
  iv        <- XtermList$variable    # index of the variable in this term
  factor    <- XtermList$factor      # fixed scale factor
  nderiv    <- XtermList$derivative  # order of derivative in this term
  Bvec      <- XtermList$parvec      # parameter vector for forcing coefficient
  estim     <- XtermList$estimate    # parameter to be estimated? (TRUE/FALSE)
  nWbasis   <- length(Bvec)          # number of basis functions
  funtype   <- !(is.basis(XtermList$fun) ||
                   is.fd(   XtermList$fun) ||
                   is.fdPar(XtermList$fun))
  return(list(index=index, iv=iv, factor=factor, Bvec=Bvec, 
              nWbasis=nWbasis, estim=estim, nderiv=nderiv, funtype=funtype))
}

#  -------------------------------------------------------------------------------

getForceTerm <- function(FtermList) {
  #  get details for a forcing term
  #  last modified 16 April 2020
  index     <- FtermList$index     # index of coefficient object in coefList
  factor    <- FtermList$factor    # fixed scale factor 
  Ufd       <- FtermList$Ufd       # functional data object for forcing term
  Avec      <- FtermList$parvec    # parameter vector for forcing coefficient
  estim     <- FtermList$estimate  # parameter to be estimated? (TRUE/FALSE)
  Ubasis    <- Ufd$basis           # functional basis object for forcing function
  Ucoef     <- Ufd$coef            # coefficient vector for forcing function
  nUbasis   <- Ubasis$nbasis       # number of basis fucntions
  nAbasis   <- length(Avec)        # number of coefficients for B-spline coeff.
  funtype   <- !(is.basis(FtermList$fun) ||
                   is.fd(   FtermList$fun) ||
                   is.fdPar(FtermList$fun))
  #  return named list
  return(list(FtermList=FtermList, funtype=funtype, 
              Ufd=Ufd, Ubasis=Ubasis, Ucoef=Ucoef, nUbasis=nUbasis,
              estim=estim, Avec=Avec, nAbasis=nAbasis, factor=factor))
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

yListCheck <- function(yList, nvar) {
  
  #  Last modified 3 June 2020
  
  if (!is.list(yList)) {
    stop("YLIST is not a list vector.")
  }
  if (!is.vector(yList)) {
    if (nvar == 1) {
      ynames <- names(yList)
      if (ynames[[1]] == "argvals" && ynames[[2]] == "y") {
        yList.tmp <- yList
        yList <- vector("list",1)
        yList[[1]] <- yList.tmp
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
      nobs <- nobs + 1
      if (ni != ydimi[1]) {
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
  
  return(list(nvec=nvec, dataWrd=dataWrd))
  
}

#  -------------------------------------------------------------------------------

inprod.Data2LD <- function(fdobj1, fdobj2=NULL, 
                           Lfdobj1=int2Lfd(0), Lfdobj2=int2Lfd(0),
                           rng=range1, wtfd=0) {
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
  #  RETURNMATRIX  If False, a matrix in sparse storage model can be returned
  #               from a call to function BsplineS.  See this function for
  #               enabling this option.
  
  #  Return:
  #  A matrix of of inner products for each possible pair
  #  of functions.
  
  #  Last modified 3 June 2020 by Jim Ramsay
  
  #  Check FDOBJ1 and get no. replications and basis object
  
  result1   <- fdchk(fdobj1)
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
    return(0)
  
  #  -----------------------------------------------------------------
  #                   loop through sub-intervals
  #  -----------------------------------------------------------------
  
  #  Set constants controlling convergence tests
  
  JMAX <- 15
  JMIN <-  5
  EPS  <- 1e-4
  
  inprodmat <- 0
  
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
    s <- matrix(0,JMAXP,1)
    sdim <- length(dim(s))
    #  the first iteration uses just the endpoints
    fx1 <- eval.fd(rngi, fdobj1, Lfdobj1)
    fx2 <- eval.fd(rngi, fdobj2, Lfdobj2)
    #  multiply by values of weight function if necessary
    if (!is.numeric(wtfd)) {
      wtd <- eval.fd(rngi, wtfd, 0)
      fx2 <- matrix(wtd,dim(wtd)[1],dim(fx2)[2]) * fx2
    }
    s[1,1] <- width*crossprod(fx1,fx2)/2
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
      fx1 <- eval.fd(x, fdobj1, Lfdobj1)
      fx2 <- eval.fd(x, fdobj2, Lfdobj2)
      if (!is.numeric(wtfd)) {
        wtd <- eval.fd(wtfd, x, 0)
        fx2 <- matrix(wtd,dim(wtd)[1],dim(fx2)[2]) * fx2
      }
      chs <- width*matrix(crossprod(fx1,fx2),1)/tnm
      s[iter,1] <- (s[iter-1,1] + chs)/2
      if (iter >= 5) {
        ind <- (iter-4):iter
        ya <- s[ind,1]
        ya <- matrix(ya,5,1)
        xa <- h[ind]
        absxa <- abs(xa)
        absxamin <- min(absxa)
        ns <- min((1:length(absxa))[absxa == absxamin])
        cs <- ya
        ds <- ya
        y  <- ya[ns,1]
        ns <- ns - 1
        for (m in 1:4) {
          for (i in 1:(5-m)) {
            ho      <- xa[i]
            hp      <- xa[i+m]
            w       <- (cs[i+1,1] - ds[i,1])/(ho - hp)
            ds[i,,] <- hp*w
            cs[i,,] <- ho*w
          }
          if (2*ns < 5-m) {
            dy <- cs[ns+1,1]
          } else {
            dy <- ds[ns,1]
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
      s[iter+1,1] <- s[iter,1]
      h[iter+1]   <- 0.25*h[iter]
      if (iter == JMAX) warning("Failure to converge.")
    }
    inprodmat <- inprodmat + ss
    
  }
  
  if (length(dim(inprodmat)) == 2) {
    #  coerce inprodmat to be nonsparse
    return(as.matrix(inprodmat))
  } else {
    #  allow inprodmat to be sparse if it already is
    return(inprodmat)
  }
  
}

#  -------------------------------------------------------------------------------

make.termmat <- function(x, fdobj, termList, ndim, userbeta) {
  #  Matrix termmat is constructed by pointwise multiplication of a matrix
  #  of values of basis functions (or functions) and a matrix of values
  #  of a coefficient function.
  #  the point-wise multiplication requires that the values be formatted
  #  so that the two matrices have the same dimensions.
  funobj <- termList$funobj
  bvec   <- termList$parvec
  bvec   <- as.matrix(bvec)
  if (is.basis(fdobj)) {
    if (userbeta) {
      betamat <- funobj$fd(x, bvec, funobj$more) %*% matrix(1,1,ndim)
    } else {
      betafd   <- funobj
      if (is.basis(betafd)) betafd <- fd(bvec, betafd) 
      if (is.fd(betafd))    betafd <- betafd  
      if (is.fdPar(betafd)) betafd <- fd(bvec, betafd$fd$basis)
      betamat <- eval.fd(x, betafd)   %*% matrix(1,1,ndim)
    }
    basismat <- eval.basis(x, fdobj)
    termmat <- basismat*betamat
  } else {
    if (userbeta) {
      betamat <- funobj$fd(x, bvec, funobj$more) %*% matrix(1,1,ndim)
    } else {
      betafd   <- funobj
      if (is.basis(betafd)) betafd <- fd(bvec, betafd) 
      if (is.fd(betafd))    betafd <- betafd  
      if (is.fdPar(betafd)) betafd <- fd(bvec, betafd$fd$basis)
      betamat <- eval.fd(x, betafd)   %*% matrix(1,1,ndim)
    }
    basismat <- eval.fd(x, fdobj)
    termmat  <- basismat*betamat
  }
  return(termmat)
}

#  -------------------------------------------------------------------------------

make.Dtermarray <- function(x, fdobj, termList, ndim, userbeta) {
  #  Matrix termarray is constructed by pointwise multiplication of a matrix
  #  of values of basis functions (or functions) and a matrix of values
  #  of a coefficient function derivatives.
  #  The point-wise multiplication requires that the values be formatted
  #  so that the two matrices have the same dimensions.
  nx      <- length(x)
  npar    <- length(termList$parvec)
  eyenpar <- diag(rep(1,npar))
  funobj  <- termList$funobj
  bvec    <- termList$parvec
  if (is.basis(fdobj)) {
    if (userbeta) {
      Dbetamat <- funobj$Dfd(x, bvec, funobj$more)
    } else {
      betaDfd <- funobj
      if (is.basis(betaDfd)) betaDfd <- fd(eyenpar, betaDfd) 
      if (is.fd(betaDfd))    betaDfd <- betaDfd  
      if (is.fdPar(betaDfd)) betaDfd <- fd(eyenpar, betaDfd$fd$basis)
      Dbetamat <- eval.fd(x, betaDfd) %*% matrix(1,1,ndim)
    }
    basismat <- eval.basis(x, fdobj)
  } else {
    if (userbeta) {
      Dbetamat <- funobj$Dfd(x, bvec, funobj$more)
    } else {
      betaDfd <- funobj
      if (is.basis(betaDfd)) betaDfd <- fd(eyenpar, betaDfd) 
      if (is.fd(betaDfd))    betaDfd <- betaDfd  
      if (is.fdPar(betaDfd)) betaDfd <- fd(eyenpar, betaDfd$fd$basis)
      Dbetamat <- eval.fd(x, betaDfd) %*% matrix(1,1,ndim)
    }
    basismat <- eval.fd(x, fdobj)
  }
  Dtermarray <- array(0,c(nx,ndim,npar))
  for (ipar in 1:npar) {
    Dtermarray[,,ipar] <- basismat*matrix(Dbetamat[,ipar],nx,ndim)
  }
  return(Dtermarray)
}

#  -------------------------------------------------------------------------------

find.dim <- function(fdobj) {
  #  Determine ndim, the number of coefficient functions.
  #  ndim will normally be one.
  
  errwrd <- FALSE
  if (is.basis(fdobj)) {
    # if fdobj is a basis, ndim is the number of basis functions
    basis <- fdobj
    ndim  <- basis$nbasis
  } else {
    if (is.fd(fdobj)) {
      #  if fdobj is a  <- function(al data object, ndim is 1
      basis <- fdobj$basis
      ndim  <- 1
    } else {
      errwrd <- TRUE
      print('First argument is neither a basis object nor an fd object.')
    }
  }
  return(list(ndim=ndim, basis=basis, errwrd=errwrd))
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
  basisobj <- fdobj$basis
  
  return(list(fdobj))
  
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


