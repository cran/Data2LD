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

fdchk <- function(basis) {
  
  #  check the class of basis and extract coefficient matrix
  
  if (is.fd(basis)) coef  <- basis$coefs
  else
    if (is.basis(basis)) {
      coef  <- diag(rep(1,basis$nbasis))
      basis <- fd(coef, basis)
    }
  else stop("basis is not an FD object.")
  
  #  extract the number of replications and basis object
  
  coefd <- dim(as.matrix(coef))
  if (length(coefd) > 2) stop("Functional data object must be univariate")
  nrep     <- coefd[2]
  basis <- basis$basis
  
  return(list(nrep, basis))
  
}

#  -------------------------------------------------------------------------------

knotmultchk <- function(basis, knotmult) {
  type <- basis$type
  if (type == "bspline") {
    # Look for knot multiplicities in first basis
    params  <- basis$params
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


