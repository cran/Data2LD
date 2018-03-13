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


