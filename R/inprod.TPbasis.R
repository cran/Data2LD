inprod.TPbasis <- function(basis1,  basis2,  basis3,  basis4,      
                           Lfdobj1=int2Lfd(0), Lfdobj2=int2Lfd(0), 
                           Lfdobj3=int2Lfd(0), Lfdobj4=int2Lfd(0), 
                           rng=c() , wtfd=0, EPS=1E-5, 
                           JMAX=16, JMIN=5) {
  #  INPROD.TPBASIS  Computes vectorized four-tensor of inner products of the
  #    tensor product of values of respective linear differential operators
  #    applied to four bases.
  #    The inner products are approximated by numerical integration using
  #    Romberg integration with the trapezoidal rule.
  #    The inner s can be over a reduced range in RNG and use a
  #    scalar weight function in WTFD.
  #
  #  Arguments:
  #  BASIS1, BASIS2, BASIS3 and BASIS4   these are basis objects,
  #            and the inner products of all quadruples of
  #            BASIS1, BASIS2, BASIS3 and BASIS4 are computed.
  #  LFDOBJ1, LFDOBJ2, LFDOBJ3 and LFDOBJ4  linear differential operators
  #                    for the corresponding basis functions.
  #  RNG    Limits of integration
  #  WTFD   A functional data object defining a weight
  #  EPS    A convergence criterion, defaults to 1e-4.
  #  JMAX   Maximum number of Richardson extrapolation iterations.
  #            Defaults to 16.
  #  JMIN   Minimum number of Richardson extrapolation iterations.
  #            Defaults to 5.
  #
  #  Return:
  #  An order NBASIS1*NBASIS2*NBASIS3*NBASIS4 matrix SS of inner products
  #  for each possible pair quadruple of basis functions.
  
  #  Last modified 2 January 2018
  
  #  check LFD objects
  
  Lfdobj1 <- int2Lfd(Lfdobj1)
  Lfdobj2 <- int2Lfd(Lfdobj2)
  Lfdobj3 <- int2Lfd(Lfdobj3)
  Lfdobj4 <- int2Lfd(Lfdobj4)
  
  #  check WTFD
  
  if (is.fd(wtfd)) {
    coefw <- wtfd$coef
    coefd <- dim(coefw)
    if (coefd[2] > 1) {
      stop('Argument WTFD is not a single function')
    } #
  } #
  
  #  check basis function objects
  
  if (!(is.basis(basis1) && is.basis(basis2) && 
        is.basis(basis3) && is.basis(basis4))) {
    stop ('The four first arguments are not basis objects.')
  } #
  
  #  determine NBASIS1 and NASIS2, and check for common range
  
  nbasis1 <- basis1$nbasis - length(basis1$dropind)
  nbasis2 <- basis2$nbasis - length(basis2$dropind)
  nbasis3 <- basis3$nbasis - length(basis3$dropind)
  nbasis4 <- basis4$nbasis - length(basis4$dropind)
  
  ncum  <- cumprod(c(nbasis1, nbasis2, nbasis3, nbasis4))
  nprod <- ncum[4]
  
  range1  <- basis1$rangeval
  range2  <- basis2$rangeval
  range3  <- basis3$rangeval
  range4  <- basis4$rangeval
  
  if (is.null(rng)) {
    rng <- range1    
  } #
  if (rng[1] < range1[1] || rng[2] > range1[2] || 
      rng[1] < range2[1] || rng[2] > range2[2] || 
      rng[1] < range3[1] || rng[2] > range3[2] || 
      rng[1] < range4[1] || rng[2] > range4[2]) {
    stop('Limits of integration are inadmissible.')
  } #
  
  #  set up first iteration using only boundary values
  
  iter  <- 0
  width <- rng[2] - rng[1]
  JMAXP <- JMAX + 1
  h     <- matrix(1,JMAXP,1)
  h[2]  <- 0.25
  s     <- matrix(0,JMAXP,nprod)
  if (!is.numeric(wtfd)) {
    wtvec <- eval.fd(wtfd, rng)
  } else {
    wtvec <- matrix(1,2,1)
  } #
  bmat1  <- eval.basis(rng, basis1, Lfdobj1)
  bmat2  <- eval.basis(rng, basis2, Lfdobj2)
  bmat3  <- eval.basis(rng, basis3, Lfdobj3)
  bmat4  <- eval.basis(rng, basis4, Lfdobj4)
  tensorprod <- InnerLoop(bmat1, bmat2, bmat3, bmat4, wtvec)
  chs    <- width * t(tensorprod) / 2
  s[1,]  <- chs
  tnm    <- 0.5
  
  #  now iterate to convergence
  
  for (iter  in  2:JMAX) {
    tnm <- tnm * 2
    del <- width / tnm
    x   <- rng[1] + (seq(del/2,rng[2],by <- del))
    nx  <- length(x)
    if ( !is.numeric(wtfd)) {
      wtvec <- eval.fd(wtfd, x)
    } else {
      wtvec <- matrix(1,nx,1)
    }
    bmat1  <- eval.basis(x, basis1, Lfdobj1)
    bmat2  <- eval.basis(x, basis2, Lfdobj2)
    bmat3  <- eval.basis(x, basis3, Lfdobj3)
    bmat4  <- eval.basis(x, basis4, Lfdobj4)
    tensorprod <- InnerLoop(bmat1, bmat2, bmat3, bmat4, wtvec)
    # print(max(abs(tensorprod)))
    chs <- width * t(tensorprod) / tnm
    chsold <- s[iter - 1,]
    s[iter,] <- (chsold + chs) / 2
    if (iter >= 5) {
      ind   <- (iter - 4):iter
      ya    <- matrix(s[ind,],length(ind),nprod)
      xa    <- h[ind]
      absxa <- abs(xa)
      absxamin <- min(absxa)
      ns <- which(absxa == absxamin)
      cs <- ya
      ds <- ya
      y  <- ya[ns,]
      ns <- ns - 1
      for (m  in  1:4) {
        for (i  in  1:(5 - m)) {
          ho      <- xa[i]
          hp      <- xa[i + m]
          w       <- (cs[i + 1,] - ds[i,])/(ho - hp)
          ds[i,] <- hp * w
          cs[i,] <- ho * w
        }
        if (2 * ns < 5 - m) {
          dy <- cs[ns + 1,]
        } else {
          dy <- ds[ns,]
          ns <- ns - 1
        }
        y <- y + dy
      }
      ss <- y
      errval <- max(abs(dy))
      ssqval <- max(abs(ss))
      if (all(ssqval > 10 * 1e-15)) {
        crit <- errval / ssqval
      } else {
        crit <- errval
      }
      if (crit < EPS && iter >= JMIN) {
        ss <- t(ss)
        return(ss)
      }
    }
    s[iter + 1,] <- s[iter,]
    h[iter + 1]  <- 0.25*h[iter]
  }
  print(paste('No convergence after ',JMAX,' steps in INPROD.BASIS.'))
  return(ss)
}
