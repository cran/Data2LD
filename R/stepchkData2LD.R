stepchkData2LD = function(stepval, cvec, deltac, limwrd, ind, 
                          climit=50*matrix(rep(1,2),2,1) %*% matrix(1,1,npar), 
                          active=1:npar, dbgwrd=FALSE) {
  #STEPCHK checks the step along a line for producing parameters within the
  #  limits specified by BOT and TOP
  #  LIMWRD    Logical variable permitting detection that parameter
  #               was on the boundary two steps in a row
  
  #  Last modified 28 August 2017
  
  npar = length(deltac)
  
  bot = climit[1,]
  top = climit[2,]
  
  # if step is too small, return with flag ind = 1
  
  if (stepval < 1e-7) {
    ind = 1
    return(list(stepval=stepval,ind=ind,limwrd=limwrd))
  }
  
  #  ensure that step does not go beyond lower limit on parameters
  
  stepvali   = stepval*deltac
  if (any(stepvali[active] < bot[active]-cvec[active])) {
    index   = active[stepvali[active] < bot[active]-cvec[active] & 
                       deltac[active] != 0]
    stepvalnew = min((bot[index]-cvec[index])/deltac[index])
    if (dbgwrd) {
      print(paste('Lower limit reached, old step, new step: ', 
                  stepval, stepvalnew))
    }
    stepval = stepvalnew
    #  check whether lower limit has been reached twice in a row
    if (limwrd[1]) {
      ind = 1
      return(list(stepval=stepval,ind=ind,limwrd=limwrd))
    } else {
      limwrd[1] = TRUE
    }
  } else {
    limwrd[1] = FALSE
  }
  if (stepval < 1e-7) {
    ind = 1
    return(list(stepval=stepval,ind=ind,limwrd=limwrd))
  }
  
  #  ensure that step does not go beyond upper limit on parameters
  
  stepvali   = stepval*deltac
  if (any(stepvali[active] > top[active]-cvec[active])) {
    index   = active[stepvali[active] > top[active]-cvec[active] & 
                       deltac[active] != 0]
    stepvalnew = min((top[index]-cvec[index])/deltac[index])
    if (dbgwrd) {
      print(paste('Upper limit reached, old step, new step: ', 
                  stepval, stepvalnew))
    }
    stepval = stepvalnew
    #  check whether upper limit has been reached twice in a row
    if (limwrd[2]) {
      ind = 1
      return(list(stepval=stepval,ind=ind,limwrd=limwrd))
    } else {
      limwrd[2] = TRUE
    }
  } else {
    limwrd[2] = FALSE
  }
  
  return(list(stepval=stepval,ind=ind,limwrd=limwrd))
  
}