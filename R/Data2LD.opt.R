Data2LD.opt <- function(yList, XbasisList, modelList, coefList, rhoMat,
                        convcrit=1e-4, iterlim=20, dbglev=1, 
                        active, parMap=diag(rep(1,npar))) {

  #  Data2LD.opt optimizes a parameter vector theta defining a
  #    linear differential operator object used to smooth a set of data.
  #
  #  Arguments:
  #  YCELL     an array containing values of curves
  #               If the array is a matrix, rows must correspond to argument
  #               values and columns to replications, and it will be assumed
  #               that there is only one variable per observation.
  #               If Y is a three-dimensional array, the first dimension
  #               corresponds to argument values, the second to replications,
  #               and the third to variables within replications.
  #               If Y is a vector, only one replicate and variable are
  #               assumed.
  #  BASISCELL  A functional data object or a BASIS object.  If so, the
  #               smoothing parameter LAMBDA is set to 0.
  #  MODELCELL  A cell aray of length NVAR. Each cell contains a
  #                struct object with members:
  #                Xcell  cell array of length number of homogeneous terms
  #                          Each cell contains a struct object with members:
  #                          WfdPar  a fdPar object for the coefficient
  #                          variable    the index of the variable
  #                          derivative  the order of its derivative
  #                          npar  if coefficient estimated, its location
  #                                   in the composite vector
  #                Fcell  cell arrau of length number of forcing terms
  #                          Each cell contains a struct object with members:
  #                          AfdPar  an fdPar object for the coefficient
  #                          Ufd     an fd object for the forcing function
  #                          npar  if coefficient estimated, its location
  #                                   in the composite vector
  #                order      the highest order of derivative
  #                name       a  tag for the variable
  #                nallXterm  the number of homogeneous terms
  #                nallFterm  the number of forcing functions
  #  COEFCELL  ... A cell array of length NCOEF.  Each cell contaions:
  #                fdPar ... an fdPar object that defines a coefficient
  #                          function \beta(t) or \alpha(t) that multiplies a
  #                          variable derivative value or a forcing function,
  #                          respectively, and also whether the coefficient
  #                          is to be estimated or held fixed.  Note that
  #                          a specified coefficient function may appear in
  #                          more than one place in a system of differential
  #                          equations, and can multiply either a variable
  #                          derivative value or a forcing function
  #                          simultaneously.  However, the same function
  #                          cannot be fixed in place and estimated in
  #                          another.
  #                index ... a vector of integers inciding the position(s) in
  #                          the composite parameter vector \theta
  #                          occupied by the coefficients of the functional
  #                          data object (fd object) in the fdPar field.
  #  RHOMAT      A matrix with the number of rows equal to the number of values of the 
  #               smoothing parameter per variable to be used, and number of columns
  #               equal to the number of variables.  
  #               Each entry in the matrix is a value \rho in the interval [0,1)  
  #               so that 0 <= \rho < 1.  
  #               For each variable ivar in the system of equations and each optimization
  #               iopt, the data are weighted by \rho(iopt,ivar) and the
  #               roughness penalty by 1-\rho(iopt,ivar).
  #               It is expected that the values of \rho within each column will be 
  #               in ascending order, and the estimated parameters for each row
  #               are passed along as initial values to be used for the optimization
  #               defined by the values of \rho in the next row.
  #  CONVCRIT    One convergence criterion, or a vector of two criteria.
  #               The first criterion is applied to the function change,
  #               and the second is applied to the gradient norm.
  #  ITERLIM   Maximum number of iterations allowed.
  #  DBGLEV    An integer controlling amount of output per iteration.
  #               Defaults to 1, which prints summary results for
  #               each iteration.
  #  ACTIVE    A vector of indices of the parameters to be optimized.  If the index
  #               of a parameter does not appear in this vector, it is left fixed
  #               at its initial value even if its specification in coefList is that
  #               it be estimated.  If the parameters are constrained by K linear
  #               constraints, ACTIVE contains integers from 1 to NPAR - K, where
  #               NPAR is the number of parameters in the model.  In other words,
  #               ACTIVE specifies the constrained parameters to be estimated.
  #  PARMAP    A rectangular matrix with number of rows equal to the
  #               number of parameters to be estimated as defined in
  #               BWTCELL, and number of columns equal to the number of
  #               parameters less the number of linear constraints on the
  #               estimated parameters.  The columns of PARMAP must be
  #               orthnormal so that t(PARMAP) %*% PARMAP is an identity matrix.
  #               t(PARMAP) %*% THETA maps unconstrained parameters and the
  #               corresponding gradient into   constrained parameter space.
  #               PARMAP %*% THETA  maps constrained parameters and the
  #               corresponding gradient into unconstrained parameter space.
  #               PARMAP will usually be set up using the full QR
  #               decomposition of a linear constraint coefficient matrix t(A)
  #               where the constraints are of the form A P <- B, A and B
  #               being known matrices.  An example of such a constraint
  #               that arises often is one where two estimated coefficients
  #               are constrained to be equal.  For example, if a vvariable
  #               X involved in an equation the form a(x - x.0), where x.0
  #               is a fixed set point or defined target level for variable
  #               X, then this would be set up as a.1 x + a.2 x.0, where
  #               coefficients a.1 and a.2 are constrained to be equal in
  #               magnitude but opposite in sign, or a.1 + a.2 <- 0.

  #  Returns:  A named list object with these fields:
  #  THETA.OPT    The optimal parameter values.
  #  BWTCELL.OPT  The corresponding optimal coefficients for the
  #                  homogeneous terms in the differential equation
  #  AWTCELL.OPT  The corresponding optimal coefficients for the
  #                  forcing terms in the differential equation

  #  Last modified 22 January 2019

  theta  <- modelList2Vec(modelList, coefList)
  ntheta <- length(theta)
  npar   <- length(theta)

  nvar <- length(yList)

  thetaCon <- t(parMap) %*% theta

  nparCon <- length(thetaCon)
  climit <- matrix(c(-1000,1000),2,1) %*% matrix(1,1,nparCon)
  active  <- 1:nparCon

  #  check rhoMat and get nopt
  
  if (!is.matrix(rhoMat)) {
    rhoMat = matrix(rhoMat,1,length(rhoMat))
  }
  rhodim <- dim(rhoMat)
  if (rhodim[2] != nvar) {
    stop(paste("The second dimension of RHOMAT is not equal to",
               "the number of variables"))
  }
  nopt   <- rhodim[1]

  #   Lists to contain results if more than a single set of rho's are involved

  thetastore   <- matrix(0,nopt,npar)
  dfstore      <- matrix(0,nopt,1)
  gcvstore     <- matrix(0,nopt,1)
  coefList.opt <- vector("list", nopt)

  coefList.opti <- coefList

  #  ------------------------------------------------------------------------
  #                  loop through rho vectors
  #  ------------------------------------------------------------------------

  for (iopt in 1:nopt) {

    rhoVeci <- rhoMat[iopt,]

    #  compute initial criterion value, gradient and hessian
    
    Data2LDList <- Data2LD(yList, XbasisList, modelList, coefList.opti, rhoVeci,
                             summary=FALSE)
    
    fvec    <- Data2LDList$MSE
    grad    <- Data2LDList$DpMSE
    gradCon <- t(parMap) %*% grad

    fvecsum <- sum(fvec)

    norm <- sqrt(mean(gradCon^2))

    #  evaluate the initial update vector for correcting the initial theta

    deltac <- -gradCon

    #  initialize iteration status arrays

    iternum <- 0
    status <- c(iternum, fvec, norm)
    if (dbglev >= 1) {
      cat("\nIter.    Criterion   Grad Length")
      cat("\n")
      # # cat(format(c(iternum,status[2:3]),digits=4))
      # cat(iternum)
      # cat("        ")
      # cat(format(status[2],nsmall=4,digits=4))
      # cat("      ")
      # cat(format(status[3],nsmall=4,digits=4))
      # cat("\n")
      cat(iternum)
      cat("        ")
      cat(round(status[2],6))
      cat("      ")
      cat(round(status[3],6))
    } else {
      cat(".")
    }
    cat("\n")
    iterhist <- matrix(0,iterlim+1,length(status))
    iterhist[1,]  <- status
    if (iterlim == 0) {
      return ()
    }

    #  -------  Begin main iterations  -----------

    MAXSTEPITER <- 5
    MAXSTEP     <- 1000
    trial       <- 0.1
    reset       <- 0
    linemat     <- matrix(0,3,5)
    thetaoldCon <- thetaCon
    fvecsumold  <- fvecsum
    gradoldCon  <- gradCon
    dbgwrd      <- dbglev >= 2

    #  ---------------  beginning of optimization loop  -----------

    for (iter in 1:iterlim) {
      iternum <- iternum + 1
      #  set logical parameters
      dblwrd <- c(FALSE,FALSE)
      limwrd <- c(FALSE,FALSE)
      stpwrd <- FALSE
      ind    <- 0
      ips    <- 0
      #  compute slope
      linemat[2,1] <- sum(deltac*gradoldCon)
      #  normalize search direction vector
      sdg          <- sqrt(sum(deltac^2))
      deltac       <- deltac/sdg
      linemat[2,1] <- linemat[2,1]/sdg
      # initialize line search vectors
      linemat[,1:4] <- matrix(c(0, linemat[2,1], fvecsum),3,1) %*% matrix(1,1,4)
      stepiter  <- 0
      if (dbglev >= 2) {
        cat(paste("                 ", stepiter, "  "))
        cat(format(round(t(linemat[,1]),6)))
        cat("\n")
      }
      #  return with error condition if initial slope is nonnegative
      if (linemat[2,1] >= 0) {
        if (dbglev >= 2) {
          cat("Initial slope ")
          cat(format(round(linemat[2,1])))
          cat("  is nonnegative.\n")
          cat("\n")
        }
        deltac        <- -gradCon
        linemat[2,1]  <- -sum(gradCon^2)
        sdg           <- sqrt(sum(deltac^2))
        deltac        <- deltac/sdg
        linemat[,1:4] <- matrix(c(0, linemat[2,1], fvecsum),3,1) %*% matrix(1,1,4)
        linestore     <- matrix(linemat[,1],3,1)
        if (dbglev > 2) {
          cat(paste("                 ", stepiter, "  "))
          cat(format(round(t(linemat[,1]),6)))
          cat("\n")
        }
      }
      #  return successfully if initial slope is very small
      if (linemat[2,1] >= -min(c(1e-3,convcrit))) {
        if (dbglev >= 2) {
          cat("Initial slope too small.\n")
          status <- c(iternum, fvecsum, norm)
          cat("\n")
        }
        coefListnew <- coefList
        break
      }
      #  first step set to trial
      if (is.nan(trial)) stop("is.nan(trial)")
      linemat[1,5] <- trial
      #  ------------  begin line search iteration loop  ----------
      thetanewCon <- thetaCon
      for (stepiter in 1:MAXSTEPITER) {
        #  check the step size and modify if limits exceeded
        stepval <- linemat[1,5]
        result <- stepchkData2LD(linemat[1,5], thetaCon, deltac, limwrd, ind,
                                 climit, active, dbgwrd)
        linemat[1,5] <- result[[1]]
        ind          <- result[[2]]
        limwrd       <- result[[3]]
        # break if limit hit twice in a row
        if (ind == 1) {
          thetanew   <- theta
          fvecsumnew <- fvecsum
          break
        }
        #  break if current step size too small
        if (linemat[1,5] <= 1e-7) {
          if (dbglev >= 2) {
            cat("\nStepsize too small:  ")
            cat(linemat[1,5])
            cat("\n")
          }
          thetanew   <- theta
          fvecsumnew <- fvecsum
          break
        }
        #  update parameter vector
        thetanewCon <- thetaCon + linemat[1,5]*deltac
        #  ---------  update function, gradient and hessian  -----------
        thetanew    <- parMap %*% thetanewCon
        coefListnew <- modelVec2List(thetanew, coefList)
        Data2LDList <- Data2LD(yList, XbasisList, modelList, coefListnew,
                                  rhoVeci, summary=FALSE)
        fvecnew       <- Data2LDList$MSE
        gradnew       <- Data2LDList$DpMSE
        hessmatnew    <- Data2LDList$D2ppMSE
        gradnewCon    <- t(parMap) %*% gradnew
        hessmatnewCon <- t(parMap) %*% hessmatnew %*% parMap
        fvecsumnew    <- sum(fvecnew)
        #  -------------------------------------------------------------
        linemat[3,5] <- fvecsumnew
        #  compute new directional derivative
        linemat[2,5] <- sum(deltac*gradnewCon)
        if (dbglev >= 2) {
          cat(paste("                 ", stepiter, "  "))
          cat(format(round(t(linemat[,5]),6)))
          cat("\n")
        }
        #  compute next line search step, also testing for convergence
        result  <- stepit(linemat, ips, dblwrd, MAXSTEP)
        linemat <- result[[1]]
        ips     <- result[[2]]
        ind     <- result[[3]]
        dblwrd  <- result[[4]]
        trial   <- linemat[1,5]
        #  ind == 0 implies convergence
        if (ind == 0 || ind == 5) {
          break
        }
      }
      #  ------------  end line search iteration loop  ----------
      thetaCon <- thetanewCon
      theta    <- thetanew
      fvecsum  <- fvecsumnew
      gradCon  <- gradnewCon
      #  test for function value made worse
      if (fvecsum > fvecsumold) {
        #  Function value worse  warn and terminate
        if (dbglev >= 2) {
          cat("\n")
          cat("Criterion increased: ")
          cat(format(round(c(fvecsumold, fvecsum),4)))
          cat("\n")
        }
        #  reset parameters and fit
        thetaCon <- thetaoldCon
        fvecsum  <- fvecsumold
        deltac   <- gradCon
        if (reset == 1) {
          #  This is the second time in a row that this
          #     has happened   quit
          if (dbglev >= 2) {
            cat("Reset twice, terminating.\n")
            cat("Convergence not attained.\n")
          }
          #  return current status of optimization
          break
        } else {
          reset <- 1
        }
      } else {
        #  function value has not increased,  check for convergence
        RMSgrad <- sqrt(mean(gradCon^2))
        if (length(convcrit) > 1) {
          convtest <- fvecsumold - fvecsum < convcrit[1] && RMSgrad < convcrit[2]
        } else {
          convtest <- fvecsumold - fvecsum < convcrit
        }
        if (convtest) {
          norm   <- sqrt(mean(gradCon^2))
          status <- c(iternum, fvecsum, norm)
          if (dbglev >= 1) {
            cat(iternum)
            cat("        ")
            cat(round(status[2],6))
            cat("      ")
            cat(round(status[3],6))
            cat("\n")
            cat("Convergence reached.\n")
          }
          #  return current status of optimization
          break
        }
        #  update old parameter vectors and fit structure
        thetaoldCon <- thetaCon
        fvecsumold  <- fvecsum
        gradoldCon  <- gradCon
        hessmatCon  <- t(parMap) %*% hessmatnew %*% parMap
        #  update the line search direction vector
        deltac      <- -solve(hessmatCon,gradCon)
        reset       <- 0
      }
      norm   <- sqrt(mean(gradCon^2))
      status <- c(iternum, fvecsum, norm)
      iterhist[iter+1,] <- status
      if (dbglev >= 1) {
        cat(iternum)
        cat("        ")
        cat(round(status[2],6))
        cat("      ")
        cat(round(status[3],6))
        cat("\n")
      }
      if (dbglev >= 1 && iter == iterlim) {
        cat(
          "\nMaximum iterations reached but convergence not attained.\n")
      }
    }

    #  ---------------  end of optimization loop  -----------

    cat("\n")

    #  evaluate the solution at the final optimal value of theta
    
    Data2LDList <- Data2LD(yList, XbasisList, modelList, coefListnew,
                           rhoVeci, summary=TRUE)
    
    #  store the optimal values for this level of the \rho's the 
    #  storage objects
    
    thetastore[iopt,]    <- theta
    dfstore[iopt]        <- Data2LDList$df
    gcvstore[iopt]       <- Data2LDList$gcv
    coefList.opt[[iopt]] <- coefListnew

  }

  #  -------------------  end of rho loop  ------------------

  return(list(thetastore=thetastore, dfstore=dfstore, gcvstore=gcvstore,
              coefList.opt=coefList.opt))

}

#  ----------------------------------------------------------------------------

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
    # if (dbgwrd) {
    #   Rprint(paste('Lower limit reached, old step, new step: ',
    #                stepval, stepvalnew))
    # }
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
    # if (dbgwrd) {
    #   Rprint(paste('Upper limit reached, old step, new step: ',
    #               stepval, stepvalnew))
    # }
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

#  ----------------------------------------------------------------------------

stepit <- function(linemat, ips, dblwrd, MAXSTEP) {
  #STEPIT computes next step size in line search algorithm
  #  Arguments:
  #  LINEMAT:  Row 1 contains step values
  #            Row 2 contains slope values
  #            Row 3 contains function values
  #  IPS:      If 1, previous slope was positive
  #  DBLWRD:   Vector of length 2:  dblwrd[1] T means step halved
  #                                 dblwrd[2] T means step doubled
  #  MAXSTEP:  maximum size of step

  #  Last modified 4 March 2014 by Jim Ramsay

  #  Wolfe condition 1
  test1.1 = linemat[3,5] <= linemat[3,1] + linemat[1,5]*linemat[2,1]/20
  #  Wolfe condition 2
  test1.2 = abs(linemat[2,5]) <= abs(linemat[2,1])/10
  # disp([test1.1, test1.2])
  test1 = test1.1 && test1.2
  # test1 = test1.2
  test2 = linemat[3,5] > linemat[3,1]
  test3 = linemat[2,5] > 0
  if ((test1 || !test3) && test2) {
    #  ************************************************************
    #  function is worse and either slope is satisfory or negative
    ips = 0        #  step is halved
    if (dblwrd[2]) {
      ind = 5
      return(list(linemat = linemat, ips = ips, ind = ind, dblwrd = dblwrd))
    }
    linemat[1,5] = min(c(linemat[1,5]/2, MAXSTEP))
    linemat[,2] = linemat[,1]
    linemat[,3] = linemat[,1]
    dblwrd = c(1, 0)
    ind = 2
    return(list(linemat = linemat, ips = ips, ind = ind, dblwrd = dblwrd))
  }
  #  *********************************************************
  if (test1) {
    #  test1 means successful convergence
    ind = 0
    return(list(linemat = linemat, ips = ips, ind = ind, dblwrd = dblwrd))
  }
  #  **********************************************************
  if (test3) {
    #  Current slope is positive
    ips = 1
    linemat[,4] = linemat[,5]
    deltaf = linemat[3,3] - linemat[3,5]
    z = (3/(linemat[1,5] - linemat[1,3]))*deltaf + linemat[2,3] + linemat[2,5]
    w = z * z - linemat[2,3] * linemat[2,5]
    if (abs(linemat[2,3] + linemat[2,5] + 2 * z) >= 1e-05 && w > 0) {
      if (w < 0) stop("w < 0")
      w = sqrt(w)
      denom = linemat[2,5] - linemat[2,3] + 2 * w
      if (denom == 0) stop("denom == 0")
      linemat[1,5] = linemat[1,3] + (1 - ((linemat[2,5] + w - z)/
            (denom))) * (linemat[1,5] - linemat[1,3])
      if (is.nan(linemat[1,5])) stop("first test3: is.nan(linemat[1,5])")
    } else {
      #  linear interpolation necessary
      aerror = linemat[1,3]
      if (linemat[1,5] > linemat[1,3]) {
        aerror = linemat[1,5]
      }
      linemat[1,5] = linemat[1,3] - linemat[2,3] * ((linemat[1,5] - linemat[1,3])/
           (linemat[2,5] - linemat[2,3]))
      if (linemat[1,5] > 2 * aerror) {
        linemat[1,5] = 2 * aerror
      }
    }
    linemat[1,5] = min(c(linemat[1,5], MAXSTEP))
    dblwrd = c(0,0)
    ind = 2
    return(list(linemat = linemat, ips = ips, ind = ind, dblwrd = dblwrd))
  }
  #  *************************************************************
  #  Current slope is negative or zero
  linemat[,2] = linemat[,3]
  linemat[,3] = linemat[,5]
  if (ips == 1) {
    #  *****************************************************
    #  previous slope is positive
    deltaf = linemat[3,5] - linemat[3,4]
    z = c(3/(linemat[1,4] - linemat[1,5])) * deltaf +
      linemat[2,5] + linemat[2,4]
    w = z * z - linemat[2,5] * linemat[2,4]
    if (abs(linemat[2,5] + linemat[2,4] + 2 * z) >= 1e-05 && w > 0) {
      w = sqrt(w)
      linemat[1,5] = linemat[1,5] + (1 - ((linemat[2,4] + w - z)/
                                            (linemat[2,4] - linemat[2,5] +
                                               2 * w))) * (linemat[1,4] - linemat[1,5])
    } else {
      aerror = linemat[1,5]
      if (linemat[1,4] > linemat[1,5]) {
        aerror = linemat[1,4]
      }
      linemat[1,5] = linemat[1,5] - linemat[2,5] *
        ((linemat[1,4] - linemat[1,5])/
           (linemat[2,4] - linemat[2,5]))
      if (linemat[1,5] > 2 * aerror) {
        linemat[1,5] = 2 * aerror
      }
    }
    linemat[1,5] = min(c(linemat[1,5], MAXSTEP))
    dblwrd = c(0,0)
    ind = 2
    return(list(linemat = linemat, ips = ips, ind = ind, dblwrd = dblwrd))
  }
  #  ******************************************************
  if ((linemat[2,3] - linemat[2,2]) * (linemat[1,3] - linemat[1,2]) > 0) {
    #  previous slope is negative
    z = c(3/(linemat[1,3] - linemat[1,2])) * (linemat[3,2] - linemat[3,3]) +
      linemat[2,2] + linemat[2,3]
    w = z * z - linemat[2,2] * linemat[2,3]
    if (abs(linemat[2,2] + linemat[2,3] + 2 * z) >= 1e-05 && w > 0) {
      w = sqrt(w)
      linemat[1,5] = linemat[1,2] +
        (1 - ((linemat[2,3] + w - z)/(linemat[2,3] - linemat[2,2] +
                                        2 * w))) * (linemat[1,3] - linemat[1,2])
    } else {
      linemat[1,5] = linemat[1,2] - linemat[2,2] *
        ((linemat[1,3] - linemat[1,2])/
           (linemat[2,3] - linemat[2,2]))
    }
    linemat[1,5] = min(c(linemat[1,5], MAXSTEP))
    dblwrd = c(0,0)
    ind = 2
    return(list(linemat = linemat, ips = ips, ind = ind, dblwrd = dblwrd))
  } else {
    #  previous slope also negative but not as much
    if (dblwrd[1]) {
      ind = 5
      return(
        list(linemat = linemat, ips = ips, ind = ind, dblwrd = dblwrd))
    } else {
      linemat[1,5] = 2 * linemat[1,5]
      linemat[1,5] = min(c(linemat[1,5], MAXSTEP))
      dblwrd = c(0,1)
      ind = 2
      return(
        list(linemat = linemat, ips = ips, ind = ind, dblwrd = dblwrd))
    }
  }
  ind = 2

}

