coefCheck <- function(coefList) {
  #  COEFCHECK checks the list of coefficient functions defining a
  #  possibly forced linear system of llinear differential equations.
  #
  #  COEFCHECK goes through the following steps:
  #
  #  1, each list member is checked for being a list object.
  #  If this test is not passed, the checking terminates.
  #
  #  2, the presence and class of required fields are checked.
  #  These are:
  #  Name       Class
  #  parvec     real vector
  #  estimate   numeric, and converted to logical
  #  coeftype   string,  {'alpha', 'beta', 'force', 'homo', 'homog',
  #                       'homogeneous'}
  #  fun        basis, fd, fdPar, list
  #  If this test is not passed, the checking terminates.
  #
  #  3. for each coefList, the index field indicating estimated parameters
  #  is set up, with the indices of the estimated homogeneous or beta
  #  coefficients coming first.
  #
  #  4. for any general coefficient functions, the fields and their objects
  #  are the checked.
  #  The required fields are:
  #  Name       Class
  #  fd         function
  #  difdip     function
  #  more       any object
  #  
  
  #  Last modified 2 January 2018
  
  #  -------------------------------  Step 1  -------------------------------
  
  #  Check that coefList is a List object and obtain number of coefficient
  #  functions
  
  if (!is.list(coefList)) {
    stop("coefList is not a List object.")
  }
  
  ncoef <- length(coefList)
  
  #  check that each List contains a list object
  
  errwrd <- FALSE
  for (icoef in 1:ncoef) {
    coefListi <- coefList[[icoef]]
    if (!is.list(coefListi)) {
      print(paste("List ", icoef, " does not contain a list object."))
      errwrd <- TRUE
    }
  }
  
  if (errwrd) {
    stop("Errors found, cannot continue.")
  }
  
  #  -------------------------------  Step 2  -------------------------------
  
  #  check that that required fields are present and that their contents
  #  are in the right class
  
  errwrd <- FALSE
  coefListNew <- vector("list",length(coefList))
  for (icoef in 1:ncoef) {
    coefListi <- coefList[[icoef]]
    
    #  check that all struct objects contain a field named "parvec"
    
    if (is.null(coefListi$parvec)) {
      warning(paste("List object for member ", icoef,  
                    " does not have a parvec field."))
      errwrd <- TRUE
    } else {
      parvec <- as.matrix(coefListi$parvec)
      if (!is.numeric(parvec)) {
        warning(paste("Field parvec for member ", icoef, " is not numeric."))
        errwrd <- TRUE
      } else {
        parvecsize <- dim(parvec)
        if (length(parvecsize) != 2) {
          warning(paste("Field parvec for member ", icoef,
                        " is not a matrix."))
          errwrd <- TRUE
        } else {
          if (!(parvecsize[1] == 1 || parvecsize[2] == 1)) {
            warning(paste("Field parvec for member ", icoef,
                          " is not a vector."))
            errwrd <- TRUE
          } else {
            coefListi$parvec <- parvec
          }
        }
      }
    }
    
    #  check that all struct objects contain a field named "estimate"
    
    if (is.null(coefListi$estimate)) {
      warning(paste("Listuct object for member ", icoef,
                    " does not have a estimate field."))
      errwrd <- TRUE
    } else {
      estimate <- coefListi$estimate
      if (!is.numeric(estimate) && !is.logical(estimate)) {
        warning(paste("Field estimate for member ", icoef,
                      " is neither numeric or logical."))
        errwrd <- TRUE
      } else {
        if (length(estimate) != 1) {
          warning(paste("Field estimate for member ", icoef,
                        " is not of length 1."))
          errwrd <- TRUE
        } else {
          coefListi$estimate <- estimate
        }
      }
    }
    
    #  check that all list objects contain a member named "coeftype"
    
    if (is.null(coefListi$coeftype)) {
      warning(paste("Listuct object for member ", icoef,
                    " does not have a coeftype field."))
      errwrd <- TRUE
    } else {
      coeftype <- coefListi$coeftype
      if (!is.character(coeftype)) {
        warning(paste("coeftype for member ", icoef,
                      " is not a string."))
        errwrd <- TRUE
      }
      if (!(coeftype == "beta"        ||
            coeftype == "homo"        ||
            coeftype == "homogeneous" ||
            coeftype == "alpha"       ||
            coeftype == "force"       ||
            coeftype == "forcing")) {
        warning(paste("Field coeftype for member ", icoef, " is not one of: ",
                      "alpha, beta, homo, homogeneous, force, forcing."))
        errwrd <- TRUE
      }
    }
    
    #  check that all struct objects contain a field named "fun"
    
    if (is.null(coefListi$fun)) {
      warning(paste("List object for member ", icoef, " does not have a fun field."))
      errwrd <- TRUE
    } else {
      fun <- coefListi$fun
      if (!(is.basis(fun) || is.fd(fun) || is.fdPar(fun) || is.list(fun))) {
        warning(paste("Field fun for member ", icoef,
                      " is not a class basis, fd, fdPar, or list"))
        errwrd <- TRUE
      }
    }
    
    coefListNew[[icoef]] <- coefListi
    
  }
  
  if (errwrd) stop("Errors found, cannot continue.")
  
  #  -------------------------------  Step 3  -------------------------------
  
  #  generate index fields sequentially
  
  ntheta <- 0
  theta  <- NULL
  m2 <- 0
  for (icoef in 1:ncoef) {
    coefListi <- coefList[[icoef]]
    if (coefListi$estimate &&
        (coefListi$coeftype == "beta"        ||
         coefListi$coeftype == "homo"        ||
         coefListi$coeftype == "homogeneous")) {
      parveci <- coefListi$parvec
      npari <- length(parveci)
      m1 <- m2 + 1
      m2 <- m2 + npari
      coefListi$index <- m1:m2
      ntheta <- ntheta + npari
      theta  <- c(theta, parveci)
      coefListNew[[icoef]] <- coefListi
    }
  }
  
  for (icoef in 1:ncoef) {
    coefListi <- coefList[[icoef]]
    if (coefListi$estimate &&
        (coefListi$coeftype == "alpha"        ||
         coefListi$coeftype == "force"        ||
         coefListi$coeftype == "forcing")) {
      nbasisi <- length(coefListi$parvec)
      m1 <- m2 + 1
      m2 <- m2 + nbasisi
      coefListi$index <- m1:m2
      ntheta <- ntheta + nbasisi
      theta  <- c(theta, parveci)
      coefListNew[[icoef]] <- coefListi
    }
  }
  
  #  -------------------------------  Step 4  -------------------------------
  
  #  if the fun field is a list object, check the object
  
  errwrd <- FALSE
  for (icoef in 1:ncoef) {
    coefListi <- coefList[[icoef]]
    if (!(is.fdPar(coefListi$fun) || is.fd(coefListi$fun) || is.basis(coefListi$fun))) {
      fdList <- coefListi$fun
      if (!is.null(fdList$fd)) {
        fun.fd <- fdList$fd
        if (!(class(fun.fd) =="function")) {
          warning(paste("Field fun$fd for member ", coefListi,
                        " does not contain a function object."))
          errwrd <- TRUE
        }      
      } else {
        warning(paste("Field fun for member ", icoef," does not have a field fd."))
        errwrd <- TRUE
      }
      if (!is.null(fdList$Dfd)) {
        fun.Dfd <- fdList$Dfd
        if (!(class(fun.Dfd) == "function")) {
          warning(paste("Field fun for member ", icoef,
                         " does not contain a function object."))
          errwrd <- TRUE
        }
      } else {
        warning(paste("Field fun for member ", coefListi,
                       " does not have a field Dfd."))
        errwrd <- TRUE
      }
    }
  }
  
  if (errwrd) {
    stop("Errors found, cannot continue.")
  }
  
  return(list(coefList=coefListNew, theta=theta, ntheta=ntheta))
  
}
  
  
  