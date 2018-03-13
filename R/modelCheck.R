modelCheck <- function(modelList, coefList) {
  # Check the modelList structure, requiring at a minimum:
  #   listarray with a member for each variable, containing:
  #     a list array of length equal number of homogeneous terms present.  
  #     Each list has:
  #       variable index i and/or tag
  #       order of derivative in the  term
  #       fdPar object for coefficient function
  #   list array of length either number of forcing functions or empty.  
  #     If not empty, lists have:
  #       fd object for forcing function(s)
  #       fdPar object for coefficient function
  
  #  Last modified 2 January 2018
  
  #  number of coefficients and number of variables
  
  ncoef <- length(coefList)
  nvar  <- length(modelList)
  
  if (is.null(modelList)) {
    stop('Argument MODELLIST is empty.')
  }
  
  if (!is.list(modelList)) {
    stop('Argument MODELLIST is not a list object.')
  }
  
  #  --------------------------------------------------------------------------
  #  check that each list is a list object and that it has the necessary
  #  fields or, if not, that these are set to default values
  #  --------------------------------------------------------------------------
  
  errwrd <- FALSE
  for (ivar in 1:nvar) {
    #  check modelList[[ivar]] is not empty and that it is a list object
    modelListi <- modelList[[ivar]]
    if(is.null(modelListi)) {
      warning(paste('No list object in modelList[[',ivar,']]'))
      errwrd <- TRUE
    }
    if (!is.list(modelListi)) {
      warning(paste('Object in modelList[[',ivar,']] is not a list object.'))
      errwrd <- TRUE
    }
  }
  
  if (errwrd) stop("One or more variables not defined.")
  
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    #  check modelList[[ivar]] for required fields, and assign default if needed
    #  check member "XList"
    #  check member "order"
    if (is.null(modelListi$order)) {
      #  assign default value
      modelListi$order <- 1
    }
    #  check member "name"
    if (is.null(modelListi$name)) {
      #  construct a name 
      modelListi$name <- paste("x",ivar,sep="")
    }
    #  check member "weight"
    if (is.null(modelListi$weight)) {
      #  construct unit weight 
      modelListi$weight <- 1
    }
    #  update containing list
    modelList[[ivar]] <- modelListi
  }
  
  #  terminate of an error has been identified
  
  if (errwrd) {
    stop('One or more terminal errors encountered.')
  }

  #  -------------------  check the values of each field  --------------------- 
  
  errwrd <- FALSE
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    #  check value of XList
    if (!is.null(modelListi$XList) && !is.list(modelListi$XList)) {
      warning(paste('XList field in modelList[[',ivar, 
                    ']] does not contain a list object.'))
      errwrd <- TRUE
    } 
    #  check value of FList
    if (!is.null(modelListi$FList) && !is.list(modelListi$FList)) {
      warning(paste('FList field in modelList[[',ivar, 
                    ']] does not contain a list object.'))
      errwrd <- TRUE
    }
    #  check value of order
    if (!is.numeric(modelListi$order)) {
      warning(paste('order field in modelList[[',ivar, 
                    ']] is not a number.'))
      errwrd <- TRUE
    }
    #  check value of name
    if (!is.character(modelListi$name)) {
      warning(paste('name field in modelList[[',ivar, 
                    ']] is not a character string.'))
      errwrd <- TRUE
    }
    #  check value of nallXterm
    if (is.null(modelListi$nallXterm)) {
      if (is.null(modelListi$XList)) {
        modelListi$nallXterm <- 0
      } else {
        modelListi$nallXterm <- length(modelListi$XList)
      }
    } else {
      if (!is.numeric(modelListi$nallXterm)) {
        warning(paste('nallXterm field in modelList[[',ivar, 
                      ']] is not a number.'))
        errwrd <- TRUE
      }
    }
    #  check value of nallFterm
    if (is.null(modelListi$nallFterm)) {
      if (is.null(modelListi$FList)) {
        modelListi$nallFterm <- 0
      } else {
        modelListi$nallFterm <- length(modelListi$FList)
      }
    } else {
      if (!is.numeric(modelListi$nallFterm)) {
        warning(paste('nallFterm field in modelList[[',ivar, 
                      ']] is not a number.'))
        errwrd <- TRUE
      }
    }
    modelList[[ivar]] <- modelListi
  }
    
  #  terminate of an error has been identified
  
  if (errwrd) {
    stop('One or more terminal errors encountered.')
  }
  
  #  -----------------------  check members of XList  -------------------------
  
  errwrd <- FALSE
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    if (!is.null(modelListi$XList[[1]])) {
      nXterm <- length(modelListi$XList)
      for (iterm in 1:nXterm) {
        XListij <- modelListi$XList[[iterm]]
        if (!is.null(XListij) && !is.list(XListij)) {
          warning(paste('XList[[',iterm,']] in modelList[[',ivar, 
                            ']] is not a list object.'))
          errwrd <- TRUE
        } else {
          #  check member ncoef
          if (is.null(XListij$ncoef)) {
            warning(paste('List object in XList[[',iterm, ']] in modelList[[', 
                          ivar,']] does not have a member ncoef.'))
            errwrd <- TRUE
          } else {
            ncoef <- XListij$ncoef
            if (!is.numeric(ncoef)) {
              warning(paste('XList[[',iterm, ']]$ncoef in modelList[[', 
                            ivar,']] is not a number.'))
              errwrd <- TRUE
            } else {
              if (round(ncoef) != ncoef) {
                warning(paste('XList[[',iterm, ']]$ncoef in modelList[[', 
                              ivar,']] is not an integer.'))
                errwrd <- TRUE
              }
            }
          }
          #  check member variable   
          if (is.null(XListij$variable)) {
            warning(paste('List object in XList[[',iterm,']] in modelList[[', 
                                ivar,']] does not have a member variable.'))
            errwrd <- TRUE
          } else {
            varij <- XListij$variable
            if (!is.numeric(varij)) {
              warning(paste('XList[[',iterm, ']]$variable in modelList[[', 
                            ivar,']] is not a number.'))
              errwrd <- TRUE
            } else {
              if (round(varij) != varij) {
                warning(paste('XList[[',iterm, ']]$ncoef in modelList[[', 
                              ivar,']] is not an integer.'))
                errwrd <- TRUE
              } else {
                if (varij > nvar) {
                  warning(paste('List object in XList[[',iterm, 
                                ']] in modelList[[',ivar, 
                                ']] has variable field > NVAR.'))
                  errwrd <- TRUE
                }
              }
            }
          }
          #  check member derivative
          if (is.null(XListij$derivative)) {
            warning(paste('List object in XList[[',iterm, 
                          ']] in modelList[[',ivar, 
                          ']] does not have a member derivative.'))
            errwrd <- TRUE
          } else {
            derivative <- XListij$derivative
            if (!is.numeric(derivative)) {
              warning(paste('XList[[',iterm, ']]$derivative in modelList[[', 
                            ivar,']] is not a number.'))
              errwrd <- TRUE
            } else {
              if (round(derivative) != derivative) {
                warning(paste('XList[[',iterm, ']]$derivative in modelList[[', 
                              ivar,']] is not an integer.'))
                errwrd <- TRUE
              } else {
                if (XListij$derivative >= modelListi$order) {
                  warning(paste('List object in XList[[',iterm, 
                                ']] in modelList[[',ivar, 
                                ']] has derivative member >= ORDER.'))
                  errwrd <- TRUE
                }
              }
            }
          }
        }
      }
    }
  }
   
  #  ----------------------  check members of FList  --------------------------
  
  errwrd <- FALSE
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    if (!is.null(modelListi$FList)) {
      FList <- modelListi$FList
      nFterm <- length(FList)
      for (jterm in 1:nFterm) {
        FListij <- FList[[jterm]]
        if (!is.null(FListij) && !is.list(FListij)) {
          warning(paste('FList[[',jterm,']] in modelList[[',ivar, 
                        ']] is not a list object.'))
          errwrd <- TRUE
        } else {
          if (is.null(FListij$ncoef)) {
            warning(paste('List object in FList[[',jterm,
                          ']] in modelList[[',ivar,
                          ']] does not have a ncoef field.'))
            errwrd <- TRUE
          } else {
            ncoef <- FListij$ncoef
            if (!is.numeric(ncoef)) {
              warning(paste('FList[[',iterm, ']]$ncoef in modelList[[', 
                            ivar,']] is not a number.'))
              errwrd <- TRUE
            } else {
              if (round(ncoef) != ncoef) {
                warning(paste('FList[[',iterm, ']]$ncoef in modelList[[', 
                              ivar,']] is not an integer.'))
                errwrd <- TRUE
              } else {
                if (is.null(FListij$Ufd)) {
                  warning(paste('List object in FList[[',jterm, 
                                ']] in modelList[[',ivar, 
                                ']] does not have an Ufd field.'))
                  errwrd <- TRUE
                } else {
                  if (!is.fd(FListij$Ufd)) {
                    warning(paste('Ufd field in FList[[',jterm, 
                                  ']] in modelList[[',ivar, 
                                  ']] is not an fd object.'))
                    errwrd <- TRUE
                  }
                }
              }
            }
          }
        }
      }
    }
  }
      
  if (errwrd) {
    stop('One or more terminal errors encountered.')
  }
    
  #  --------------------------------------------------------------------------
  #  check that list objects contain a field named 'factor', and, if
  #  not, or if empty, replace by factor <- 1.
  #  --------------------------------------------------------------------------
  
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    nallXterm <- modelListi$nallXterm
    if (nallXterm > 0) {
      for (iterm in 1:nallXterm) {
        XListij <- modelListi$XList[[iterm]]
        if (is.null(XListij$factor)) {
          XListij$factor <- 1
          modelListi$XList[[iterm]] <- XListij
        }
      }
    }
    nallFterm <- modelListi$nallFterm
    if (nallFterm > 0) {
      for (iterm in 1:nallFterm) {
        FListij <- modelListi$FList[[iterm]]
        if (is.null(FListij$factor)) {
          FListij$factor <- 1
          modelListi$FList[[iterm]] <- FListij
        }
      }
    }
    modelList[[ivar]] <- modelListi
  }
  
  #  --------------------------------------------------------------------------
  #  Check the fields conMap, conVal and orthoMap
  #  Note that modelCheck has no way of knowing the number of basis functions
  #  for this variable.  Conformity of conMap and orthoMap to this will be
  #  checked in function D2LD.
  #  --------------------------------------------------------------------------
  
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    if (is.null(modelListi$conMap)) {
      modelListi$conMap <- NULL
    }
    if (is.null(modelListi$conMap)) {
      ncon <- 0
    } else {
      conMap <- modelListi$conMap
      ncon <- dim(modelListi$conMap)[1]
    }
    if (is.null(modelListi$conVal)) {
      modelListi$conVal <- NULL
    }
    if (!is.null(modelListi$conVal)) {
      modelListi$conVal <- matrix(0,ncon,1)
    }
    if (is.null(modelListi$orthoMap)) {
      modelListi$orthoMap <- NULL
    }
    if (!is.null(modelListi$orthoMap)) {
      ncoef <- dim(conMap)[2]
      qrList <- qr(t(conMap))
      Q <- qrList$Q
      T <- qrList$R
      modelListi$Q1       <- Q[,1:ncon]
      modelListi$orthoMap <- Q[,(ncon+1):ncoef]
      modelListi$Tmat     <- T[1:ncon,1:ncon]
    }
    modelList[[ivar]] <- modelListi
  }
    
  modelListnew <- modelList
    
  return(modelListnew)
    
}
  
  