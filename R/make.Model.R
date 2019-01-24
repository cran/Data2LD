make.Model <- function(XbasisList, variableList, coefList, report=TRUE) {
  # Check the variableList structure, requiring at a minimum:
  #   listarray with a member for each variable, containing:
  #     a list array of length number of homogeneous terms present.  
  #     Each list has:
  #       variable index i and/or tag
  #       order of derivative in the  term
  #       fdPar object for coefficient function
  #   list array of length either number of forcing functions or empty.  
  #     If not empty, lists have:
  #       fd object for forcing function(s)
  #       fdPar object for coefficient function
  # Then construct a model object.
  
  #  Last modified 18 January 2019
  
  #  number of coefficients and number of variables
  
  ncoef <- length(coefList)
  nvar  <- length(variableList)
  
  if (is.null(variableList)) {
    stop('Argument variableList is empty.')
  }
  
  if (!is.list(variableList)) {
    stop('Argument variableList is not a list object.')
  }
  
  #  check that each cell is a struct object
  
  errwrd <- FALSE
  for (ivar in 1:nvar) {
    #  check variableList[[ivar]] is not empty and that it is a list object
    variableListi <- variableList[[ivar]]
    if(is.null(variableListi)) {
      warning(paste('No list object in variableList[[',ivar,']]'))
      errwrd <- TRUE
    }
    if (!is.list(variableListi)) {
      warning(paste('Object in variableList[[',ivar,
                    ']] is not a list object.'))
      errwrd <- TRUE
    }
  }
  
  if (errwrd) stop("One or more variables not defined.")
  
  #  -------------------------------------------------------------------
  #  check that each list is a list object and that it has the necessary
  #  fields or, if not, that these are set to default values
  #  -------------------------------------------------------------------
  
  for (ivar in 1:nvar) {
    variableListi <- variableList[[ivar]]
    #  check variableList[[ivar]] order, name and weight members, 
    #  and assign default if needed
    #  check member "order"
    if (is.null(variableListi$order)) {
      #  assign default value
      variableListi$order <- 1
    }
    #  check member "name"
    if (is.null(variableListi$name)) {
      #  construct a name 
      variableListi$name <- paste("x",ivar,sep="")
    }
    #  check member "weight"
    if (is.null(variableListi$weight)) {
      #  construct unit weight 
      variableListi$weight <- 1
    }
    #  update containing list
    variableList[[ivar]] <- variableListi
  }
  
  #  terminate of an error has been identified
  
  if (errwrd) {
    stop('One or more terminal errors encountered.')
  }

  #  ------------  check the values of each field  --------------------- 
  
  errwrd <- FALSE
  for (ivar in 1:nvar) {
    variableListi <- variableList[[ivar]]
    #  check value of XList
    if (!is.null(variableListi$XList)) { 
      if (!is.list(variableListi$XList)) {
        warning(paste('XList field in variableList[[',ivar, 
                      ']] does not contain a list object.'))
        errwrd <- TRUE
      }
    } 
    #  check value of FList
    if (!is.null(variableListi$FList)) { 
      if (!is.list(variableListi$FList)) {
        warning(paste('FList field in variableList[[',ivar, 
                    ']] does not contain a list object.'))
        errwrd <- TRUE
      }
    }
    #  check value of order
    if (!is.numeric(variableListi$order)) {
      warning(paste('order field in variableList[[',ivar, 
                    ']] is not a number.'))
      errwrd <- TRUE
    }
    #  check value of name
    if (!is.character(variableListi$name)) {
      warning(paste('name field in variableList[[',ivar, 
                    ']] is not a character string.'))
      errwrd <- TRUE
    }
    #  check value of nallXterm
    if (is.null(variableListi$nallXterm)) {
      if (is.null(variableListi$XList)) {
        variableListi$nallXterm <- 0
      } else {
        variableListi$nallXterm <- length(variableListi$XList)
      }
    } else {
      if (!is.numeric(variableListi$nallXterm)) {
        warning(paste('nallXterm field in variableList[[',ivar, 
                      ']] is not a number.'))
        errwrd <- TRUE
      }
    }
    #  check value of nallFterm
    if (is.null(variableListi$nallFterm)) {
      if (is.null(variableListi$FList)) {
        variableListi$nallFterm <- 0
      } else {
        variableListi$nallFterm <- length(variableListi$FList)
      }
    } else {
      if (!is.numeric(variableListi$nallFterm)) {
        warning(paste('nallFterm field in variableList[[',ivar, 
                      ']] is not a number.'))
        errwrd <- TRUE
      }
    }
    variableList[[ivar]] <- variableListi
  }
    
  #  terminate of an error has been identified
  
  if (errwrd) {
    stop('One or more terminal errors encountered.')
  }
  
  #  ----------------  check members of XList  -------------------------
  
  errwrd <- FALSE
  for (ivar in 1:nvar) {
    variableListi <- variableList[[ivar]]
    if (!is.null(variableListi$XList)) {
      nXterm <- length(variableListi$XList)
      for (iterm in 1:nXterm) {
        XListij <- variableListi$XList[[iterm]]
        if (!is.null(XListij) && !is.list(XListij)) {
          warning(paste('XList[[',iterm,']] in variableList[[',ivar, 
                            ']] is not a list object.'))
          errwrd <- TRUE
        } else {
          #  check member ncoef
          if (is.null(XListij$ncoef)) {
            warning(paste('List object in XList[[',iterm, 
                          ']] in variableList[[', 
                          ivar,']] does not have a member ncoef.'))
            errwrd <- TRUE
          } else {
            ncoef <- XListij$ncoef
            if (!is.numeric(ncoef)) {
              warning(paste('XList[[',iterm, ']]$ncoef in variableList[[', 
                            ivar,']] is not a number.'))
              errwrd <- TRUE
            } else {
              if (round(ncoef) != ncoef) {
                warning(paste('XList[[',iterm, 
                              ']]$ncoef in variableList[[', 
                              ivar,']] is not an integer.'))
                errwrd <- TRUE
              }
            }
          }
          #  check member variable   
          if (is.null(XListij$variable)) {
            warning(paste('List object in XList[[',iterm,
                          ']] in variableList[[', 
                          ivar,']] does not have a member variable.'))
            errwrd <- TRUE
          } else {
            varij <- XListij$variable
            if (!is.numeric(varij)) {
              warning(paste('XList[[',iterm, 
                            ']]$variable in variableList[[', 
                            ivar,']] is not a number.'))
              errwrd <- TRUE
            } else {
              if (round(varij) != varij) {
                warning(paste('XList[[',iterm, 
                              ']]$ncoef in variableList[[', 
                              ivar,']] is not an integer.'))
                errwrd <- TRUE
              } else {
                if (varij > nvar) {
                  warning(paste('List object in XList[[',iterm, 
                                ']] in variableList[[',ivar, 
                                ']] has variable field > NVAR.'))
                  errwrd <- TRUE
                }
              }
            }
          }
          #  check member derivative
          if (is.null(XListij$derivative)) {
            warning(paste('List object in XList[[',iterm, 
                          ']] in variableList[[',ivar, 
                          ']] does not have a member derivative.'))
            errwrd <- TRUE
          } else {
            derivative <- XListij$derivative
            if (!is.numeric(derivative)) {
              warning(paste('XList[[',iterm, 
                            ']]$derivative in variableList[[', 
                            ivar,']] is not a number.'))
              errwrd <- TRUE
            } else {
              if (round(derivative) != derivative) {
                warning(paste('XList[[',iterm, 
                              ']]$derivative in variableList[[', 
                              ivar,']] is not an integer.'))
                errwrd <- TRUE
              } else {
                if (XListij$derivative >= variableListi$order) {
                  warning(paste('List object in XList[[',iterm, 
                                ']] in variableList[[',ivar, 
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
   
  #  ---------------  check members of FList  --------------------------
  
  errwrd <- FALSE
  for (ivar in 1:nvar) {
    variableListi <- variableList[[ivar]]
    if (!is.null(variableListi$FList)) {
      FList <- variableListi$FList
      nFterm <- length(FList)
      for (jterm in 1:nFterm) {
        FListij <- FList[[jterm]]
        if (!is.null(FListij) && !is.list(FListij)) {
          warning(paste('FList[[',jterm,']] in variableList[[',ivar, 
                        ']] is not a list object.'))
          errwrd <- TRUE
        } else {
          if (is.null(FListij$ncoef)) {
            warning(paste('List object in FList[[',jterm,
                          ']] in variableList[[',ivar,
                          ']] does not have a ncoef field.'))
            errwrd <- TRUE
          } else {
            ncoef <- FListij$ncoef
            if (!is.numeric(ncoef)) {
              warning(paste('FList[[',iterm, ']]$ncoef in variableList[[', 
                            ivar,']] is not a number.'))
              errwrd <- TRUE
            } else {
              if (round(ncoef) != ncoef) {
                warning(paste('FList[[',iterm, 
                              ']]$ncoef in variableList[[', 
                              ivar,']] is not an integer.'))
                errwrd <- TRUE
              } else {
                if (is.null(FListij$Ufd)) {
                  warning(paste('List object in FList[[',jterm, 
                                ']] in variableList[[',ivar, 
                                ']] does not have an Ufd field.'))
                  errwrd <- TRUE
                } else {
                  if (!is.fd(FListij$Ufd)) {
                    warning(paste('Ufd field in FList[[',jterm, 
                                  ']] in variableList[[',ivar, 
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
    
  #  -------------------------------------------------------------------
  #  check that list objects contain a field named 'factor', and, if
  #  not, or if empty, replace by factor <- 1.
  #  -------------------------------------------------------------------
  
  for (ivar in 1:nvar) {
    variableListi <- variableList[[ivar]]
    nallXterm <- variableListi$nallXterm
    if (nallXterm > 0) {
      for (iterm in 1:nallXterm) {
        XListij <- variableListi$XList[[iterm]]
        if (is.null(XListij$factor)) {
          XListij$factor <- 1
          variableListi$XList[[iterm]] <- XListij
        }
      }
    }
    nallFterm <- variableListi$nallFterm
    if (nallFterm > 0) {
      for (iterm in 1:nallFterm) {
        FListij <- variableListi$FList[[iterm]]
        if (is.null(FListij$factor)) {
          FListij$factor <- 1
          variableListi$FList[[iterm]] <- FListij
        }
      }
    }
    variableList[[ivar]] <- variableListi
  }
  
  #  -------------------------------------------------------------------
  #  If report is TRUE, print out a summary of the model
  #  -------------------------------------------------------------------

  if (report) {
    for (ivar in 1:nvar) {
      variableListi <- variableList[[ivar]]
      cat("Model Summary\n")
      cat("---------------------------\n")
      cat("Variable",ivar,",",variableListi$name,
          "with derivative order",variableListi$order,
          "and weight",variableListi$weight,"\n")
      XList = variableListi$XList
      if (!is.null(XList)) {
        cat("Homogeneous term(s) in the equation:\n")
        nXList = length(XList)
        for (j in 1:nXList) {
          XListj = XList[[j]]
          cat("    Term",j,XListj$name,"\n")
          cat("    has coefficient number",XListj$ncoef,
                  "and constant factor",XListj$factor,".\n")
          cat("    They multiply variable",XListj$variable,
                  "with derivative order",XListj$derivative,".\n")
        }
      }
      FList = variableListi$FList
      if (!is.null(FList)) {
        cat("Forcing term(s) in the equation:\n")
        nFList = length(FList)
        for (j in 1:nFList) {
          FListj = FList[[j]]
          cat("    Term",j,FListj$name,"\n")
          cat("    has coefficient number",FListj$ncoef,
                   "and constant factor",FListj$factor,".\n")
        }
      }
    }
    cat("---------------------------\n")
  }

  #  -------------------------------------------------------------------
  #  Set up the four-way tensors and save them as fields 
  #  -------------------------------------------------------------------
  
  BtensorList  <-   Btensorfn(XbasisList, variableList, coefList)
  BAtensorList <-  BAtensorfn(XbasisList, variableList, coefList)
  AtensorList  <-   Atensorfn(            variableList, coefList)
  
  for (ivar in 1:nvar) {
    variableListi <- variableList[[ivar]]
    variableListi$Btens  <-  BtensorList[[ivar]]
    variableListi$BAtens <- BAtensorList[[ivar]]
    variableListi$Atens  <-  AtensorList[[ivar]]
    variableList[[ivar]] <- variableListi
  }
  
  modelList <- variableList
  
  return(modelList)
    
}
  
  