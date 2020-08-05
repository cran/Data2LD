checkModel <- function(XbasisList, modelList, summarywrd=TRUE) {
  # Check the modelList Listure, requiring at a minimum:
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
  # Then conList a model object.
  
  #  Last modified 17 April 2020
  
  #  number of coefficients and number of variables
  
  if (is.null(modelList)) {
    stop("Argument modelList is empty.")
  }
  
  if (!is.list(modelList)) {
    stop("Argument modelList is not a list object.")
  }
  
  #  check that within each modelList is a List object
  
  nvar  <- length(modelList)
  
  errwrd <- FALSE
  for (ivar in 1:nvar) {
    #  check modelList[[ivar]] is not empty and that it is a list object
    modelListi <- modelList[[ivar]]
    if(is.null(modelListi)) {
      warning(paste("No list object in modelList[[",ivar,"]]"))
      errwrd <- TRUE
    }
    if (!is.list(modelListi)) {
      warning(paste("Object in modelList[[",ivar,
                    "]] is not a list object."))
      errwrd <- TRUE
    }
  }
  
  if (errwrd) stop("One or more variables not defined.")
  
  #  -------------------------------------------------------------------
  #  check that each list has the necessary
  #  fields or, if not, that these are set to default values
  #  -------------------------------------------------------------------
  
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    #  check modelList[[ivar]] order, name and weight members, 
    #  and assign default if needed
    #  check member "order"
    if (is.null(modelListi$order)) {
      #  assign default value
      modelListi$order <- 1
    }
    #  check member "name"
    if (is.null(modelListi$name)) {
      #  conList a name 
      modelListi$name <- paste("x",ivar,sep="")
    }
    #  check member "weight"
    if (is.null(modelListi$weight)) {
      #  conList unit weight 
      modelListi$weight <- 1
    }
    #  update containing list
    modelList[[ivar]] <- modelListi
  }
  
  #  terminate of an error has been identified
  
  if (errwrd) {
    stop("One or more terminal errors encountered.")
  }
  
  #  ------------  check the values of each field  --------------------- 
  
  errwrd <- FALSE
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    order <- modelListi$order
    #  check order field
    if (!is.numeric(order)) {
      warning(paste("order field in modelList[[",ivar,
                    "]] does not contain a numeric object."))
      errwrd <- TRUE
    } else {
      if (order != floor(order)) {
        warning(paste("Field order in modelList[[",ivar,
                      "]] is not an integer."))
        errwrd <- TRUE
      }
    }
    #  check value of XList
    if (!is.null(modelListi$XList)) { 
      if (!is.list(modelListi$XList)) {
        warning(paste("XList field in modelList[[",ivar, 
                      "]] does not contain a list object."))
        errwrd <- TRUE
      }
    } 
    #  check value of FList
    if (!is.null(modelListi$FList)) { 
      if (!is.list(modelListi$FList)) {
        warning(paste("FList field in modelList[[",ivar, 
                      "]] does not contain a list object."))
        errwrd <- TRUE
      }
    }
    #  check value of name
    if (!is.character(modelListi$name)) {
      warning(paste("name field in modelList[[",ivar, 
                    "]] is not a character string."))
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
        warning(paste("nallXterm field in modelList[[",ivar, 
                      "]] is not a number."))
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
        warning(paste("nallFterm field in modelList[[",ivar, 
                      "]] is not a number."))
        errwrd <- TRUE
      }
    }
    modelList[[ivar]] <- modelListi
  }
  
  #  terminate of an error has been identified
  
  if (errwrd) {
    stop("One or more terminal errors encountered.")
  }
  
  #  ----------------  check members of XList  -------------------------
  
  errwrd <- FALSE
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    if (!is.null(modelListi$XList)) {
      nXterm <- length(modelListi$XList)
      for (iterm in 1:nXterm) {
        XListij <- modelListi$XList[[iterm]]
        if (!is.null(XListij) && !is.list(XListij)) {
          warning(paste("XList[[",iterm,"]] in modelList[[",ivar, 
                        "]] is not a list object."))
          errwrd <- TRUE
        } else {
          #  check member variable   
          if (is.null(XListij$variable)) {
            warning(paste("List object in XList[[",iterm,
                          "]] in modelList[[", 
                          ivar,"]] does not have a member variable."))
            errwrd <- TRUE
          } else {
            varij <- XListij$variable
            if (!is.numeric(varij)) {
              warning(paste("XList[[",iterm, 
                            "]]$variable in modelList[[", 
                            ivar,"]] is not a number."))
              errwrd <- TRUE
            } else {
              if (round(varij) != varij) {
                warning(paste("XList[[",iterm, 
                              "]]$ncoef in modelList[[", 
                              ivar,"]] is not an integer."))
                errwrd <- TRUE
              } else {
                if (varij > nvar) {
                  warning(paste("List object in XList[[",iterm, 
                                "]] in modelList[[",ivar, 
                                "]] has variable field > NVAR."))
                  errwrd <- TRUE
                }
              }
            }
          }
          #  check member derivative
          if (is.null(XListij$derivative)) {
            warning(paste("List object in XList[[",iterm, 
                          "]] in modelList[[",ivar, 
                          "]] does not have a member derivative."))
            errwrd <- TRUE
          } else {
            derivative <- XListij$derivative
            if (!is.numeric(derivative)) {
              warning(paste("XList[[",iterm, 
                            "]]$derivative in modelList[[", 
                            ivar,"]] is not a number."))
              errwrd <- TRUE
            } else {
              if (round(derivative) != derivative) {
                warning(paste("XList[[",iterm, 
                              "]]$derivative in modelList[[", 
                              ivar,"]] is not an integer."))
                errwrd <- TRUE
              } else {
                if (XListij$derivative >= modelListi$order) {
                  warning(paste("List object in XList[[",iterm, 
                                "]] in modelList[[",ivar, 
                                "]] has derivative member >= ORDER."))
                  errwrd <- TRUE
                }
              }
            }
          }
          # check member funobj
          funobj <- XListij$funobj
          if (is.null(funobj)) {
            warning(paste("List object in XList[[",iterm,
                          "]] in modelList[[",ivar,
                          "]] does not have a member funobj."))
            errwrd <- TRUE
          } else {
            if (!(is.basis(funobj) || is.fd(funobj) ||
                  is.fdPar(funobj) || is.list(funobj))) {
              warning(paste("XList[[",iterm,
                            "]]$funobj in modelList[[",ivar,
                            "]] is not of class ",
                            "basis, fd, fdPar or function."))
              errwrd <- TRUE
            }
          }
          # check member parvec
          parvec <- XListij$parvec
          if (is.null(parvec)) {
            warning(paste("List object in XList[[",iterm,
                          "]] in modelList[[",ivar,
                          "]] does not have a member parvec."))
            errwrd <- TRUE
          } else {
            if (!is.numeric(parvec)) {
              warning(paste("XList[[",iterm,
                            "]]$parvec in modelList[[",ivar,
                            "]] is not numeric."))
              errwrd <- TRUE
            }
          }
          # check member index
          index <- XListij$index
          if (!is.null(index)) {
            if (!is.numeric(index)) {
              warning(paste("XList[[",iterm,
                            "]]$index in modelList[[",ivar,
                            "]] is not numeric."))
              errwrd <- TRUE
            } else {
              if (any(index != floor(index))) {
                warning(paste("XList[[",iterm,
                              "]]$index in modelList[[",ivar,
                              "]] is not all integers."))
                errwrd <- TRUE
              }
            }
          }
          # check member estimate
          estimate <- XListij$estimate
          if (is.null(estimate)) {
            warning(paste("List object in XList[[",iterm,
                          "]] in modelList[[",ivar,
                          "]] does not have a member estimate"))
            errwrd <- TRUE
          } else {
            if (!(is.logical(estimate))) {
              warning(paste("XList[[",iterm,
                            "]]$estimate in modelList[[",ivar,
                            "]] is not logical."))
              errwrd <- TRUE
            }
          }
        }
      }
    }
  }
  
  #  ---------------  check members of FList  --------------------------
  
  errwrd <- FALSE
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    if (!is.null(modelListi$FList)) {
      FList <- modelListi$FList
      nFterm <- length(FList)
      for (jterm in 1:nFterm) {
        FListij <- FList[[jterm]]
        if (!is.null(FListij) && !is.list(FListij)) {
          warning(paste("FList[[",jterm,"]] in modelList[[",ivar, 
                        "]] is not a list object."))
          errwrd <- TRUE
        } else {
          Ufd <- FListij$Ufd
          if (is.null(Ufd)) {
            warning(paste("List object in FList[[",jterm, 
                          "]] in modelList[[",ivar, 
                          "]] does not have an Ufd field."))
            errwrd <- TRUE
          } else {
            if (!is.fd(Ufd)) {
              warning(paste("Ufd field in FList[[",jterm, 
                            "]] in modelList[[",ivar, 
                            "]] is not an fd object."))
              errwrd <- TRUE
            }
          }
        }
      }
      # check member funobj
      funobj <- FListij$funobj
      if (is.null(funobj)) {
        warning(paste("List object in FList[[",iterm,
                      "]] in modelList[[",ivar,
                      "]] does not have a member funobj."))
        errwrd <- TRUE
      } else {
        if (!(is.basis(funobj) || is.fd(funobj) ||
              is.fdPar(funobj) || is.list(funobj))) {
          warning(paste("FList[[",iterm,
                        "]]$funobj in modelList[[",ivar,
                        "]] is not of class ",
                        "basis, fd, fdPar or function."))
          errwrd <- TRUE
        }
      }
      # check member parvec
      parvec <- FListij$parvec
      if (is.null(parvec)) {
        warning(paste("List object in FList[[",iterm,
                      "]] in modelList[[",ivar,
                      "]] does not have a member parvec."))
        errwrd <- TRUE
      } else {
        if (!is.numeric(parvec)) {
          warning(paste("FList[[",iterm,
                        "]]$parvec in modelList[[",ivar,
                        "]] is not numeric."))
          errwrd <- TRUE
        }
      }
      # check member index
      index <- FListij$index
      if (!is.null(index)) {
        if (!is.numeric(index)) {
          warning(paste("FList[[",iterm,
                        "]]$index in modelList[[",ivar,
                        "]] is not numeric."))
          errwrd <- TRUE
        } else {
          if (any(index != floor(index))) {
            warning(paste("FList[[",iterm,
                          "]]$index in modelList[[",ivar,
                          "]] is not all integers."))
            errwrd <- TRUE
          }
        }
      }
      # check member estimate
      estimate <- FListij$estimate
      if (is.null(estimate)) {
        warning(paste("List object in FList[[",iterm,
                      "]] in modelList[[",ivar,
                      "]] does not have a member estimate"))
        errwrd <- TRUE
      } else {
        if (!(is.logical(estimate))) {
          warning(paste("FList[[",iterm,
                        "]]$funobj in modelList[[",ivar,
                        "]] is not logical."))
          errwrd <- TRUE
        }
      }
    }
  }
  
  if (errwrd) {
    stop("One or more terminal errors encountered.")
  }
  
  #  -------------------------------------------------------------------
  #  check that list objects contain a field named "factor", and, if
  #  not, or if empty, replace by factor <- 1.
  #  -------------------------------------------------------------------
  
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
  
  #  compute total number of parameters
  
  #  homogeneous term portion of parameter vector
  
  m2 <- 0
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    if (modelListi$nallXterm > 0) {
      for (iw in 1:modelListi$nallXterm) {
        modelListiw <- modelListi$XList[[iw]]
        parveciw    <- modelListiw$parvec
        m1 <- m2 + 1
        m2 <- m2 + length(parveciw)
        modelListiw$index <- m1:m2
        modelListi$XList[[iw]] <- modelListiw
      }
    }
    modelList[[ivar]] <- modelListi
  }
  
  #  forcing term portion of parameter vector
  
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    if (modelListi$nallFterm > 0) {
      for (jforce in 1:modelListi$nallFterm) {
        modelListj <- modelListi$FList[[jforce]]
        parvecj    <- modelListj$parvec
        m1 <- m2 + 1
        m2 <- m2 + length(parvecj)
        modelListj$index <- m1:m2
        modelListi$FList[[jforce]] <- modelListj
      }
    }
    modelList[[ivar]] <- modelListi
  }
  
  nparam <- m2
  
  #  -------------------------------------------------------------------
  #  Set up the four-way tensors and save them as fields 
  #  -------------------------------------------------------------------
  
  #  first check XbasisList
  
  if (!is.list(XbasisList)) {
    stop("Argument XbasisList is not a list object")
  }
  
  if (length(XbasisList) != nvar) {
    stop("Length of argument XbasisList is not number of variables.")
  }
  
  errwrd = FALSE
  for (ivar in 1:nvar) {
    if (!is.basis(XbasisList[[ivar]])) {
      if (is.fd(XbasisList[[ivar]])) {
        XbasisList[[ivar]] = XbasisList[[ivar]]$basis
      } else {
        if (is.fdPar(XbasisList[[ivar]])) { 
          XbasisList[[ivar]] = XbasisList[[ivar]]$fd$basis
        } else {
          errwrd = TRUE
        }
      }
    }
  }
  
  if (errwrd) {
    stop(paste("One or more objects in argument XbasisList ", 
               "are not basis, fd or fdPar objects."))
  }
  
  # proceed with computing tensors
  
  BtensorList  <-   Btensorfn(XbasisList, modelList)
  BAtensorList <-  BAtensorfn(XbasisList, modelList)
  AtensorList  <-   Atensorfn(            modelList)
  
  for (ivar in 1:nvar) {
    modelListi <- modelList[[ivar]]
    modelListi$Btens  <-  BtensorList[[ivar]]
    modelListi$BAtens <- BAtensorList[[ivar]]
    modelListi$Atens  <-  AtensorList[[ivar]]
    modelList[[ivar]] <- modelListi
  }
  
  modelList1 = modelList[[1]]$Btens
  modelList <- modelList
  
  if (summarywrd) printModel(modelList)
  
  return(list(modelList=modelList, nparam=nparam))
  
}

