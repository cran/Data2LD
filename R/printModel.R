printModel <- function(modelList, titlestr=NULL) {
  if (!is.null(titlestr)) {
    print(" ")
    print(titlestr)
    print(" ")
  }
  nvar <- length(modelList)
  for (ivar in 1:nvar) {
    #  print variable information
    modelListi <- modelList[[ivar]]
    cat("----------------------------------------------------------------------\n")
    cat(paste("Variable",ivar,  modelListi$name,", order =",
              modelListi$order, ", weight =", modelListi$weight,"\n"))
    #  print information for (each homogeneous term
    cat("    Homogenous term(s):\n") 
    if (modelListi$nallXterm > 0) {
      for (iw in 1:modelListi$nallXterm) {
        modelListiw <- modelListi$XList[[iw]]
        #  print number, variable, derivative and factor
        cat(paste("    Homogenous term",iw,
                  ", variable",    modelListiw$variable,
                  ", derivative =",modelListiw$derivative,
                  ", factor =",    modelListiw$factor,"\n")) 
        #  print positions in parameter vector
        nindex <- length(modelListiw$index)
        cat("    Parameter vector position(s):\n")
        if (nindex == 1) {
          cat(paste("    ",modelListiw$index))
        } else {
          for (k in 1:nindex) {
            cat(paste("    ",modelListiw$index[k],", "))
          }
        }
        cat("\n")
        #  print parvec values
        cat("    Parameter vector value(s):\n")
        if (length(unique(modelListiw$parvec)) == 1) {
          cat(paste("    All values =", modelListiw$parvec[1]))
        } else {
          for (k in 1:nindex) {
            cat(paste("    ",round(modelListiw$parvec[k],3),", "))
          }
        }
        cat("\n")
        funobj = modelListiw$fun
        if (!(is.basis(funobj) || is.fd(funobj) || is.fdPar(funobj))) {
          cat("    Coefficient type is general\n")
        } else {
          cat("    Coefficient type is functional data object\n")
        }
        estimate <- modelListiw$estimate        
        if (sum(estimate) == nindex) {
          cat("    All parameters are estimated.\n")
        }
        if (sum(estimate) == 0) {
          cat("    All parameters are fixed.\n")
        }
        if (nindex > 1 && sum(estimate) > 0 && sum(estimate) < nindex) {
          cat("    Parameter estimation status:\n")
          for (k in 1:nindex) {
            cat(paste(modelListiw$estimate[k],", "))
          }
          cat("\n")
          cat("\n")
        }     
      }
    }
    
    #  print information for each forcing term
    if (modelListi$nallFterm > 0) {
      cat("    ----------------------------------------\n")
      cat(paste("    Forcing terms:\n"))
      for (j in 1:modelListi$nallFterm) {
        modelListij <- modelListi$FList[[j]]
        #  print number, variable, derivative and factor
        cat(paste("    Forcing term",j,", factor =",modelListij$factor,"\n"))
        #  print positions in parameter vector
        nindex <- length(modelListij$index)
        cat("    Parameter vector position(s):\n")
        if (nindex == 1) {
          cat(paste("    ",modelListij$index))
        } else {
          for (k in 1:nindex) {
            cat(paste("  ",modelListij$index[k]))
          }
        }
        cat("\n")
        #  print parvec values
        cat("    Parameter vector value(s):\n")
        if (length(unique(modelListij$parvec)) == 1) {
          cat(paste("    All values = ", modelListij$parvec[1]))
        } else {
          cat(paste("    ", modelListij$parvec[1]))
          for (k in 2:nindex) {
            cat(paste("  ",modelListij$parvec[k]))
          }
          cat("\n")
        }
        cat("\n")
        funobj = modelListij$fun
        if (!(is.basis(funobj) || is.fd(funobj) || is.fdPar(funobj))) {
          cat("    Coefficient type is general\n")
        } else {
          cat("    Coefficient type is functional data object\n")
        }
        estimate <- modelListij$estimate        
        if (sum(estimate) == nindex) {
          cat("    All parameters are estimated.\n")
        }
        if (sum(estimate) == 0) {
          cat("    All parameters are fixed.\n")
        }
        if (nindex > 1 && sum(estimate) > 0 && sum(estimate) < nindex) {
          cat("    Parameter estimation status:\n")
          cat(paste("    ", modelListij$estimate[1]))
          for (k in 2:nindex) {
            cat(paste("  ",modelListij$estimate[k]))
          }
          cat("\n")
          cat("\n")
        }
      }
    }
  }
  cat("----------------------------------------------------------------------\n")
  cat("\n")
}