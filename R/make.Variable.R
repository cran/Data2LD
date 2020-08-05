make.Variable <- function(name="", order=1, XList=NULL, FList=NULL, weight=1) {
  #  make_Variable assembles four arguments into a struct object that is
  #  used by function make_Model  to set up a linear dynamicsystem object.
  #  Arguments are:
  #  NAME  ... A string to be used as the name of the variable
  #  ORDER ... The order of the derivative on the left side
  #  XCELL ... A cell array containing specifications for the 
  #            homogeneous terms.  See make_Xterm for more details.
  #  FCELL ... A cell array containing specifications for tthe 
  #            forcing terms.  See make_Fterm for more details.
  
  #  Last modified 16 April 2020
  
  if (floor(order) != order || order < 1) {
    stop("Argument ORDER is not a positive integer.")
  }
  if (!is.null(XList) && !is.null(FList)) 
    termList <- list(name=name, order=order, XList=XList, FList=FList, weight=weight)
  if (!is.null(XList) &&  is.null(FList)) 
    termList <- list(name=name, order=order, XList=XList, weight=weight)
  if ( is.null(XList) && !is.null(FList)) 
    termList <- list(name=name, order=order, FList=FList, weight=weight)
  if ( is.null(XList) &&  is.null(FList)) 
    termList <- list(name=name, order=order, weight=weight)
  return(termList)
}