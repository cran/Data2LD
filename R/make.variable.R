make.variable <- function(name="", order=1, XList=NULL, FList=NULL) {
  if (floor(order) != order || order < 1) {
    stop("Argument ORDER is not a positive integer.")
  }
  if (!is.null(XList) && !is.null(FList)) 
    termList <- list(name=name, order=order, XList=XList, FList=FList)
  if (!is.null(XList) &&  is.null(FList)) 
    termList <- list(name=name, order=order, XList=XList)
  if ( is.null(XList) && !is.null(FList)) 
    termList <- list(name=name, order=order, FList=FList)
  if ( is.null(XList) &&  is.null(FList)) 
    termList <- list(name=name, order=order)
  return(termList)
}