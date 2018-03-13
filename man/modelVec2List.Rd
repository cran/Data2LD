\name{modelVec2List}
\alias{modelVec2List}
\title{
  Put values in a parameter vector into the coefficient specifications within a model 
  List object.
}
\description{
  A system of linear differential equations is defined by a list of length equal to 
  the number of variables in the system.  Each member of this list defines a single 
  linear differential equation.  Within each equation there are typically one or more 
  terms involving a coefficient function multiplying either a derivative of one of the 
  variables, or multiplying a forcing function.  The coefficient functions in term are 
  defined by one or more parameters.  For some coefficients, its parameters are regarded
  as fixed, but for others, they require estimation from the data.  This  function 
  replaces all the parameters in the coefficient functions in a coefficient list object 
  requiring estimation by values within a one-dimensional vector.
}
\usage{
modelVec2List(thetavec, coefList)
}
\arguments{
  \item{thetavec}{A numeric vector containing values for all the parameters in a 
  linear differential system that requires estimation.}
  \item{coefList}{A list object containing the specifications of one or more 
  coefficient functions.}                                      
}
\details{
  A coefficient specification can be made manually, or can set up in a by a single 
  invocation of function \code{make.coef}.  
}
\value{
  A numerical vector object containing values of the parameters to be estimated.  
  The order of these parameters is the order in which they occur as each member of a 
  model list object in turn is processed.
}
\references{
    J. O. Ramsay and G. Hooker (2017) \emph{Dynamic Data Analysis}. Springer.
}
\seealso{
  \code{\link{modelList2Vec}}
}
