\name{modelList2Vec}
\alias{modelList2Vec}
\title{
  Extract the master parameter vector from a model list object.
}
\description{
  A system of linear differential equations is defined by a list of length equal to the 
  number of variables in the system.  Each member of this list defines a single linear 
  differential equation.  Within each equation there are typically one or more terms 
  involving a coefficient function multiplying either a derivative of one of the 
  variables, or multiplying a forcing function.  
  The coefficient functions in term are defined by one or more parameters.  For some
  coefficients, its parameters are regarded as fixed, but for others, they require 
  estimation from the data.  This function extracts all the parameters requiring 
  estimation and returns them as a one-dimensional vector.
}
\usage{
modelList2Vec(modelList, nparam)
}
\arguments{
  \item{modelList}{A list object containing the specification of a Data2LD model. 
  Each member of this list contains a list object that defines a single linear 
  differential equation.}
  \item{nparam}{The total number of parameters in the master parameter vector,
  which contains all parameters, whether estimated or not.}                                      
}
\details{
  The parameters in the master parameter vector are organized as follows.  
  First come the parameters defining the homogeneous coefficient functions, in the
  order in which the appear in the equation system.  Next come the parameters
  defining the forcing coefficient functions, also in the order in which they 
  appear.  
  
  Normally the user will have a limited need to deal with the master parameter vector.  
  Both the initial and final parameter values are displayed by the function
  \code{printModel} in a much more readable form.
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
  \code{\link{modelVec2List}}
}
