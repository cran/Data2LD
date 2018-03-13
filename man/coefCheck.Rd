\name{coefCheck}
\alias{coefCheck}
\title{
  Check a list object containing coefficient specifications.
}
\description{
  A system of linear differential equations is defined by a list of length equal to the 
  number of variables in the system.  Each member of this list defines a single linear 
  differential equation.  Within each equation there are typically one or more terms 
  involving a coefficient function multiplying either a derivative of one of the 
  variables, or multiplying a forcing function.  The coefficient functions in term 
  are defined by one or more parameters.  For some coefficients, its parameters are 
  regarded as fixed, but for others, they require estimation from the data.  This 
  function checks that each of the coefficient specifications has been properly set 
  up using values that are appropriate.
  In addition, it also assembles a numeric vector containing the values of of the 
  estimated parameters, and returns this as well as the number of values in a named list.  
  It is easy to make mistakes in setting up lists of coefficient objects,
  and it is strongly suggested to invoke this function prior to invoking any function 
  using this list.
}
\usage{
  coefCheck(coefList)
}
\arguments{
  \item{coefList}{A list object containing the specifications of one or more coefficient 
  functions.}                                      
}
\details{
  A coefficient specification can be made manually, or can set up in a by a single 
  invocation of function \code{make.coef}. This function checks that all coefficient 
  specifications encountered when processing a model List, whether or not the 
  coefficients are to be estimated, are consistent with what is required to 
  specify a coefficient.  If certain values are not supplied, default values are 
  assigned.  For example, if the member with name \code{estimate} is not supplied, 
  it is set to TRUE.  If one or more parameter values are not specified, they are 
  set to zeros.
}
\value{
  A named list object. Member \code{coefList} is an image of that argument where all 
  default values are assigned.  This should be defined as the version to be used in 
  subsequent calculations.  Member \code{theta} is a numeric vector containing the 
  values of all parameters requiring estimation.  Member \code{ntheta} is the length 
  of this parameter vector.
}
\references{
    J. O. Ramsay and G. Hooker (2017) \emph{Dynamic Data Analysis}. Springer.
}
\seealso{
  \code{\link{make.coef}},
  \code{\link{make.Xterm}},
  \code{\link{make.Fterm}},
  \code{\link{make.variable}},
  \code{\link{modelCheck}}
}
\examples{
#  For examples of the use of this function, 
#  see these examples in the description of the function \code{make.variable}.
#
#  Example 1:  The refinery data
#
#
#  Example 2:  The X coordinate of the fda script data
#
#
#  Example 3:  The X coordinate of the fda script data
#
#
#  Example 4  Average temperature for Montreal
#
}
