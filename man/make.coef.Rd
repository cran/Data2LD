\name{make.coef}
\alias{make.coef}
\title{
  Check each specification in list of linear differential equation specifications.
}
\description{
  A system of linear differential equations is defined by a list of length equal to the 
  number of variables in thesystem.  Each member of this list defines a single linear 
  differential equation.  Within each equation there are typically one or more terms 
  involving a coefficient function multiplying either a derivative of one of the 
  variables, or multiplying a forcing function.  This function sets up a list object 
  that specifies the structure of a single coefficient function.
}
\usage{
make.coef(funobj, parvec, estimate=TRUE, coeftype="beta")
}
\arguments{
  \item{funobj}{This is the definition coefficient function.  This may be specified 
  in two ways:  most commonly as a functional data or a functional parameter object, 
  or alternative by a list containing as members a function object for evaluating the
  function and a function object for evaluating its gradient with respect to the 
  parameters defining it.}
  \item{parvec}{A numeric constant or vector of values of parameters defining the 
  function defined by argument \code{funobj.}}
  \item{estimate}{A logical value indicating whether the parameters are to be 
  estimated (TRUE) or regarded as fixed (FALSE).}
  \item{coeftype}{Whether the coefficient multiplies a variable in the system ("beta")
  or multiplies a forcing term ("alpha").  This information is not used in this 
  release of the package, but may be needed in subsequent releases.}
}
\details{
  This function checks that all supplied conform to what is required.  In the case 
  where a user-defined coefficient function is provided, argument \code{funobj} 
  is a list containing four named members: \code{fd} is a function object 
  returning the value of the coefficient function, \code{Dfd} returns its gradient with 
  respect to the parameters, \code{theta} is a numeric vector containing values for the 
  parameters defining the function, and \code{more} is a list containing any other 
  objects required by \code{fd} and \code{Dfd}.  
}
\value{A named list object defining a coefficient object.
}
\references{
    J. O. Ramsay and G. Hooker (2017) \emph{Dynamic Data Analysis}. Springer.
}
\seealso{
  \code{\link{make.Xterm}},
  \code{\link{make.Fterm}},
  \code{\link{make.variable}},
  \code{\link{coefCheck}}
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
