\name{make.Xterm}
\alias{make.Xterm}
\title{
  Check homogeneous term specifications for a linear differential equation.
}
\description{
  A system of linear differential equations is defined by a list of length equal to the 
  number of variables in the system.  Each member of this list defines a single linear 
  differential equation.  Within each equation there are typically one or more terms 
  involving a coefficient function multiplying either a derivative of one of the 
  variables in the system.  These terms often involve a fixed constant multiplier, 
  which is frequently either 1 or -1.
  This function sets up a list object that specifies the structure of a single term.
}
\usage{
make.Xterm(variable, ncoef, derivative=0, factor=1)
}
\arguments{
  \item{variable}{An integer specifying the variable within which this is a term.}
  \item{ncoef}{An integer specifying the position of the coefficient specification 
  in a coefficient list.}
  \item{derivative}{The order of the derivative of the variable.}
  \item{factor}{A real number that is treated as fixed. For example, it is frequently 
  the case that a variable will appear in two or more places in a system of equations, 
  and sometimes multipled by -1.
  }                                      
}
\details{
  This function checks that all supplied terms conform to what is required.
}  
\value{
  A named list object defining a homogeous term in a linear differential equation.
}
\references{
    J. O. Ramsay and G. Hooker (2017) \emph{Dynamic Data Analysis}. Springer.
}
\seealso{
  \code{\link{make.coef}},
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
