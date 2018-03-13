\name{make.Fterm}
\alias{make.Fterm}
\title{
  Check each specification in list of forcing term specifications for a 
  linear differential equation.
}
\description{
  A system of linear differential equations is defined by a list of length equal to the 
  number of variables in the system.  Each member of this list defines a single linear 
  differential equation.  Within each equation there may be one or more terms involving 
  a coefficient function multiplying either a known fixed function, called a 
  forcing function.  Forcing terms often involve a fixed constant multiplier, 
  which is frequently either 1 or -1.
  This function sets up a list object that specifies the structure of a single term.
}
\usage{
make.Fterm(ncoef, Ufd, factor=1)
}
\arguments{
  \item{ncoef}{An integer specifying the position of the coefficient specification 
  in a coefficient list.}
  \item{Ufd}{A known forcing function.  This may be specified as a functional data 
  or functional parameter object, but may also be specified as a user-defined function.  
  See function \code{make.coef} for more details.}
  \item{factor}{A real number that is treated as fixed. For example, it is frequently 
  the case that a variable will appear in two or more places in a system of equations, 
  and sometimes multipled by -1.}                                      
}
\details{
  This function checks that all supplied terms conform to what is required.
}  
\value{
  A named list object defining a forcing term in a linear differential equation.
}
\references{
    J. O. Ramsay and G. Hooker (2017) \emph{Dynamic Data Analysis}. Springer.
}
\seealso{
  \code{\link{make.coef}},
  \code{\link{make.Xterm}},
  \code{\link{make.variable}},
  \code{\link{modelCheck}}
}
\examples{
#  For examples of the use of this function, 
#  see these examples in the description of the function make.variable.
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
