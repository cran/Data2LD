\name{modelCheck}
\alias{modelCheck}
\title{
  Check linear differential equation specifications.
}
\description{
  A system of linear differential equations is defined by a list of length equal 
  to the number of variables in the system.  Each member of this list defines a single 
  linear differential equation.  Within each equation there are typically one or more 
  terms involving a coefficient function multiplying either a derivative of one of the 
  variables, or multiplying a forcing function.  This function checks that all linear 
  differential equation specifications conform to what is required.  
  Default values are supplied where needed.  New named members are also supplied.  
  It is easy to make mistakes in setting up lists of variable specifications,
  and it is ESSENTIAL to invoke this function prior to invoking any function using 
  this list.
}
\usage{
modelCheck(modelList, coefList)
}
\arguments{
  \item{modelList}{A list object containing the specification of a Data2LD model. 
  Each member of this list contains a list object that defines a single linear 
  differential equation.}
  \item{coefList}{A list object containing the specifications of one or more 
  coefficient functions.}                                      
}
\details{
  A model list object containing variable specifications can be made manually, 
  or it can set up in a by a single invocation of function \code{make.variable}. 
  This function checks that all variable specifications encountered when processing 
  a model List, whether or not the coefficients are to be estimated, are consistent 
  with what is required to specify a linear differential equation.
  If certain values are not supplied, default values are assigned.  Two new named 
  members are also supplied for each variable.  Member \code{nallXterm} is the number 
  of terms in the right side homogeneous portion of the equation, which
  involve a derivative of a variable in the system multiplied by a coefficient function.  
  Member \code{nallFterm} is the number of forcing terms in the equation.  
  These two new members are essential for function \code{Data2LD}, and,
  although \code{modelCheck} is also invoked inside \code{Data2LD} therefore a use of 
  \code{modelCheck} is essential prior to using this function.
}
\value{
  A named list object with the same structure as the argument along with added 
  default values and members \code{nallXterm} and \code{nallFterm}.
}
\references{
    J. O. Ramsay and G. Hooker (2017) \emph{Dynamic Data Analysis}. Springer.
}
\seealso{
  \code{\link{make.coef}},
  \code{\link{make.Xterm}},
  \code{\link{make.Fterm}},
  \code{\link{make.variable}},
  \code{\link{coefCheck}},
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
