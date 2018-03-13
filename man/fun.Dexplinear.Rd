\name{fun.Dexplinear}
\alias{fun.Dexplinear}
\title{
  A function to evaluate the derivatives with respect to the coefficients 
  of the exponential of a functional data object.
}
\description{
  Coefficient functions in linear differential systems are often required to be strictly 
  positive.  This can be achieved by defining the coefficient function to be the 
  exponential of a functional data object.  This function evaluates the partial 
  derivatives of the coefficient function with respect to the coefficients in the 
  functional data object.  It is used as the member \code{fun.Dfd} of a user-defined 
  coefficient function.
}
\usage{
  fun.Dexplinear(tvec, bvec, Bbasisobj)
}
\arguments{
  \item{tvec}{A vector of times at which the function is to be evaluated.}
  \item{bvec}{A vector or a single column of coefficient values for the functional 
  data object.}
  \item{Bbasisobj}{A functional basis object.}
}
\details{
  The length of argument \code{bvec} must match the number of basis functions in 
  argument \code{Bbasisobj}.
}
\value{
  A matrix with the number of rows equal to the length of argument \code{tvec} and 
  number of columns equal to the number of basis functions.  
  This matrix contains the partial derivatives with respect to the coefficients 
  defining the functional data object evaluated at the time values in argument 
  \code{tvec}.
}
\references{
    J. O. Ramsay and G. Hooker (2017) \emph{Dynamic Data Analysis}. Springer.
}
\seealso{
  \code{\link{fun.explinear}}
}
