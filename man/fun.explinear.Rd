\name{fun.explinear}
\alias{fun.explinear}
\title{
A function to evaluate the exponential of a functional data object.
}
\description{
  Coefficient functions in linear differential systems are often required to be strictly 
  positive.  This can be achieved by defining the coefficient function to be the 
  exponential of a functional data object, and using this as the user-defined coefficient
  function \code{fun$fd} that evaluates the coefficient.
}
\usage{
 fun.explinear(tvec, bvec, Bbasisobj)
}
\arguments{
  \item{tvec}{A vector of times at which the function is to be evaluated.}
  \item{bvec}{A vector or a single column of coefficient values for the functional data 
  object.}
  \item{Bbasisobj}{A functional basis object.}
}
\details{
  The length of argument \code{bvec} must match the number of basis functions in 
  argument \code{Bbasisobj}.
}
\value{
  A numeric vector of the same length as that of argument \code{tvec} containing 
  the exponentials of the values of the functional basis object.
}
\references{
    J. O. Ramsay and G. Hooker (2017) \emph{Dynamic Data Analysis}. Springer.
}
\seealso{
  \code{\link{fun.explinear}}
}
