\name{getForceTerm}
\alias{getForceTerm}
\title{
Retrieve the specification of a forcing term in a linear differential equation.
}
\description{
A forcing term in a linear differential equation is a term on the right side of the equation that consists of an external variable that is multiplied by a coefficient function, and, possibly, is also multiplied by a fixed scalar constant.  The coefficient function may be constant, or may vary with time itself, or with any external variable.  It may not, however, depend on the value of either any of the variables in the system or on any derivative of any variable.  Function \code{getHomoTerm} retrieves and returns as a named list the objects that specify the homogeneous term. The named list is detailed below.
}
\usage{
getForceTerm(FtermList, coefList)
}
\arguments{
  \item{FtermList}{ A named list containing the values defining the term, with fields and properties as follows:
  \enumerate{
    \item{\code{ncoef}:  the the index in the coefficient list object of the variable in the term.}, 
    \item{\code{factor}: a fixed scalar multiplier, such as -1, for the term.  The default is 1.}
    \item{\code{derivative}:  the order of the derivative of the variable in the term.  The default is 0.}
  }
  }
  \item{coefList}{A named list containing specifications all of the coefficient functions used in the model.  These
  specifications can be set up using function \code{make.coef}.}
}
\value{
A named list containing the argument values.
}
\references{
    J. O. Ramsay and G. Hooker (2017) \emph{Dynamic Data Analysis}. Springer.
}
\author{
J. O. Ramsay
}
\seealso{
  \code{\link{coefCheck}},
  \code{\link{make.Coef}},
  \code{\link{make.Fterm}},
  \code{\link{make.Model}}
  \code{\link{make.Variable}},
  \code{\link{make.Xterm}},
}
