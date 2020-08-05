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
make.Xterm(funobj, parvec, estimate, variable, derivative=0, factor=1)
}
\arguments{
  \item{funobj}{A specification of the coefficient function.  If \code{funobj}
  is of class \code{basisobj}, \code{fd} or \code{fdPar}, the coefficient function
  is a functional data object for a single function.  If the class is \code{list},
  its members define a function, its derivative, and other information that may
  be needed in the definition.}
  \item{parvec}{The vector of parameters defining coefficient function.}
  \item{estimate}{A vector of the same length as \code{parvec} containing
  logical values defining which parameteres are to be estimated (TRUE) or kept
  fixed (FALSE).  If all values are the same, a single TRUE or FALSE is sufficient.}
  \item{variable}{An integer specifying the variable within which this is a term.}
  \item{derivative}{The order of the derivative of the variable.}
  \item{factor}{A real number that is treated as fixed. For example, it is frequently 
  the case that a variable will appear in two or more places in a system of equations, 
  and sometimes multipled by -1.}   
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
  \code{\link{make.Fterm}},
  \code{\link{make.Variable}},
  \code{\link{printModel}}
}
\examples{
#  Here is the code that sets up the single homogeneous term for the 
#  refinery data examnple:
conbasis  <- create.constant.basis(c(0,193))
betafdPar <- fdPar(conbasis)
#  Xterm Fields:    funobj     parvec  estimate  variable deriv. factor
XTerm <- make.Xterm(betafdPar, 0.04,  TRUE,      1,       0,     -1)
#  Enter this list object in a list of length one.
XList <- vector("list", 1)
XList[[1]] <- XTerm
#  Here is the code that sets up the single homogeneous term involving a 
#  non-functional data function for the average temperature for Montreal.
#  Functions \code{fun.explinear} and \code{fun.Dexplinear} are in the 
#  package, and compute the positive value of the exponential transform
#  of a functional data object and its derivative.  The basis object
#  Set up the list object for the positive coefficient for 
#  the homogeneous term.  
#  \code{Wbasisobj} is also supplied in the list object \code{linfun}.
nWbasis   <- 7
Wbasisobj <- create.fourier.basis(c(0,365), nWbasis)
linfun <- list(fd=fun.explinear, Dfd=fun.Dexplinear, more=Wbasisobj)
estimate = rep(TRUE,7)
estimate[7] <- FALSE
parvec <- matrix(0,7,1)
Xterm <- make.Xterm(funobj=linfun, parvec=parvec, estimate=estimate,
                    variable=1, derivative=0, factor= -1)
XList <- vector("list", 1)
XList[[1]] <- Xterm
#
#  For other examples of the use of this function, 
#  see these examples in the description of the function \code{make.Variable}.
#
#  The single homogeneous term of the head impact data
#
#  The single zero derivative terms in the two second order equations
#  for the "fda" script model.
#
#  The speed and control terms of the cruise control data
#
}
