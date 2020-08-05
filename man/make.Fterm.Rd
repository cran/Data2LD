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
make.Fterm(funobj, parvec, estimate, Ufd, factor=1)
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
  \code{\link{checkModel}},
  \code{\link{make.Variable}},
  \code{\link{make.Xterm}},
  \code{\link{printModel}}
}
\examples{
#  Here is the code that sets up the single forcing term for the 
#  refinery data examnple:
TimeData <- RefineryData[,1]
TrayData <- RefineryData[,2]
ValvData <- RefineryData[,3]
conbasis <- create.constant.basis(c(0,193))
alphafdPar <- fdPar(conbasis)
Valvbreaks <- c(0,67,193)
Valvnbasis <- 2
Valvnorder <- 1
Valvbasis  <- create.bspline.basis(c(0,193), Valvnbasis, Valvnorder, Valvbreaks)
Valvfd     <- smooth.basis(TimeData, ValvData, Valvbasis)$fd
#  Fterm Fields:    funobj      parvec  estimate  Ufd      factor
FTerm <- make.Fterm(alphafdPar, 1.0,    TRUE,     Valvfd,  1)
FList <- vector("list", 1)
FList[[1]] <- FTerm
#
#  For other examples of the use of this function, 
#  see these examples in the description of the function \code{make.Variable}.
#
#  The single forcing term of the head impact data
#
#  The cosine and constant forcing terms for the average temperature in
#  Montreal.
#
#  The 20-step forcing function and its coefficient for the "fda" script data.
#
#  The set point forcing term for control variable of the cruise control data.
#
}
