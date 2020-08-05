\name{make.Variable}
\alias{make.Variable}
\title{
  Check a list of linear differential equation specifications.
}
\description{
  A system of linear differential equations is defined by a list of length equal 
  to the number of variables in the system.  Each member of this list defines a 
  single linear differential equation.  Within each equation there are typically 
  one or more terms involving a coefficient function multiplying either a derivative 
  of one of the variables, or multiplying a forcing function.  This function sets up 
  a list object that specifies the structure of a single coefficient function.
}
\usage{
make.Variable(name = "", order = 1, XList = NULL, FList = NULL, weight=1)
}
\arguments{
  \item{name}{A string to be used as a name for the equation or variable.  
  By default, it is "xn" where "n" is the position in the list.}
  \item{order}{A nonnegative integer specifying the order of derivative on the left side
  of the equation.}
  \item{XList}{A list each member of which specifies one of the terms in the homogeneous 
  portion of the right side of the equation.  These can be constructed using function 
  \code{make.Xterm}.}
  \item{FList}{A list each member of which specifies one of the forcing terms in the 
  forcing portion of the right side of the equation.  These can be constructed using 
  function \code{make.Fterm}.}
  \item{weight}{A positive weight that can be applied when variables are differently
  weighted.}
}
\details{
  This function checks that all supplied terms conform to what is required.  
  Use of the functions \code{make.Xterm} and \code{make.Fterm} is recommended.
}
\value{
  A named list object defining a linear differential equation.
}
\references{
    J. O. Ramsay and G. Hooker (2017) \emph{Dynamic Data Analysis}. Springer.
}
\seealso{
  \code{\link{make.Xterm}},
  \code{\link{make.Fterm}},
  \code{\link{printModel}}
}
\examples{
#
#  Example 1:  The refinery data
#
N <- 194
TimeData <- RefineryData[,1]
TrayData <- RefineryData[,2]
ValvData <- RefineryData[,3]
#  Define the List array containing the tray data
TrayList <- list(argvals=RefineryData[,1], y=RefineryData[,2])
TrayDataList  <- vector("list",1)
TrayDataList[[1]] <- TrayList
# construct the basis object for tray variable
#  Order 5 spline basis with four knots at the 67
#   to allow discontinuity in the first derivative
#  and 15 knots between 67 and 193
Traynorder <- 5
Traybreaks <- c(0, rep(67,3), seq(67, 193, len=15))
Traynbasis <- 22
TrayBasis  <- create.bspline.basis(c(0,193), Traynbasis, Traynorder, 
                                   Traybreaks)
#  Set up the basis list for the tray variable
TrayBasisList    <- vector("list",1)
TrayBasisList[[1]] <- TrayBasis
#  Both alpha and beta coefficient functions will be constant.
#  Define the constant basis
conbasis <- create.constant.basis(c(0,193))
betafdPar  <- fdPar(conbasis)
alphafdPar <- fdPar(conbasis)
#  Construct a step basis (order 1) and valve level
Valvbreaks <- c(0,67,193)
Valvnbasis <- 2
Valvnorder <- 1
Valvbasis  <- create.bspline.basis(c(0,193), Valvnbasis, Valvnorder, Valvbreaks)
#  smooth the valve data
Valvfd <- smooth.basis(TimeData, ValvData, Valvbasis)$fd
#  Set up the model list for the tray variable
#  Define single homogeneous term
#  Xterm Fields:    funobj     parvec  estimate  variable deriv. factor
XTerm <- make.Xterm(betafdPar, 0.04,  TRUE,     1,       0,     -1)
XList <- vector("list", 1)
XList[[1]] <- XTerm
#  Define the single forcing term
#  Fterm Fields:    funobj      parvec  estimate  Ufd      factor
FTerm <- make.Fterm(alphafdPar, 1.0,    TRUE,     Valvfd,  1)
FList <- vector("list", 1)
FList[[1]] <- FTerm
#  Define the single differential equation in the model
TrayVariable <- make.Variable(XList=XList, FList=FList, 
                              name="Tray level", order=1)
#  Set up the model List array
TrayModelList <- vector("list",1)
TrayModelList[[1]] <- TrayVariable
#  Run a check on TrayVariableList and make the modelList object
checkModelList <- checkModel(TrayBasisList, TrayModelList, summarywrd=TRUE)
TrayModelList  <- checkModelList$model
nparam         <- checkModelList$nparam
#
#  Example 2:  The head impact data
#
ImpactTime <- HeadImpactData[,2]  #  time in milliseconds
ImpactData <- HeadImpactData[,3]*100  #  acceleration in millimeters/millisecond^2
ImpactRng  <- c(0,60) # Define range time for estimated model
#  Define the List array containing the data
ImpactDataList1 <- list(argvals=ImpactTime, y=ImpactData)
ImpactDataList <- vector("list",1)
ImpactDataList[[1]] <- ImpactDataList1
# Define a unit pulse function located at times 14-15
Pulsebasis <- create.bspline.basis(ImpactRng, 3, 1, c(0,14,15,60))
Pulsefd    <- fd(matrix(c(0,1,0),3,1),Pulsebasis)
#  Construct basis for output x(t), 
#  multiple knots at times 14 and 15.
#  Order 6 spline basis, with three knots at the impact point and 
#  three knots at impact + delta to permit discontinuous first 
#  derivatives at these points
knots     <- c(0,14,14,14,15,15,seq(15,60,len=11))
norder    <- 6
nbasis    <- 21
ImpactBasis <- create.bspline.basis(ImpactRng,nbasis,norder,knots)
#  plot the basis
ImpactTimefine <- seq(0,60,len=201)
ImpactBasisfine <- eval.basis(ImpactTimefine, ImpactBasis)
#  set up basis list
ImpactBasisList <- vector("list",1)
ImpactBasisList[[1]] <- ImpactBasis
# Set up the constant basis, make three coefficients and check them
conbasis <- create.constant.basis(ImpactRng)
# Define the terms in the second order linear equation
# Define the two terms in the homogeneous part of the equation
#  Xterm Fields:     funobj     parvec  estimate  variable deriv. factor
Xterm1 <- make.Xterm(conbasis,  1,      TRUE,     1,       0,     -1)
Xterm2 <- make.Xterm(conbasis,  1,      TRUE,     1,       1,     -1)
# Set up the XList vector of length two
XList <- vector("list",2)
XList[[1]] <- Xterm1
XList[[2]] <- Xterm2
# Define the forcing term
#  Fterm Fields:    funobj    parvec  estimate  Ufd      factor
Fterm <- make.Fterm(conbasis, 1.0,   TRUE,     Pulsefd, 1)
# Set up the forcing term list of length one
FList      <- vector("list",1) 
FList[[1]] <- Fterm
#  Define the single differential equation in the model
ImpactVariable <- make.Variable(name="Impact", order=2, XList=XList, FList=FList)
#  Define list of length one containing the equation definition
ImpactModelList=vector("list",1)
ImpactModelList[[1]] <- ImpactVariable
#  Check the object for internal consistency and 
#  set up the model object
#  Run a check on TrayVariableList and make the modelList object
checkModelList  <- checkModel(ImpactBasisList, ImpactModelList, summarywrd=TRUE)
ImpactModelList <- checkModelList$model
nparam          <- checkModelList$nparam
#
#  Example 3:  The X coordinate of the fda script data
#
fdascriptX <- as.matrix(fdascript[,1,1])
fdascriptY <- as.matrix(fdascript[,1,2])
#  Define the observation times in 100ths of a second
centisec <- seq(0,2.3,len=1401)*100
fdarange <- range(centisec)
fda.dataList = vector("list", 2)
fda.dataList[[1]]$argvals <- centisec
fda.dataList[[1]]$y       <- fdascriptX
fda.dataList[[2]]$argvals <- centisec
fda.dataList[[2]]$y       <- fdascriptY
#  --------------------------------------------------------------------
#   Set up the basis object for representing the X and Y coordinate 
#   functions
#  --------------------------------------------------------------------
#  The basis functions for each coordinate will be order 6 B-splines,
#  32 basis functions per second, 4 per 8 herz cycle
#  2.3 seconds 2.3*32=73.6.
#  Order 6 is used because we want to estimate a smooth second
#  derivative.
#  This implies norder + no. of interior knots <- 74 - 1 + 6 <- 79 
#  basis functions.  
nbasis   <- 79
norder   <-  6 #  order 6 for a smooth 2nd deriv.
fdabasis <- create.bspline.basis(fdarange, nbasis, norder)
fdaBasisList <- vector("list",2)
fdaBasisList[[1]] <- fdabasis
fdaBasisList[[2]] <- fdabasis
#  --------------------------------------------------------------------
#   Use the basis to estimate the curve functions and their 
#   second derivatives
#  --------------------------------------------------------------------
ResultX <- smooth.basis(centisec, fdascriptX, fdabasis)
fdascriptfdX <- ResultX$fd
D2fdascriptX <- eval.fd(centisec, fdascriptfdX, 2)
ResultY <- smooth.basis(centisec, fdascriptY, fdabasis)
fdascriptfdY <- ResultY$fd
D2fdascriptY <- eval.fd(centisec, fdascriptfdY, 2)
#  --------------------------------------------------------------------
#   Set up the basis objects for representing the coefficient 
#   functions
#  --------------------------------------------------------------------
#  Define a constant basis and function over the 230 centisecond range 
conbasis <- create.constant.basis(fdarange)
confd    <- fd(1,conbasis)
confdPar <- fdPar(confd)
#  Define the order one Bspline basis for the coefficient
#  of the forcing coefficients
nevent    <- 20
nAorder   <- 1
nAbasis   <- nevent
nAknots   <- nAbasis + 1
Aknots    <- seq(fdarange[1],fdarange[2],len=nAknots)
Abasisobj <- create.bspline.basis(fdarange,nAbasis,nAorder,Aknots)
#  --------------------------------------------------------------------
#  Set up single homogeneous terms in the homogeneous portion of the 
#  equations  D^2x = -beta_x x and D^2y = -beta_y y.  
#  Each term involves a scalar constant multiplier -1.
#  --------------------------------------------------------------------
Xterm <- make.Xterm(confdPar, 0.04, TRUE, variable=1, derivative=0, factor=-1)
Yterm <- make.Xterm(confdPar, 0.04, TRUE, variable=2, derivative=0, factor=-1)
#  Set up to list arrays to hold these term specifications
XListX <- vector("list",1)
XListX[[1]] <- Xterm
XListY <- vector("list",1)
XListY[[1]] <- Yterm
#  --------------------------------------------------------------------
#  Set up coefficients for the forcing terms 
#  \alpha_x(t)*1 and \alpha_y(t)*1.  
#  The forcing functions are the unit constant function.
#  --------------------------------------------------------------------
parvecA <- matrix(0,nAbasis,1)
FtermX <- make.Fterm(Abasisobj, parvecA, TRUE, Ufd=confd, factor=1)
FtermY <- make.Fterm(Abasisobj, parvecA, TRUE, Ufd=confd, factor=1)
#  Set up list arrays to hold forcing terms specs
FListX <- vector("list", 1)
FListX[[1]] <- FtermX
FListY <- vector("list", 1)
FListY[[1]] <- FtermY
#  --------------------------------------------------------------------
#  make variables for X and Y coordinates and system object fdaVariableList
#  --------------------------------------------------------------------
Xvariable <- make.Variable(name="X", order=2, XList=XListX, FList=FListX)
Yvariable <- make.Variable(name="Y", order=2, XList=XListY, FList=FListY)
#  List array for the whole system
fdaVariableList <- vector("list",2)
fdaVariableList[[1]] <- Xvariable
fdaVariableList[[2]] <- Yvariable
#  check the system specification for consistency,  This is a mandatory
#  step because objects are set up in the output list that are needed
#  in the analysis.  A summary is displayed.
checkModelList <- checkModel(fdaBasisList, fdaVariableList, TRUE)
fda.modelList <- checkModelList$modelList
nparam        <- checkModelList$nparam
#
#  Example 4:  Average temperature in Montreal
#
daytime73 <- seq(2.5,362.5,len=73)
dayrange  <- c(0,365)
#  set up block averages for Canadian temperature data,
#  available from the fda package
tempav <- CanadianWeather$dailyAv[,,"Temperature.C"]
tempav73 <- matrix(0,73,35)
m2 <- 0
for (  i in 1 : 73 ) {
  m1 <- m2 + 1
  m2 <- m2 + 5
  tempavi <- apply(tempav[m1:m2,1:35],2,mean)
  tempav73[i,] <- tempavi
} 
#  center the time of the year on winter
winterind73  <- c( 37:73,1: 36)
tempav73 <- tempav73[winterind73,]
#  select the data for Montreal
station <- 12 
#  set up TempDataList
TempDataList1 <- list(argvals=as.matrix(daytime73), y=as.matrix(tempav73[,station]))
TempDataList  <- vector("list",1)
TempDataList[[1]]    <- TempDataList1
#  basis for 5-day block averages
norder    <- 5
daybreaks <- seq(0,365,5)
nbreaks   <- length(daybreaks)
nbasis    <- norder + nbreaks - 2
daybasis  <- create.bspline.basis(dayrange, nbasis)
TempBasisList <- vector("list",1)
TempBasisList[[1]] <- daybasis
#  Define the two forcing functions:
#  constant function Ufd for constant forcing
Uconbasis <- create.constant.basis(dayrange)
Uconfd    <- fd(1, Uconbasis)
Uconvec   <- matrix(1,73,1)
#  cosine function for solar radiation forcing
uvec      <- -cos((2*pi/365)*(daytime73+10+182))
Ucosbasis <- create.fourier.basis(dayrange, 3)
Ucosfd    <- smooth.basis(daytime73, uvec, Ucosbasis)$fd
#  Define three basis functional objects for coefficients
#  A fourier basis for defining the positive homogeneous 
#  coefficient
nWbasis   <- 7
Wbasisobj <- create.fourier.basis(dayrange, nWbasis)
#  constant forcing coefficient for constant forcing
nAbasisC   <- 1
AbasisobjC <- create.constant.basis(dayrange)
#  fourier coefficint function for radiative forcing
nAbasisF   <- 1
AbasisobjF <- create.constant.basis(dayrange)
#  Set up the list object for the positive coefficient for 
#  the homogeneous term.  
linfun <- list(fd=fun.explinear, Dfd=fun.Dexplinear, more=Wbasisobj)
#  Define the variable for temperature
#  define homogeneous term and list container
estimate = rep(TRUE,7)
estimate[7] <- FALSE
parvec <- matrix(0,7,1)
# parvec[7] <- 3
Xterm <- make.Xterm(funobj=linfun, parvec=parvec, estimate=estimate,
                    variable=1, derivative=0, factor= -1)
XList <- vector("list", 1)
XList[[1]] <- Xterm
#  define two forcing terms
Fterm1 <- make.Fterm(funobj=AbasisobjC, parvec=1.0, estimate=TRUE,
                     Ufd=Uconfd, factor=1)
Fterm2 <- make.Fterm(funobj=AbasisobjF, parvec=1.0, estimate=TRUE,
                     Ufd=Ucosfd, factor=1)
#  forcing term container
FList <- vector("list", 2)
FList[[1]] <- Fterm1
FList[[2]] <- Fterm2
#  set up variable list object
TempVariableList1 <- make.Variable(name="Temperature", order=1, XList=XList, FList=FList)
# set up TempVariableList
TempModelList <- vector("list",1)
TempModelList[[1]] <- TempVariableList1
# define the model list object
TempModelcheck <- checkModel(TempBasisList, TempModelList, TRUE)
TempModelList <- TempModelcheck$modelList
nparam <- TempModelcheck$nparam
}
