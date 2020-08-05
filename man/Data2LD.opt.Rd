\name{Data2LD.opt}
\alias{Data2LD.opt}
\title{
  Optimize the mean of squared errors data-fitting criterion for a system of linear 
  differential equations.
}
\description{
  This function calls function \code{Data2LD} in each iteration of a quasi-Newton 
  type optimization algorithm, allowing for a set of linear restrictions on the 
  parameter vector.
}  
\usage{
Data2LD.opt(yList, XbasisList, modelList, rhoMat, 
            convcrit = 1e-4, iterlim = 20, dbglev = 1, parMap = diag(rep(1,nparam)))
}
\arguments{
  \item{yList}{A list of length NVAR.  Each list contains in turn a struct object with fields:
    \describe{
      \item{argvals}{is a vector of length \eqn{n_i} of observation times}
      \item{y}{is a matrix with \eqn{n_i} rows and \emph{N} columns.
      The number of columns must be the same for all variables,
      except that, if a list is empty, that variable is taken to be not observed.}
    }
  }
  \item{XbasisList}{A list array of length \emph{d}.  Each member contains in turn
   a functional data object or a functional basis object.}
  \item{modelList}{A list each member of which specifies one of linear differential 
  equations in the system.  These can be    
   constructed using function \code{make.Variable}.  See the help file for that 
   function for further details.}
  \item{rhoMat}{A matrix with the number of rows equal to the number of values of the 
  smoothing parameter per variable to be used in the optimization, and number of columns
  equal to the number of variables.  
  Each entry in the matrix is a value rho in the interval [0,1)  
  so that 0 <= rho < 1.  
  For each variable ivar in the system of equations and each optimization
  iopt, the data are weighted by rho(iopt,ivar) and the
  roughness penalty by 1-rho(iopt,ivar).
  It is expected that the values of rho within each column will be 
  in ascending order, and the estimated parameters for each row
  are passed along as initial values to be used for the optimization
  defined by the values of rho in the next row.}
  \item{convcrit}{A one- or two-part convergence criterion.  The first part is 
  for testing the convergence of the criterion values, and the second part, if 
  present in a vector of length 2, is for testing the norm of the gradient.}
  \item{iterlim}{Maximum number of iterations allowed.}
  \item{dbglev}{Control the output on each iteration.  If 0, no output is displayed.  
  If 1, the criterion and gradient norm at each iteration.  
  If 2, the stepsize, slope and criterion at each step in the line search.
  if 3, displays, in addition to 2, a plot of the criterion values and slope after each 
  line search followed by a pause.}
  \item{parMap}{A rectangular matrix with number of rows equal to the
  number of parameters to be estimated, and number of columns equal to the number of 
  parameters less the number of linear constraints on the estimated parameters.  
  The columns of \code{parMap} must be orthnormal so that the crossproduct is an identity 
  matrix. The crossproduct of \code{t(parMap)} and \code{THETA} maps unconstrained 
  parameters and the corresponding gradient into  constrained parameter space.
  \code{parMap} will usually be set up using the full QR decomposition of a linear 
  constraint coefficient matrix \code{t(A)} where the constraints are of the form 
  A P = B, A and B being known matrices.  An example of such a constraint
  that arises often is one where two estimated coefficients are constrained to be equal.  
  For example, if a variable X involved in an equation the form a(x - x.0), where x.0
  is a fixed set point or defined target level for variable X, then this would be set 
  up as a.1 x + a.2 x.0, where coefficients a.1 and a.2 are constrained to be equal in
  magnitude but opposite in sign, or a.1 + a.2 = 0. }
}                                      
\value{A named list object containing these results of the analysis:\
  \itemize{
    \item{theta.opt:} {A vector containing the optized parameter values.}
    \item{modelList.opt:} {A list object defining the  optimized linear equation
    system.  This may be displayed using function \code{printModel}.}
  }
}
\references{
    J. O. Ramsay and G. Hooker (2017) \emph{Dynamic Data Analysis}. Springer.
}
\seealso{
  \code{\link{Data2LD}}
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
#  Evaluate the fit to the data given the initial parameter estimates (0 and 0)
#  This also initializes the four-way tensors so that they are not re-computed
#  for subsquent analyses.
rhoVec <- 0.5
Data2LDList <- Data2LD(TrayDataList, TrayBasisList, TrayModelList, rhoVec)
MSE  <- Data2LDList$MSE    # Mean squared error for fit to data
DMSE <- Data2LDList$DpMSE  #  gradient with respect to parameter values
#  These parameters control the optimization.
dbglev   <-  1    #  debugging level
iterlim  <- 50    #  maximum number of iterations
convrg   <- c(1e-8, 1e-4)  #  convergence criterion
#  This sets up an increasing sequence of rho values
gammavec <- 0:7
rhoMat   <- as.matrix(exp(gammavec)/(1+exp(gammavec)),length(gammavec),1)
nrho     <- length(rhoMat)
dfesave  <- matrix(0,nrho,1)
gcvsave  <- matrix(0,nrho,1)
MSEsave  <- matrix(0,nrho,1)
thesave  <- matrix(0,nrho,nparam)
#  This initializes the list object containing coefficient estimates
TrayModelList.opt <- TrayModelList
# Loop through values of rho.  
# for test purposes, only the first value rho = 0.5 is used here
for (irho in 1:1) {
  rhoi <- rhoMat[irho,]
  print(paste("rho <- ",round(rhoi,5)))
  Data2LDResult <- 
    Data2LD.opt(TrayDataList, TrayBasisList, TrayModelList.opt, 
                rhoi, convrg, iterlim, dbglev)
  theta.opti <- Data2LDResult$thetastore
  TrayModelList.opt <- modelVec2List(TrayModelList.opt, theta.opti)
  Data2LDList <- Data2LD(TrayDataList, TrayBasisList, TrayModelList.opt, rhoi)
  MSE       <- Data2LDList$MSE 
  df        <- Data2LDList$df
  gcv       <- Data2LDList$gcv 
  ISE       <- Data2LDList$ISE 
  Var.theta <- Data2LDList$Var.theta
  thesave[irho,] <- theta.opti
  dfesave[irho]  <- df
  gcvsave[irho]  <- gcv
  MSEsave[irho]  <- MSE
}
}
