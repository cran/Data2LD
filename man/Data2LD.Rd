\name{Data2LD}
\alias{Data2LD}
\title{
  Evaluate the mean of squared errors data-fitting criterion and its gradient for a 
  system of linear differential equations.
}
\description{
Data2LD ... stands for "Data to Linear Dynamics."
The function approximates the data in argument YLIST by one or smooth 
functions \eqn{latex}{x_i, i=1,\ldots,d.}.  This approximation is defined by a set
of \emph{d} linear differential or algebraic equations defined by a set of
parameters some of which require estimation from the data. 
The approximation minimizes the mean of squared residuals, expressed as 
follows:
\deqn{latex}{H(\theta) = \sum_i w_i \sum_j^{n_i} 
             \sum_{\ell}^N x [y_{ij \ell} - x_i(t_j)]^2}
where:
\itemize{
  \item \eqn{latex}{i \in (1,\ldots,d)} indexes equations in a system differential and 
  algebraic equations.
  \item \eqn{latex}{j = 1,\ldots,n_i} indexes times of observation of a variable.
  \item \eqn{latex}{\ell = 1,\ldots,N} indexes replications of observations.
  \item \eqn{latex}{w_i} is a positive fixed real number weighting the contribution of 
  a variable to the fitting criterion.
}
Only a subset of the variables may be observed, so that not all values of index i are 
actually used.  The number and location of times of observation \eqn{latex}{t_j}  can 
vary from one observed variable to another.
The fitting functions \eqn{latex}{x_i(t)} in turn, here assumed to be defined over
the interval \emph{[0,T]}, are defined by basis function expansions:
\deqn{x_i(t_j) = \sum_k^{K_i} c_{ik}(\theta|\rho) \phi_{ik}(t_j)}
where \eqn{latex}{\phi_{ik}(t)} is the \emph{k}th function in a system of 
\eqn{latex}{K_i} basis functions used to approximate the ith variable, and 
\eqn{latex}{c_{ik}} is the coefficientin the expansion of \eqn{latex}{\phi_{ik}(t)}.
The number of \eqn{latex}{K_i} basis functions and the type of basis function system
can vary from one variable to another.  This information is contained in argument 
\code{XbasisList} described below.
The coefficients \eqn{latex}{c_{ik}(\theta|\rho)} defining the smoothing functions  
are functions of the unknown parameters in vector \eqn{latex}{\theta} that define
the differential and algebraic equations and that require estimation
from the data.  The smoothing parameter \eqn{latex}{\rho} is a value in the interval
\eqn{[0,1)}.  The coefficient functions \eqn{latex}{c_{ik}(\theta|\rho)} minimize the 
inner least squares criterion, expressed here for simplicity for a single  
variable or equation:
\deqn{latex}{J(c|\theta) = (1-\rho) \sum_j^n [y_j - c_k^K \phi_k(t_j)]^2/n +
                \rho \int_0^T {L(\theta)x(t)]^2 dt/T}}
Linear differential operator \eqn{latex}{L(\theta)} is a linear differential or
algebraic equation rewritten as an operator by subtracting the right
side of the equation from the left, so that a solution x of the 
equation satisfies  \eqn{latex}{Lx = 0}.
Each L in a system of equation depends on one or more parameters that 
need to be estimated and are contained in parameter vector \eqn{latex}{\theta}.
The linear differential equation is of the form
    \deqn{latex}{D^m x_i = sum_k^d sum_j^{m_k} beta_{kj}(t) D^{j-1} x_k + 
              sum_f^{F_i} \alpha_{fi}(t) u_{f,i}, 
    i=1,...,d,  f=1,...,F_i}
where each coefficient is expanded in terms of its own number of
B-spline basis functions:
    \deqn{latex}{\beta_{ij}(t)  = \sum_k^{K_{ij}} b_{ijk} \phi_{ijk}(t),}
    \deqn{latex}{\alpha_{fi}(t) = \sum_k^{K_{fi}} a_{fik} \psi_{fik}(t)}
As smoothing parameter \eqn{latex}{\rho} increases toward its upper limit of 1,
the second roughness penalty term in J is more and more emphasized
and the fit to the data less and less emphasized, so that \eqn{latex}{x_i} is
required to approach a solution to its respective equation.
The highest order of derivative can vary from one equation to another,
and in particular can be larger than 1.
Each right side may contain contributions from all variables in the
equation.  Morever, these contributions may be from all permissible
derivatives of each equation, wher "permissible" means up to one less 
than the highest order in the variable.
In addition, each equation may be forced by a number of forcing terms
that can vary from equation to equation.
This version approximates the integrals in the penalty terms by using 
inprod_basis to compute the cross-product matrices for the  
\eqn{latex}{\beta}-coefficient basis functions and the corresponding derivative of 
the x-basis functions,and the cross-product matrices for the 
\eqn{latex}{\alpha}-coefficients and the corresponding U functions.  
These are computed upon the first call to D2LD, and then retained 
for subsequent calls by using the persistent command.  This technique
for speeding up computaton is called memoization.  
The structure of the model is defined in list \code{modelList}, which is
described below.
This version disassociates coefficient functions from equation 
definitions to allow some coefficients to be used repeatedly and for
both homogeneous and forcing terms.  It requires an extra argument
COEFLIST that contains the coefficients and the position of their
coefficient vectors in vector \eqn{latex}{\theta}.
}
\usage{
Data2LD(yList, XbasisList, modelList, rhoVec = 0.5*rep(1,nvar), 
        summary=TRUE)
}
\arguments{
  \item{yList}{A list of length NVAR.  Each list contains  a list object with these 
  fields:
    \itemize{
     \item{argvals:}{ a vector of length \eqn{latex}{n_i} of observation times}
     \item{y:}{ a matrix with \eqn{latex}{n_i} rows and \emph{N} columns.
      The number of columns must be the same for all variables,
      except that, if a list is empty, that variable is taken to be not observed.}
    }
  }
  \item{XbasisList}{A list array of length \emph{d}.  Each member contains in turn
   a functional data object or a functional basis object.}
  \item{modelList}{A list each member of which specifies one of linear differential 
  equations in the system.  These can be constructed using function \code{make.Variable}.  
  See the help file for that function for further details.}
  \item{rhoVec}{A real number in the half-open interval \eqn{[0,T)]}.  
  Its value determines the relative emphasis on fitting the data versus being a 
  solution of the differential equation.  The smaller the value, the more the emphasis 
  on data fitting.}
  \item{summary}{A logical constant.  If TRUE, a variety of summeary measures are 
  computed, including degrees of freedom (\code{df}), the generalized cross-validation 
  index (\code{gcv)} and an estimate of the sampling variance-covariance matrix.  
  However, these can in some circumstances be computationally demanding, and 
  \code{summary=FALSE} limits the output to the fitting criterion, gradient and 
  hessian matrix, which can greatly speed up computation.  Function \code{Data2LD} is 
  called in this way within function \code{Data2LD.opt}, for example.}
}
\value{A named list object containing these results of the analysis:
  \item{MSE}{The weighted mean squared errors computed over the variables with data.}
  \item{DpMSE}{The gradient of the objective function MSE with respect to the estimated 
  parameters.}
  \item{D2ppMSE}{A square symmetric matrx of that contains the second partial 
  derivatives of the objective function with respect to the parameter being estimated.}
  \item{XfdParList}{A list of length \emph{d} containing functional parameter
  objects for the estimated functions \eqn{latex}{x_i(t)}.}
  \item{df}{An equivalent degrees of freedom value 
  \eqn{latex}{df = trace(2 YM - YM*YM)} where \emph{YM} is the matrix \code{y2cMap} 
  described below.}
  \item{gcv}{The generalized cross-validation measure GCV.  The value of
  \eqn{latex}{\rho} corresponding to the minimim of GCV across values of
  smoothing parameter \eqn{latex}{\rho} is often chose for an automatic 
  data-driven level of smoothing.}
  \item{ISE}{The sum across variables of the integrated squared value of
  the differential operator.  This value multiplied by \eqn{latex}{\rho}
  and divided by \emph{T}, the width of the domain, is the second
  term the objective function.}
  \item{Var.theta}{A square symmetric matrix containing estimates of the sampling 
  variances and covariances for the estimated parameter values.}
  \item{Rmat}{A square symmetric matrix associated with the homogeneous portions of 
  the equations that has many useful applications, including constructions of more 
  compact and meaningful basis function expansions for the \eqn{latex}{x_i}'s.}
  \item{Smat}{A matrix containing information on the covariation of the basis functions 
  and the forcing functions.}
}
\references{
    J. O. Ramsay and G. Hooker (2017) \emph{Dynamic Data Analysis}. Springer.
}
\seealso{
  \code{\link{checkModel}},
  \code{\link{make.Variable}},
  \code{\link{make.Xterm}},
  \code{\link{make.Fterm}}
  \code{\link{checkModel}}
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
# An evaluation of the criterion at the initial values
rhoVec <- 0.5  #  light smoothing
Data2LDResult <- Data2LD(ImpactDataList, ImpactBasisList, ImpactModelList, rhoVec)
MSE        <- Data2LDResult$MSE        # Mean squared error for fit to data
DpMSE      <- Data2LDResult$DpMSE      #  gradient with respect to parameter values
D2ppMSE    <- Data2LDResult$D2ppMSE    #  Hessian matrix
}
