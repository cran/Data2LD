
#  The data are analyzed in the Ramsay and Silverman books, and consist of 
#  of the script version of "fda" written by Jim Ramsay.  Because script
#  is cursive, there are two variables, "X" being the motion across the 
#  front of the body, and "Y" being fore-and-aft motion.

#  The two second order linear differential equations analyzed here are:
#  $D^2X(t) = -b_X X(t) + a_X(t)$
#  $D^2Y(t) = -b_Y Y(t) + a_Y(t)$
#  where $b_X$ and $b_Y$ are scalar constants, called the "reaction speed"
#  for the coefficient.  The coefficient functions $a_X(t)$ and $a_Y(t)$
#  are step functions, each having 20 steps.  They can be thught of as
#  coefficient functions multiplying a dummy forcing function U(t) which
#  only takes the value one.
#  Thus, each equation there are 21 parameters to estimate from the data,
#  the coordinate's reaction speed and the 20 heights of the steps.
#  Although there are no shared terms in these two equations, so that
#  they could be estimated indendently, for the purposes of neatness and
#  clarity, we shall treat them as a system of two equations with a 
#  total of 42 parameters to be estimated.
#  Each estimated step function is a spline function of order one.

#  --------------------------------------------------------------------
#        Load the fda script data and define the time values
#  --------------------------------------------------------------------

#  select out the two 1401 by 20 matrices of coordinate values

fdascriptX <- as.matrix(fdascript[,1,1])
fdascriptY <- as.matrix(fdascript[,1,2])

#  Define the observation times in 100ths of a second

centisec <- seq(0,2.3,len=1401)*100
fdarange <- range(centisec)

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

#  Define the order one Bspline basis for the coefficient
#  of the forcing coefficients

nevent    <- 20
nAorder   <- 1
nAbasis   <- nevent
nAknots   <- nAbasis + 1
Aknots    <- seq(fdarange[1],fdarange[2],len=nAknots)
Abasisobj <- create.bspline.basis(fdarange,nAbasis,nAorder,Aknots)

#  --------------------------------------------------------------------
#  Define the two pairs of coefficient functions, save them in a List 
#  object, and check the list object
#  --------------------------------------------------------------------

coef1 <- make.Coef(conbasis,  0.04,                TRUE)
coef2 <- make.Coef(Abasisobj, matrix(1,nAbasis,1), TRUE)
coef3 <- make.Coef(conbasis,  0.04,                TRUE)
coef4 <- make.Coef(Abasisobj, matrix(1,nAbasis,1), TRUE)

# Set up a list array containing the coefficient lists

coefList <- vector("list",4)
coefList[[1]] <- coef1
coefList[[2]] <- coef2
coefList[[3]] <- coef3
coefList[[4]] <- coef4

# Check the coefficients, a mandatory step since coefCheck also
# sets up the number of parameters to estimate.    A summary is displayed.

coefResult <- coefCheck(coefList)
coefList   <- coefResult$coefList
ntheta     <- coefResult$ntheta
print(paste("ntheta = ",ntheta))

#  --------------------------------------------------------------------
#  Set up single homogeneous terms in the homogeneous portion of the 
#  equations  D^2x = -beta_x x and D^2y = -beta_y y.  
#  Each term involves a scalar constant multiplier -1.
#  --------------------------------------------------------------------

Xterm <- make.Xterm(variable=1, derivative=0, ncoef=1, factor=-1)
Yterm <- make.Xterm(variable=2, derivative=0, ncoef=3, factor=-1)

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

FtermX <- make.Fterm(ncoef=2, Ufd=confd, factor=1)
FtermY <- make.Fterm(ncoef=4, Ufd=confd, factor=1)

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

fda.modelList <- make.Model(fdaBasisList, fdaVariableList, coefList)

#  --------------------------------------------------------------------
#  Set up the data lists for X- and Y-coordinates 
#  usiing only the first replication.
#  --------------------------------------------------------------------

fda.dataListX <- list(argvals=centisec, y=fdascriptX[,1])
fda.dataListY <- list(argvals=centisec, y=fdascriptY[,1])

# Set up the data list array containing these two list objects

fda.dataList <- vector("list", 2)
fda.dataList[[1]] <- fda.dataListX
fda.dataList[[2]] <- fda.dataListY

#  --------------------------------------------------------------------
#  Single evaluation in order to set up the 4-way tensors and
#  make the model list object
#  --------------------------------------------------------------------

rhoVec <- matrix(0.5,1,2)
Data2LDResult <- Data2LD(fda.dataList, fdaBasisList, fda.modelList, coefList, rhoVec)

MSE        <- Data2LDResult$MSE        #  Mean squared error for fit to data
DpMSE      <- Data2LDResult$DpMSE      #  gradient with respect to parameter values
D2ppMSE    <- Data2LDResult$D2ppMSE    #  Hessian matrix
XfdParList <- Data2LDResult$XfdParList #  List of fdPar objects for variable values 
df         <- Data2LDResult$df         #  Degrees of freedom for fit
gcv        <- Data2LDResult$gcv        #  Generalized cross-validation coefficient
ISE        <- Data2LDResult$ISE        #  Size of second term, integrated squared error
Var.theta  <- Data2LDResult$Var.theta  #  Estimate sampling variance for parameters

#  --------------------------------------------------------------------
#  Set up sequence of rho values to be used in the optimization
#  and the display level, iteration limit and convergence criterion
#  --------------------------------------------------------------------

gamvec <- 0:7
rhoMat <- matrix(exp(gamvec)/(1+exp(gamvec)),length(gamvec),2)
nrho   <- dim(rhoMat)[1]

#  values controlling optimization

dbglev   <- 1              #  display level
iterlim  <- 50             #  maximum number of iterations
convrg   <- c(1e-6, 1e-4)  #  convergence criterion

#  --------------------------------------------------------------------
#  Set up arrays to hold results over rho and 
#  Initialize coefficient list
#  --------------------------------------------------------------------

dfesave <- matrix(0,nrho,1)
gcvsave <- matrix(0,nrho,1)
MSEsave <- matrix(0,nrho,1)
ISEsave <- matrix(0,nrho,1)
thesave <- matrix(0,nrho,ntheta)

coefList.opt <- coefList

#  --------------------------------------------------------------------
#  Step through rho values, optimizing at each step, and saving
#  results along the way.
#  --------------------------------------------------------------------

for (irho in 1:nrho) {
  rhoVeci <- matrix(rhoMat[irho,],1,2)
  print(paste("rho <- ",round(rhoVeci[1],5)))
  Data2LD.optResult <- Data2LD.opt(fda.dataList, fdaBasisList, fda.modelList, coefList.opt, 
                                   rhoVeci, convrg, iterlim, dbglev)
  theta.opti    <- Data2LD.optResult$thetastore
  coefList.opt  <- modelVec2List(theta.opti, coefList.opt)
  Data2LDResult <- Data2LD(fda.dataList, fdaBasisList, fda.modelList, coefList.opt, rhoVeci)
  thesave[irho,]   <- theta.opti
  dfesave[irho,1]  <- Data2LDResult$df
  gcvsave[irho,1]  <- Data2LDResult$gcv
  MSEsave[irho,1]  <- Data2LDResult$MSE
  ISEsave[irho,1]  <- Data2LDResult$ISE
}

#  --------------------------------------------------------------------
#  Display degrees of freedom and gcv values for each step
#  --------------------------------------------------------------------

print("        rho    df      gcv")
print(round(cbind(rhoMat[,1], dfesave, gcvsave),4))

#  --------------------------------------------------------------------
#  Evaluate the fit for parameter values at the highest rho value
#  and display the numerical results
#  --------------------------------------------------------------------

irho <- 8
theta    <- thesave[irho,]
coefList <- modelVec2List(theta, coefList)
rho       <- rhoMat[irho,]
Data2LDResult <- Data2LD(fda.dataList, fdaBasisList, fda.modelList, coefList, rho)
MSE        <- Data2LDResult$MSE        #  Mean squared error for fit to data
DpMSE      <- Data2LDResult$DpMSE      #  gradient with respect to parameter values
D2ppMSE    <- Data2LDResult$D2ppMSE    #  Hessian matrix
XfdParList <- Data2LDResult$XfdParList #  List of fdPar objects for variable values 
df         <- Data2LDResult$df         #  Degrees of freedom for fit
gcv        <- Data2LDResult$gcv        #  Generalized cross-validation coefficient
ISE        <- Data2LDResult$ISE        #  Size of second term, integrated squared error
Var.theta  <- Data2LDResult$Var.theta  #  Estimate sampling variance for parameters

print(paste("MSE = ", round(MSE,4)))
print(paste("df  = ", round(df, 4)))
print(paste("gcv = ", round(gcv,4)))

print(paste("X reaction speed = ",round(theta[ 1],3)))
print(paste("Y reaction speed = ",round(theta[22],3)))

#  --------------------------------------------------------------------
#  Display the graphical results
#  --------------------------------------------------------------------

fdascriptfdX      <- Data2LDResult$XfdParList[[1]]$fd
fdascriptfdY      <- Data2LDResult$XfdParList[[2]]$fd

fdascriptTimefine <- seq(fdarange[1],fdarange[2],len=201)
fdascriptfineX    <- eval.fd(fdascriptTimefine, fdascriptfdX)
fdascriptfineY    <- eval.fd(fdascriptTimefine, fdascriptfdY)

fdascriptKnotsX    <- eval.fd(Aknots, fdascriptfdX)
fdascriptKnotsY    <- eval.fd(Aknots, fdascriptfdY)

fdascriptKnotsX <- approx(centisec, fdascriptX[,1], Aknots)
fdascriptKnotsY <- approx(centisec, fdascriptY[,1], Aknots)

#  plot fit to the data as a script

par(mfrow=c(1,1))
plot(fdascriptfineX, fdascriptfineY, type="l", 
     xlab="X coordinate", ylab="Y coordinate")
points(fdascriptX[,1], fdascriptY[,1], pch="o")
points(fdascriptKnotsX$y, fdascriptKnotsY$y, pch="X", col=2)

#  plot step heights as a script

indX <-  2:21;
indY <- 23:42;

StepX <- theta[indX]
StepY <- theta[indY]

plot(StepX, StepY, type="b", 
     xlab="X coordinate", ylab="Y coordinate")

#  plot fit to the data in two-panels

par(mfrow=c(2,1))
plot(centisec, fdascriptX[,1], type="p",
     xlab="Time (centisec)", ylab="X-coordinate")
lines(fdascriptTimefine, fdascriptfineX)

plot(centisec, fdascriptY[,1], type="p",
     xlab="Time (centisec)", ylab="Y-coordinate")
lines(fdascriptTimefine, fdascriptfineY)

#  plot fit to acceleration

D2fdascriptfineX  <- eval.fd(fdascriptTimefine, fdascriptfdX, 2)
D2fdascriptfineY  <- eval.fd(fdascriptTimefine, fdascriptfdY, 2)

par(mfrow=c(2,1))
plot(centisec, D2fdascriptX[,1], type="p",
     xlab="Time (centisec)", ylab="D2X-coordinate")
lines(fdascriptTimefine, D2fdascriptfineX)

plot(centisec, D2fdascriptY[,1], type="p",
     xlab="Time (centisec)", ylab="D2Y-coordinate")
lines(fdascriptTimefine, D2fdascriptfineY)


