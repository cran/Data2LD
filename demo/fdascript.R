
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

load("~/Documents/R/Data2LD/Data/fdascript.rda")

#  select out the two 1401 by 20 matrices of coordinate values

fdascriptX <- as.matrix(fdascript[,1,1])
fdascriptY <- as.matrix(fdascript[,1,2])

#  Define the observation times in 100ths of a second

centisec <- seq(0,230,len=1401)
fdarange <- range(centisec)

yList = vector("list", 2)
yList[[1]]$argvals <- centisec
yList[[1]]$y       <- fdascriptX
yList[[2]]$argvals <- centisec
yList[[2]]$y       <- fdascriptY

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

XbasisList <- vector("list",2)
XbasisList[[1]] <- fdabasis
XbasisList[[2]] <- fdabasis

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

Xterm <- make.Xterm(confdPar, 0.4, TRUE, variable=1, derivative=0, factor=-1)
Yterm <- make.Xterm(confdPar, 0.4, TRUE, variable=2, derivative=0, factor=-1)

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

parvecA <- matrix(1,nAbasis,1)
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

checkfdaModel <- checkModel(XbasisList, fdaVariableList, TRUE)

fdaModel <- checkfdaModel$modelList
nparam   <- checkfdaModel$nparam

print(paste("Number of parameters =",nparam))

rhoVec <- matrix(0.5,2,1)

Data2LDResult <- Data2LD(yList, XbasisList, fdaModel, rhoVec, summary=TRUE)

MSE0      <- Data2LDResult$MSE        #  Mean squared error for fit to data
DpMSE0    <- Data2LDResult$DpMSE      #  gradient with respect to parameter values
D2ppMSE0  <- Data2LDResult$D2ppMSE    #  Hessian matrix
XfdList   <- Data2LDResult$XfdParList        
df        <- Data2LDResult$df      
gcv       <- Data2LDResult$gcv    
ISE       <- Data2LDResult$ISE        
Var.theta <- Data2LDResult$Var.theta      
Rmat      <- Data2LDResult$Rmat    
Smat      <- Data2LDResult$Smat        
DRarray   <- Data2LDResult$DRarray      
DSmat     <- Data2LDResult$DSmat    

MSE0
DpMSE0

#  --------------------------------------------------------------------
#  Set up sequence of rho values to be used in the optimization
#  and the display level, iteration limit and convergence criterion
#  --------------------------------------------------------------------

gamvec <- 0:7
rhoMat <- t(matrix(exp(gamvec)/(1+exp(gamvec)),8,2))
nrho   <- dim(rhoMat)[2]

#  values controlling optimization

dbglev   <- 1     #  display level
iterlim  <- 50    #  maximum number of iterations
convrg   <- 1e-6  #  convergence criterion

#  --------------------------------------------------------------------
#  Set up arrays to hold results over rho and 
#  Initialize coefficient list
#  --------------------------------------------------------------------

dfesave <- matrix(0,nrho,1)
gcvsave <- matrix(0,nrho,1)
MSEsave <- matrix(0,nrho,1)
ISEsave <- matrix(0,nrho,1)
thesave <- matrix(0,nrho,nparam)

fdaModel.opt <- fdaModel

#  --------------------------------------------------------------------
#  Step through rho values, optimizing at each step, and saving
#  results along the way.
#  --------------------------------------------------------------------

for (irho in 1:nrho) {
  rhoVeci <- matrix(rhoMat[,irho],2,1)
  print(paste("rho <- ",round(rhoVeci[1],5)))
  Result <- Data2LD.opt(yList, XbasisList, fdaModel.opt, rhoVeci,  
                                   convrg, iterlim, dbglev)
  theta.opti       <- Result$thetastore
  fdaModel.opt     <- modelVec2List(fdaModel.opt, theta.opti)
  Data2LDResult    <- Data2LD(yList, XbasisList, fdaModel.opt, rhoVeci, summary=TRUE)
  thesave[irho,]   <- theta.opti
  dfesave[irho,1]  <- Data2LDResult$df
  gcvsave[irho,1]  <- Data2LDResult$gcv
  MSEsave[irho,1]  <- Data2LDResult$MSE
  ISEsave[irho,1]  <- Data2LDResult$ISE
}

#  Each of the optimizations is started from the converged value for
#  previous optimization, and they all converge smoothly and rapidly.
#  If you want see how this can fail, change fdaModel.opt to just
#  fdaModel, so that each iteration starts from the initial parameter
#  values.  You will see a catastrophic failure for the last value of
#  rho.

#  --------------------------------------------------------------------
#  Display degrees of freedom and gcv values for each step
#  --------------------------------------------------------------------

print("        rho    df      gcv")
print(round(cbind(rhoMat[1,], dfesave, gcvsave),4))

#  plot flow across rho values for the two homogeneous coefficients

matplot(rhoMat[1,], thesave[,1:2], type="b", pch="o",lwd=2,
xlab =  'rho', ylabel = 'beta coefficients')
legend(x=0.6, y=0.04, c("X", "Y"), lty=c(1,2), lwd=2, col=c(1,2))

#  plot flow of MSE and ISE

par(mfrow=c(2,1))
plot(rhoMat[1,], MSEsave, type="b", pch="o",lwd=2,
     xlab = "rho", ylab = "MSE")
subplot(2,1,2)
plot(rhoMat[1,], ISEsave, type="b", pch="o",lwd=2,
     xlab = "rho", ylab = "ISE")
     
#  --------------------------------------------------------------------
#  Evaluate the fit for parameter values at the highest rho value
#  and display the numerical results
#  --------------------------------------------------------------------

irho <- 8
theta    <- thesave[irho,]
fdaModel <- modelVec2List(fdaModel, theta)
rho       <- rhoMat[,irho]
Data2LDResult <- Data2LD(yList, XbasisList, fdaModel, rho, summary=TRUE)
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

print(paste("X reaction speed = ",round(theta[1],3)))
print(paste("Y reaction speed = ",round(theta[2],3)))

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

#  plot fit to the data in two-panels

par(mfrow=c(2,1))
plot(centisec, fdascriptX[,1], type="p",
     xlab="Time (centisec)", ylab="X-coordinate")
lines(fdascriptTimefine, fdascriptfineX)

plot(centisec, fdascriptY[,1], type="p",
     xlab="Time (centisec)", ylab="Y-coordinate")
lines(fdascriptTimefine, fdascriptfineY)

#  plot the step heights

fdascriptKnotsX    = eval.fd(Aknots, fdascriptfdX);
fdascriptKnotsY    = eval.fd(Aknots, fdascriptfdY);

par(mfrow=c(2,1))

plot( c(  0,  0), c(0,fdascriptKnotsX[1]), type="l", lwd=2,
      xlim=c(0,230), ylim=c(-40,40), xlab = "", ylab = "X")
lines(c(230,230), c(fdascriptKnotsX[20],0), type="l", lwd=2) 
lines(c(  0,230), c(0,0), type="l", lwd = 2, lty=2)
for (i in 1:20) {
  lines(c(Aknots[i],Aknots[i+1]), c(fdascriptKnotsX[i],fdascriptKnotsX[i]), 
       type="l", lwd=2)
}
for (i in 1:19) {
  lines(c(Aknots[i+1],Aknots[i+1]), 
        c(fdascriptKnotsX[i],fdascriptKnotsX[i+1]), type="l", lwd=2)
}

plot( c(  0,  0), c(0,fdascriptKnotsY[1]), type="l", lwd=2,
      xlim=c(0,230), ylim=c(-40,40), xlab = "Time (centiseconds", ylab = "Y")
lines(c(230,230), c(fdascriptKnotsY[20],0), type="l", lwd=2) 
lines(c(  0,230), c(0,0), type="l", lwd = 2, lty=2)
for (i in 1:20) {
  lines(c(Aknots[i],Aknots[i+1]), c(fdascriptKnotsY[i],fdascriptKnotsY[i]), 
        type="l", lwd=2)
}
for (i in 1:19) {
  lines(c(Aknots[i+1],Aknots[i+1]), 
        c(fdascriptKnotsY[i],fdascriptKnotsY[i+1]), type="l", lwd=2)
}

#  plot step heights as a script

indX <-  3:22;
indY <- 23:42;

StepX <- theta[indX]
StepY <- theta[indY]

plot(StepX, StepY, type="b", 
     xlab="X coordinate", ylab="Y coordinate")



