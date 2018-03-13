#  The model requires a constant coefficient for the reaction speed
#  and a spline function for the coefficient of the forcing function,
#  which in this model is the unit function.
#  The estimated coefficient function is a step function,
#  or a spline function of order one, with 20 steps.
#  load the fda script data
fdascriptX <- fdascript[,,1]
fdascriptY <- fdascript[,,2]
#  Define the observation times in 100ths of a second
centisec <- seq(0,2.3,len=1401)*100
fdarange <- range(centisec)
#  Define a constant basis
conbasis <- create.constant.basis(fdarange)
#  Define the order one Bspline basis for the coefficient
#  of the forcing term.
nevent    <- 20
nAorder   <- 1
nAbasis   <- nevent
nAknots   <- nAbasis + 1
Aknots    <- seq(fdarange[1],fdarange[2],len=nAknots)
Abasisobj <- create.bspline.basis(fdarange,nAbasis,nAorder,Aknots)
#  Define the two coefficient  functions in this model
coef1 <- make.coef(conbasis,  0.04,                TRUE,   "beta")
coef2 <- make.coef(Abasisobj, matrix(1,nAbasis,1), TRUE,   "alpha")
# List array containing the coefficient lists
coefList <- vector("list",2)
coefList[[1]] <- coef1
coefList[[2]] <- coef2
# Check the coefficients
coefResult <- coefCheck(coefList)
coefList   <- coefResult$coefList
ntheta     <- coefResult$ntheta
print(paste("ntheta = ",ntheta)) 
#  Set up single homogeneous term in D^2x = -beta x 
Xterm <- make.Xterm(variable=1, derivative=0, ncoef=1, factor=-1)
XList <- vector("list",1)
XList[[1]] <- Xterm
#  Set up coefficient for forcing term \alpha(t)*1
confd <- fd(1,conbasis)
Fterm <- make.Fterm(ncoef=2, Ufd=confd, factor=1)
FList <- vector("list", 1)
FList[[1]] <- Fterm
#  make variable for X coefficient of script
Xvariable <- make.variable(name="X", order=2, XList=XList, FList=FList)
#  List array for the whole system
fdaXList = vector("list",1)
fdaXList[[1]] <- Xvariable
#  check the system specification for consistency
fdaXList <- modelCheck(fdaXList, coefList)
#  Set up the data lists for X- and Y-coordinates
yListX <- list(list(argvals=centisec, y=fdascriptX[,1]))
yListY <- list(list(argvals=centisec, y=fdascriptY[,1]))
##  The basis functions for the script will be B-splines
#  32 basis functions per second, 4 per 8 herz cycle
#  2.3 seconds 2.3*32=73.6, 
#  This implies norder + no. of interior knots <- 74 - 1 + 6 <- 79 
#  basis functions.  
nXbasis  <- 79
norder   <-  6 #  order 6 for a smooth 2nd deriv.
fdabasis <- create.bspline.basis(fdarange, nXbasis, norder)
XbasisList = vector("list",1)
XbasisList[[1]] <- fdabasis
#  Single evaluation in order to set up the 4-way tensors
rhoVec <- 0.5
Data2LDResult <- Data2LD(yListX, XbasisList, fdaXList, coefList, rhoVec)
MSE        <- Data2LDResult$MSE        #  Mean squared error for fit to data
DpMSE      <- Data2LDResult$DpMSE      #  gradient with respect to parameter values
D2ppMSE    <- Data2LDResult$D2ppMSE    #  Hessian matrix
XfdParList <- Data2LDResult$XfdParList #  List of fdPar objects for variable values 
df         <- Data2LDResult$df         #  Degrees of freedom for fit
gcv        <- Data2LDResult$gcv        #  Generalized cross-validation coefficient
ISE        <- Data2LDResult$ISE        #  Size of second term, integrated squared error
Var.theta  <- Data2LDResult$Var.theta  #  Estimate sampling variance for parameters
#  set up sequence of rho values
gamvec <- 0:7
rhoVec <- exp(gamvec)/(1+exp(gamvec))
nrho   <- length(rhoVec)
#  values controlling optimization
dbglev   <- 1    #  debugging level
iterlim  <- 50    #  maximum number of iterations
convrg  <- c(1e-6, 1e-4)  #  convergence criterion
#  Set up arrays to hold results over rho
dfesaveX <- matrix(0,nrho,2)
gcvsaveX <- matrix(0,nrho,2)
MSEsaveX <- matrix(0,nrho,2)
ISEsaveX <- matrix(0,nrho,2)
thesaveX <- matrix(0,nrho,ntheta)
#  Initialize coefficient list
coefList.optX <- coefList
#  step through rho values, optimizing at each step
#  X-coordinate:
for (irho in 1:nrho) {
  rhoi <- rhoVec[irho]
  print(paste("rho <- ",round(rhoi,5)))
  Data2LD.optResult <- Data2LD.Opt(yListX, XbasisList, fdaXList, coefList.optX, 
                                   rhoi, convrg, iterlim, dbglev)
  theta.opti    <- Data2LD.optResult$thetastore
  coefList.opti <- modelVec2List(theta.opti, coefList.optX)
  Data2LDResult <- Data2LD(yListX, XbasisList, fdaXList, coefList.opti, rhoi)
  thesaveX[irho,]   <- theta.opti
  dfesaveX[irho,1]  <- Data2LDResult$df
  gcvsaveX[irho,1]  <- Data2LDResult$gcv
  MSEsaveX[irho,1]  <- Data2LDResult$MSE
  ISEsaveX[irho,1]  <- Data2LDResult$ISE
  coefList.optX     <- coefList.opti
}
#  Set up arrays to hold results over rho
dfesaveY <- matrix(0,nrho,2)
gcvsaveY <- matrix(0,nrho,2)
MSEsaveY <- matrix(0,nrho,2)
ISEsaveY <- matrix(0,nrho,2)
thesaveY <- matrix(0,nrho,ntheta)
#  Initialize coefficient list
coefList.optY <- coefList
#  step through rho values, optimizing at each step
#  Y-coordinate:
for (irho in 1:nrho) {
  rhoi <- rhoVec[irho]
  print(paste("rho <- ",round(rhoi,5)))
  Data2LD.optResult <- Data2LD.Opt(yListY, XbasisList, fdaXList, coefList.optY, 
                                   rhoi, convrg, iterlim, dbglev)
  theta.opti    <- Data2LD.optResult$thetastore
  coefList.opti <- modelVec2List(theta.opti, coefList.optY)
  Data2LDResult <- Data2LD(yListY, XbasisList, fdaXList, coefList.opti, rhoi)
  thesaveY[irho,]   <- theta.opti
  dfesaveY[irho,1]  <- Data2LDResult$df
  gcvsaveY[irho,1]  <- Data2LDResult$gcv
  MSEsaveY[irho,1]  <- Data2LDResult$MSE
  ISEsaveY[irho,1]  <- Data2LDResult$ISE
  coefList.optY     <- coefList.opti
}
# display degrees of freedom and gcv values
print("for X:  rho    df      gcv")
print(round(cbind(rhoVec, dfesaveX, gcvsaveX),4))
print("for Y:  rho    df      gcv")
print(round(cbind(rhoVec, dfesaveY, gcvsaveY),4))
# Evaluate the fit for parameter values at highest rho value
irho <- 8
thetaX    <- thesaveX[irho,]
coefListX <- modelVec2List(thetaX, coefList)
thetaY    <- thesaveY[irho,]
coefListY <- modelVec2List(thetaY, coefList)
rho       <- rhoVec[irho]
Data2LDResultX <- Data2LD(yListX, XbasisList, fdaXList, coefListX, rho)
Data2LDResultY <- Data2LD(yListY, XbasisList, fdaXList, coefListY, rho)
print(paste("MSEX = ", round(Data2LDResultX$MSE,4)))
print(paste("dfX  = ", round(Data2LDResultX$df, 4)))
print(paste("gcvX = ", round(Data2LDResultX$gcv,4)))
print(paste("MSEY = ", round(Data2LDResultY$MSE,4)))
print(paste("dfY  = ", round(Data2LDResultY$df, 4)))
print(paste("gcvY = ", round(Data2LDResultY$gcv,4)))
FDAScriptfdX      <- Data2LDResultX$XfdParList[[1]]$fd
FDAScriptfdY      <- Data2LDResultY$XfdParList[[1]]$fd
FDAScriptTimefine <- seq(fdarange[1],fdarange[2],len=201)
FDAScriptfineX    <- eval.fd(FDAScriptTimefine, FDAScriptfdX)
FDAScriptfineY    <- eval.fd(FDAScriptTimefine, FDAScriptfdY)
# plot fit to the data
plot(FDAScriptfineX, FDAScriptfineY, type="l", 
     xlab="X coordinate", ylab="Y coordinate")
points(fdascriptX[,1], fdascriptY[,1], pch="o") 
