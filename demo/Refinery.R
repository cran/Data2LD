#  Analysis  of the refinery data
#  Dx(t) = -beta x + alpha u(t) 

#  set up the 194 by 3 matrix of refinery data

#  load all objects up to but not including the first call to Data2LD

load("refinery.RData")

#  load the data to be analyzed

load("../Data/RefineryData.rda")

N <- 194
TimeData <- RefineryData[,1]
TrayData <- RefineryData[,2]
ValvData <- RefineryData[,3]

par(mfrow=c(2,1))
plot(TimeData, TrayData, type="p") 
lines(c(67,67), c(0,4.0), type="l")
plot(TimeData, ValvData, type="p")
lines(c(67,67), c(0,0.5), type="l")

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
#  plot the basis

par(mfrow=c(1,1))
plot(TrayBasis, xlab="Time", ylab="B-spline basis functions")

#  Set up the basis list for the tray variable

TrayBasisList    <- vector("list",1)
TrayBasisList[[1]] <- TrayBasis

#  Both alpha and beta coefficient functions will be constant.
#  Define the constant basis

conbasis   <- create.constant.basis(c(0,193))
confd      <- fd(0,conbasis)
betafdPar  <- fdPar(confd)
alphafdPar <- fdPar(confd)

#  Construct a step basis (order 1) and valve level

Valvbreaks <- c(0,67,193)
Valvnbasis <- 2
Valvnorder <- 1
Valvbasis  <- create.bspline.basis(c(0,193), Valvnbasis, Valvnorder, Valvbreaks)
#  smooth the valve data 
Valvfd <- smooth.basis(TimeData, ValvData, Valvbasis)$fd

par(mfrow=c(2,1))
plot(Valvbasis, xlab="Time", ylab="B-spline basis functions")
#  plot the smooth and the data
plotfit.fd(ValvData, TimeData, Valvfd, xlab="Time", ylab="Valve setting")

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

#  Run a check on TrayModelList and make the modelList object.

checkModelList <- checkModel(TrayBasisList, TrayModelList)
TrayModelList  <- checkModelList$model
nparam         <- checkModelList$nparam

print(paste("Number of parameters =",nparam))

#  Evaluate the fit to the data given the initial parameter estimates (0 and 0)
#  This also initializes the four-way tensors so that they are not re-computed
#  for subsquent analyses.

rhoVec <- 0.5

Data2LDList <- Data2LD(TrayDataList, TrayBasisList, TrayModelList, rhoVec)

MSE  <- Data2LDList$MSE    # Mean squared error for fit to data
DMSE <- Data2LDList$DpMSE  #  gradient with respect to parameter values

MSE
DMSE

## Optimization of the criterion

dbglev   <-  1    #  debugging level
iterlim  <- 50    #  maximum number of iterations
convrg   <- c(1e-8, 1e-4)  #  convergence criterion

gammavec <- 0:7
nrho <- length(gammavec)
rhoMat   <- as.matrix(exp(gammavec)/(1+exp(gammavec)),nrho,1)
dfesave  <- matrix(0,nrho,1)
gcvsave  <- matrix(0,nrho,1)
MSEsave  <- matrix(0,nrho,1)
ISEsave  <- matrix(0,nrho,1)
thesave  <- matrix(0,nrho,nparam)

TrayModelList.opt <- TrayModelList

for (irho in 1:nrho) {
  rhoi <- rhoMat[irho]
  print(paste("rho <- ",round(rhoi,5)))
  Data2LDResult <- Data2LD.opt(TrayDataList, TrayBasisList, 
                               TrayModelList,  rhoi, convrg, 
                               iterlim, dbglev)
  theta.opti <- Data2LDResult$thetastore
  TrayModelList.opti <- modelVec2List(TrayModelList, theta.opti)
  Data2LDList <- Data2LD(TrayDataList, TrayBasisList, TrayModelList, rhoi)
  MSE       <- Data2LDList$MSE 
  df        <- Data2LDList$df
  gcv       <- Data2LDList$gcv 
  ISE       <- Data2LDList$ISE 
  Var.theta <- Data2LDList$Var.theta
  thesave[irho,] <- theta.opti
  dfesave[irho]   <- df
  gcvsave[irho]   <- gcv
  MSEsave[irho]   <- MSE
  ISEsave[irho]   <- ISE
}

# print degrees of freedom and gcv values, values of
# mean squared error for data fitting, and 
# integrated squared error for fidelity to the 
# differential equation

print("    rho      df         gcv,       MSE       ISE")
for (irho in 1:nrho) {
  print(round(c(rhoMat[irho,1], dfesave[irho], gcvsave[irho], 
                MSEsave[irho], ISEsave[irho]),5))
}

#  plot the evolution of the parameters over the values of rho

par(mfrow=c(2,1))
plot(rhoMat, thesave[,1], type <- "b", lwd=2, 
     xlab="rho", ylab="rate parameter")
plot(rhoMat, thesave[,2], type <- "b", lwd=2, 
     xlab="rho", ylab="forcing parameter")

# Evaluate the fit for (parameter values at highest rho value

irho <- 8  #  evaluate solution for (highest rho value
TrayModelList <- modelVec2List(TrayModelList, thesave[irho,])
Data2LDList <- Data2LD(TrayDataList, TrayBasisList, TrayModelList, rhoMat[irho,])
MSE       <- Data2LDList$MSE 
df        <- Data2LDList$df
gcv       <- Data2LDList$gcv 
ISE       <- Data2LDList$ISE 
Var.theta <- Data2LDList$Var.theta
print(paste("MSE = ", round(MSE,5)))
print(paste("df  = ", round(df,1)))
print(paste("gcv = ", round(gcv,5)))

# plot fit of solution to data

Trayfd   <- Data2LDList$XfdParList[[1]]$fd
tfine    <- seq(0,193,len=101)
Trayfine <- eval.fd(tfine, Trayfd)
Ufine    <- eval.fd(tfine, Valvfd)

par(mfrow=c(1,1))
plot(tfine, Trayfine, type="l", lwd=2, 
     xlab="Time", ylab="Tray level", xlim=c(0,193), ylim=c(-0.5,4.5))
points(TimeData, TrayData)

# compute the derivative Dx(t) and the right side of the equation

theta <- thesave[irho,]
DTrayData <- eval.fd(TimeData, Trayfd, 1)
DTrayfit  <- -theta[1]*Trayfine + theta[2]*Ufine

# plot the left and right sides of the equation

plot(TimeData, DTrayData, type="b", lwd=2,
     xlab="Time", ylab="DTray level", xlim=c(0,193), ylim=c(-0.01,0.12))
lines(tfine, DTrayfit)

# compute twice standard error of estimates (95# confidence limits)

stderr  <- sqrt(diag(Var.theta))
thetaUP <- theta + 2*stderr
thetaDN <- theta - 2*stderr

print(round(c(theta[1], thetaDN[1], thetaUP[1], stderr[1]),4))
print(round(c(theta[2], thetaDN[2], thetaUP[2], stderr[2]),4))

