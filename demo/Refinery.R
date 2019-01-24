#  set up the 194 by 3 matrix of refinery data

N <- 194
TimeData <- RefineryData[,1]
TrayData <- RefineryData[,2]
ValvData <- RefineryData[,3]

par(mfrow=c(2,1))
plot(TimeData, TrayData, type="p") 
lines(c(67,67), c(0,4.0), type="l")
plot(TimeData, ValvData, type="p")
lines(c(67,67), c(0,0.5), type="l")

#  Construct a step basis (order 1) and valve level

Valvbreaks <- c(0,67,193)
Valvnbasis <- 2
Valvnorder <- 1
Valvbasis  <- create.bspline.basis(c(0,193), Valvnbasis, Valvnorder, Valvbreaks)

par(mfrow=c(2,1))
plot(Valvbasis, xlab="Time", ylab="B-spline basis functions")
#  smooth the valve data
Valvfd <- smooth.basis(TimeData, ValvData, Valvbasis)$fd
#  plot the smooth and the data
plotfit.fd(ValvData, TimeData, Valvfd, xlab="Time", ylab="Valve setting")

#  Define the coefficient List array 
#  Both alpha and beta coefficient functions will be constant.
#  Define the constant basis
conbasis <- create.constant.basis(c(0,193))
#  Make the two coefficient function lists and
#  store these in a coefficient function list
#  Use positive random initial values for testing purposes
TrayCoefList <- vector("list",2)
TrayCoefList[[1]] <- make.Coef(fun=conbasis, parvec=exp(rnorm(1)), 
                               estimate=TRUE)
TrayCoefList[[2]] <- make.Coef(fun=conbasis, parvec=exp(rnorm(1)), 
                               estimate=TRUE)

#  Run a check on the coefficient List array, 
#  which also counts the number of estimated parameters

TraycoefResult <- coefCheck(TrayCoefList)
TrayCoefList   <- TraycoefResult$coefList
TrayNtheta     <- TraycoefResult$ntheta

#  Set up the model list for the tray variable

#  Define single homogeneous term

XTerm <- make.Xterm(variable=1, derivative=0, ncoef=1, 
                    factor=-1, name="reaction")
XList <- vector("list", 1)
XList[[1]] <- XTerm

#  Define the single forcing term

FTerm <- make.Fterm(ncoef=2, Ufd=Valvfd, name="Valve")
FList <- vector("list", 1)
FList[[1]] <- FTerm

#  Define the single differential equation in the model

TrayVariable <- make.Variable(XList=XList, FList=FList, 
                              name="Tray level", order=1)

#  Set up the model List array

TrayVariableList <- vector("list",1)
TrayVariableList[[1]] <- TrayVariable

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

#  Run a check on TrayVariableList and make the modelList object

TrayModelList <- make.Model(TrayBasisList, TrayVariableList, TrayCoefList)

#  Define the List array containing the tray data

TrayList <- list(argvals=RefineryData[,1], y=RefineryData[,2])
TrayDataList  <- vector("list",1)
TrayDataList[[1]] <- TrayList

#  Evaluate the fit to the data given the initial parameter estimates (0 and 0)
#  This also initializes the four-way tensors so that they are not re-computed
#  for subsquent analyses.

rhoVec <- 0.5
Data2LDList <- Data2LD(TrayDataList, TrayBasisList, TrayModelList, 
                       TrayCoefList, rhoVec)
MSE  <- Data2LDList$MSE    # Mean squared error for fit to data
DMSE <- Data2LDList$DpMSE  #  gradient with respect to parameter values

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
thesave  <- matrix(0,nrho,TrayNtheta)

TrayCoefList.opt <- TrayCoefList

for (irho in 1:nrho) {
  rhoi <- rhoMat[irho]
  print(paste("rho <- ",round(rhoi,5)))
  Data2LDResult <- Data2LD.opt(TrayDataList, TrayBasisList, 
                               TrayModelList,  TrayCoefList.opt, 
                               rhoi, convrg, 
                               iterlim, dbglev)
  theta.opti <- Data2LDResult$thetastore
  TrayCoefList.opti <- modelVec2List(theta.opti, TrayCoefList)
  Data2LDList <- Data2LD(TrayDataList, TrayBasisList, TrayModelList, 
                         TrayCoefList.opti, rhoi)
  MSE       <- Data2LDList$MSE 
  df        <- Data2LDList$df
  gcv       <- Data2LDList$gcv 
  ISE       <- Data2LDList$ISE 
  Var.theta <- Data2LDList$Var.theta
  thesave[irho,] <- theta.opti
  dfesave[irho]   <- df
  gcvsave[irho]   <- gcv
  MSEsave[irho]   <- MSE
}

# print degrees of freedom and gcv values

print("    rho      df         gcv,       MSE")
for (irho in 1:nrho) {
  print(round(c(rhoMat[irho,1], dfesave[irho], gcvsave[irho], MSEsave[irho]),5))
}

## Evaluate the fit for (parameter values at highest rho value

irho <- 8  #  evaluate solution for (highest rho value
TrayCoefList <- modelVec2List(thesave[irho,], TrayCoefList)
Data2LDList <- Data2LD(TrayDataList, TrayBasisList, TrayModelList, 
                       TrayCoefList, rhoMat[irho,])
MSE       <- Data2LDList$MSE 
df        <- Data2LDList$df
gcv       <- Data2LDList$gcv 
ISE       <- Data2LDList$ISE 
Var.theta <- Data2LDList$Var.theta
print(paste("MSE = ", round(MSE,5)))
print(paste("df  = ", round(df,1)))
print(paste("gcv = ", round(gcv,5)))

Trayfd   <- Data2LDList$XfdParList[[1]]$fd
tfine    <- seq(0,193,len=101)
Trayfine <- eval.fd(tfine, Trayfd)
Ufine    <- eval.fd(tfine, Valvfd)

# plot fit of solution to data

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

##  plot the evolution of the parameters over the values of rho

par(mfrow=c(2,1))
plot(rhoMat, thesave[,1], type <- "b", lwd=2, 
     xlab="rho", ylab="rate parameter")
plot(rhoMat, thesave[,2], type <- "b", lwd=2, 
     xlab="rho", ylab="forcing parameter")
