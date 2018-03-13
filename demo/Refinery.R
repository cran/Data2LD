#  set up the 194 by 3 matrix of refinery data from a .txt file

N <- 194

#  Set up the variables and plot their data

TimeData <- RefineryData[,1]
TrayData <- RefineryData[,2]
ValvData <- RefineryData[,3]

par(mfrow=c(2,1))
plot(TimeData, TrayData, type="p") 
lines(c(67,67), c(0,4.0), type="l")
plot(TimeData, ValvData, type="p")
lines(c(67,67), c(0,0.5), type="l")

#  Construct a step basis (order 1) for smoothing valve level

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
TrayCoefList <- vector("list",2)
TrayCoefList[[1]] <- make.coef(fun=conbasis, parvec=0, estimate=TRUE, coeftype="beta")
TrayCoefList[[2]] <- make.coef(fun=conbasis, parvec=0, estimate=TRUE, coeftype="alpha")

#  Run a check on the coefficient List array, 
#  which also counts the number of estimated parameters

TraycoefResult <- coefCheck(TrayCoefList)
TrayCoefList   <- TraycoefResult$coefList
TrayNtheta     <- TraycoefResult$ntheta

#  Set up the model list for the tray variable

#  Define single homogeneous term

XTerm <- make.Xterm(variable=1, derivative=0, ncoef=1, factor=-1)
XList <- vector("list", 1)
XList[[1]] <- XTerm

#  Define the single forcing term

FTerm <- make.Fterm(ncoef=2, Ufd=Valvfd)
FList <- vector("list", 1)
FList[[1]] <- FTerm

#  Define the single differential equation in the model

TrayVariable <- make.variable(XList=XList, FList=FList, name="Tray level", order=1)

#  Set up the model List array

TrayModelList <- vector("list",1)
TrayModelList[[1]] <- TrayVariable

#  Run a check on TrayModelList

TrayModelList <- modelCheck(TrayModelList, TrayCoefList)

#  Define the List array containing the tray data

TrayData <- list(argvals=RefineryData[,1], y=RefineryData[,2])
TrayDataList  <- vector("list",1)
TrayDataList[[1]] <- TrayData

# construct the basis object for tray variable
#  Order 5 spline basis with four knots at the 67
#   to allow discontinuity in the first derivative
#  and 15 knots between 67 and 193

Traynorder <- 5
Traybreaks <- c(0, rep(67,3), seq(67, 193, len=15))
Traynbasis <- 22
TrayBasis  <- create.bspline.basis(c(0,193), Traynbasis, Traynorder, Traybreaks)
#  plot the basis
par(mfrow=c(1,1))
plot(TrayBasis, xlab="Time", ylab="B-spline basis functions")

#  Set up the basis list for the tray variable

TrayBasisList    <- vector("list",1)
TrayBasisList[[1]] <- TrayBasis

#  Set smoothing constant rho to a value specifying light smoothing
rhoVec <- 0.5 
#  Evaluate the fit to the data given the initial parameter estimates (0 and 0)
#  This also initializes the four-way tensors so that they are not re-computed
#  for subsquent analyses.

Data2LDList <- Data2LD(TrayDataList, TrayBasisList, TrayModelList, TrayCoefList, rhoVec)
MSE  <- Data2LDList$MSE    # Mean squared error for fit to data
DMSE <- Data2LDList$DpMSE  #  gradient with respect to parameter values

## Optimization of the criterion

dbglev   <-  1    #  debugging level
iterlim  <- 50    #  maximum number of iterations
convrg   <- c(1e-8, 1e-4)  #  convergence criterion

gammavec <- 0:7
rhoVec   <- exp(gammavec)/(1+exp(gammavec))
nrho    <- length(rhoVec)
dfesave <- matrix(0,nrho,1)
gcvsave <- matrix(0,nrho,1)
MSEsave <- matrix(0,nrho,1)
thesave <- matrix(0,nrho,TrayNtheta)

TrayCoefList.opt <- TrayCoefList

for (irho in 1:nrho) {
  rhoi <- rhoVec[irho]
  print(paste("rho <- ",round(rhoi,5)))
  Data2LDResult <- Data2LD.Opt(TrayDataList, TrayBasisList, TrayModelList, TrayCoefList.opt, 
                               rhoi, convrg, iterlim, dbglev)
  theta.opti <- Data2LDResult$thetastore
  TrayCoefList.opti <- modelVec2List(theta.opti, TrayCoefList)
  Data2LDList <- Data2LD(TrayDataList, TrayBasisList, TrayModelList, TrayCoefList.opti, rhoi)
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

print("    rho      df         gcv")
for (irho in 1:nrho) {
  print(round(c(rhoVec[irho], dfesave[irho], gcvsave[irho]),5))
}


## Evaluate the fit for (parameter values at highest rho value

irho <- 8  #  evaluate solution for (highest rho value
theta <- thesave[irho,]
coefList <- BAwTimeData2list(theta, coefList)
rho <- rhoVec[irho]
Data2LDList <- Data2LD(TrayDataList, TrayBasisList, TrayModelList, coefList, rho)
MSE       <- Data2LDList$MSE 
df        <- Data2LDList$df
gcv       <- Data2LDList$gcv 
ISE       <- Data2LDList$ISE 
Var.theta <- Data2LDList$Var.theta
print(paste("MSE = ", round(MSE,5)))
print(paste("df  = ", round(df,1)))
print(paste("gcv = ", round(gcv,5)))

TrayfdParList <- Data2LDList$TrayfdParList
Trayfd   <- TrayfdParList[[1]]$fd
Trayfine <- eval.fd(tfine, Trayfd)
Ufine    <- eval.fd(tfine, Ufd)

# plot fit of solution to data

plot(tfine, Trayfine, type="l", lwd=2, 
     xlab="Time", ylab="Tray level", xlim=c(0,193), ylim=c(-0.5,4.5))
points(TimeData, TrayData)
legend("x(t)", "Data")

# compute the derivative Dx(t) and the right side of the equation

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
plot(rhoVec, thesave[,1], type <- "b", lwd=2, xlab="rho", ylab="rate parameter")
plot(rhoVec, thesave[,2], type <- "b", lwd=2, xlab="rho", ylab="forcing parameter")



