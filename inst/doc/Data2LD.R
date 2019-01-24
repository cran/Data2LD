## ----} setup, include=FALSE----------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  fig.width  = 7,
  fig.height = 5
  )

## ------------------------------------------------------------------------
library("Data2LD")
N <- 194
TimeData <- RefineryData[,1]
TrayData <- RefineryData[,2]
ValvData <- RefineryData[,3]

## ------------------------------------------------------------------------
plot(TimeData, TrayData, type="p", xlab="Time (sec)", ylab="Speed (km/hr)") 
lines(c(67,67), c(0,4.0), type="l")
plot(TimeData, ValvData, type="p", xlab="Time (sec)", ylab="Control")
lines(c(67,67), c(0,0.5), type="l")

## ------------------------------------------------------------------------
  TrayList <- list(argvals=RefineryData[,1], y=RefineryData[,2])
  TrayDataList  <- vector("list",1)
  TrayDataList[[1]] <- TrayList

## ------------------------------------------------------------------------
Valvbreaks <- c(0,67,193)
Valvnbasis <- 2
Valvnorder <- 1
Valvbasis  <- create.bspline.basis(c(0,193), Valvnbasis, Valvnorder, Valvbreaks)
Valvfd     <- smooth.basis(TimeData, ValvData, Valvbasis)$fd

## ------------------------------------------------------------------------
  Traynorder <- 5
  Traybreaks <- c(0, rep(67,3), seq(67, 193, len=15))
  Traynbasis <- 22
  TrayBasis  <- create.bspline.basis(c(0,193), Traynbasis, Traynorder, 
                                   Traybreaks)

## ------------------------------------------------------------------------
  par(mfrow=c(1,1))
  plot(TrayBasis, xlab="", ylab="B-spline basis functions")
  TrayBasisList    <- vector("list",1)
  TrayBasisList[[1]] <- TrayBasis

## ------------------------------------------------------------------------
  conbasis <- create.constant.basis(c(0,193))
  TrayCoefList <- vector("list",2)
  TrayCoefList[[1]] <- make.Coef(fun=conbasis, parvec=exp(rnorm(1)), 
                               estimate=TRUE)
  TrayCoefList[[2]] <- make.Coef(fun=conbasis, parvec=exp(rnorm(1)), 
                               estimate=TRUE)

## ------------------------------------------------------------------------
TraycoefResult <- coefCheck(TrayCoefList)
TrayCoefList   <- TraycoefResult$coefList
TrayNtheta     <- TraycoefResult$ntheta
print(paste("Total number of parameters = ",TrayNtheta))

## ------------------------------------------------------------------------
  #funList       <- list(fun=fun.explinear, Dfun=fun.Dexplinear)
  #coefList1     <- make.Coef(funList, 0, TRUE, "beta")
  #coefList2     <- make.Coef(confdPar, 0, TRUE, "alpha")
  #coefList      <- vector("list",2)
  #coefList[[1]] <- coefList1
  #coefList[[2]] <- coefList2

## ------------------------------------------------------------------------
  XTerm <- make.Xterm(variable=1, derivative=0, ncoef=1, 
                      factor=-1, name="reaction")
  XList <- vector("list", 1)
  XList[[1]] <- XTerm

## ------------------------------------------------------------------------
  FTerm <- make.Fterm(ncoef=2, Ufd=Valvfd, name="Valve")
  FList <- vector("list", 1)
  FList[[1]] <- FTerm

## ------------------------------------------------------------------------
  TrayVariable <- make.Variable(XList=XList, FList=FList, 
                                name="Tray level", order=1)

## ------------------------------------------------------------------------
  TrayVariableList <- vector("list",1)
  TrayVariableList[[1]] <- TrayVariable

## ------------------------------------------------------------------------
  TrayModelList <- make.Model(TrayBasisList, TrayVariableList, TrayCoefList)

## ------------------------------------------------------------------------
  rhoVec <- 0.5
  Data2LDList <- Data2LD(TrayDataList, TrayBasisList, TrayModelList, 
                       TrayCoefList, rhoVec)
  MSE  <- Data2LDList$MSE    #  Mean squared error for fit to data
  DMSE <- Data2LDList$DpMSE  #  Gradient with respect to parameter values

## ------------------------------------------------------------------------
dbglev   <-  1    #  debugging level
iterlim  <- 50    #  maximum number of iterations
convrg   <- c(1e-8, 1e-4)  #  convergence criterion

gammavec <- 0:6
nrho     <- length(gammavec)
rhoMat   <- as.matrix(exp(gammavec)/(1+exp(gammavec)),nrho,1)
dfesave  <- matrix(0,nrho,1)
gcvsave  <- matrix(0,nrho,1)
MSEsave  <- matrix(0,nrho,1)
thesave  <- matrix(0,nrho,TrayNtheta)

## ------------------------------------------------------------------------
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

## ------------------------------------------------------------------------
#   rho   df     gcv     MSE
# 0.500  20.5  0.0043  0.00347
# 0.731  18.9  0.0043  0.00351
# 0.881  16.3  0.0043  0.00364
# 0.953  13.0  0.0046  0.00403
# 0.982   9.6  0.0054  0.00491
# 0.993   6.7  0.0066  0.00617
# 0.998   4.4  0.0076  0.00723

## ------------------------------------------------------------------------
irho <- nrho  #  evaluate solution for (highest rho value
TrayCoefList <- modelVec2List(thesave[irho,], TrayCoefList)
Data2LDList <- Data2LD(TrayDataList, TrayBasisList, TrayModelList, 
                       TrayCoefList, rhoMat[irho,])

## ------------------------------------------------------------------------
Trayfd   <- Data2LDList$XfdParList[[1]]$fd
tfine    <- seq(0,193,len=101)
Trayfine <- eval.fd(tfine, Trayfd)
par(mfrow=c(1,1))
plot(tfine, Trayfine, type="l", lwd=2, 
     xlab="Time", ylab="Tray level", xlim=c(0,193), ylim=c(-0.5,4.5))
points(TimeData, TrayData)

## ------------------------------------------------------------------------
T     <- 80  #  seconds
rng   <- c(0,T)
n     <- 41
tobs  <- seq(0,T,len=n)
tfine <- seq(0,T,len=501)

## ------------------------------------------------------------------------
steporder  <- 1  #  step function basis
stepnbasis <- 4  #  four basis functions
stepbreaks <- c(0,20,40,60,80)
stepbasis  <- create.bspline.basis(rng, stepnbasis, steporder, stepbreaks)
stepcoef   <- c(60,40,80,60)  # target speed for each step
SetPtfd    <- fd(stepcoef, stepbasis)  # define the set point function

## ------------------------------------------------------------------------
cruise0 <- function(t, y, parms) {
  DSvec    <- matrix(0,2,1)
  Uvec     <- eval.fd(t, parms$SetPtfd)
  DSvec[1] <- -y[1] + y[2]/4
  DSvec[2] <-  Uvec - y[1]
  return(list(DSvec=DSvec))
}

## ------------------------------------------------------------------------
y0 <- matrix(0,2,1)
parms = list(SetPtfd=SetPtfd)
ytrue = lsoda(y0, tfine[1:500], cruise0, parms)
ytrue = rbind(ytrue,matrix(c(80,60,240),1,3))

## ------------------------------------------------------------------------
speedObs    <- approx(tfine, ytrue[,2], seq(0,80,len=41))$y
controlObs  <- approx(tfine, ytrue[,3], seq(0,80,len=41))$y

## ------------------------------------------------------------------------
sigerr <- 2
yobs     <- matrix(0,length(tobs),2)
yobs[,1] <- as.matrix(  speedObs + rnorm(41)*sigerr)
yobs[,2] <- as.matrix(controlObs + rnorm(41)*sigerr*4)

## ------------------------------------------------------------------------
plot(tfine, ytrue[,2], type="l", xlab="Time (secs)", ylab="Speed")
lines(c(0,T), c(60,60), lty=3)
points(tobs, yobs[,1], pch="o")
plot(tfine, ytrue[,3], type="l", xlab="Time (secs)", ylab="Control level")
lines(c(0,T), c(240,240), lty=3)
points(tobs, yobs[,2], pch="o")

## ------------------------------------------------------------------------
Sdata <- list(argvals=tobs, y=yobs[,1])
Cdata <- list(argvals=tobs, y=yobs[,2])
cruiseDataList <- vector("list",2)
cruiseDataList[[1]] <- Sdata
cruiseDataList[[2]] <- Cdata

## ------------------------------------------------------------------------
cruiseBasisList <- vector("list",2)
delta <- 2*(1:10)
breaks   <- c(0, delta, 20, 20+delta, 40, 40+delta, 60, 60+delta)
nbreaks  <- length(breaks)
nSorder <- 5
nSbasis <- nbreaks + nSorder - 2
Sbasis  <- create.bspline.basis(c(0,80), nSbasis, nSorder, breaks)
cruiseBasisList[[1]] <- Sbasis
nCorder <- 4
nCbasis <- nbreaks + nCorder - 2
Cbasis  <- create.bspline.basis(c(0,80), nCbasis, nCorder, breaks)
cruiseBasisList[[2]] <- Cbasis

## ------------------------------------------------------------------------
conbasis    <- create.constant.basis(rng)
confdPar    <- fdPar(conbasis)

## ------------------------------------------------------------------------
cruiseCoefListS.S <- make.Coef(confdPar,    1, TRUE)
cruiseCoefListS.C <- make.Coef(confdPar,  1/4, TRUE)
cruiseCoefListC.S <- make.Coef(confdPar,    1, TRUE)
cruiseCoefListC.C <- make.Coef(confdPar,    0, FALSE)

cruiseCoefList <- list(4,1)
cruiseCoefList[[1]] <- cruiseCoefListS.S
cruiseCoefList[[2]] <- cruiseCoefListS.C
cruiseCoefList[[3]] <- cruiseCoefListC.S
cruiseCoefList[[4]] <- cruiseCoefListC.C

## ------------------------------------------------------------------------
Result   <- coefCheck(cruiseCoefList)
cruiseCoefList <- Result[[1]] 
 theta         <- Result[[2]]
ntheta         <- Result[[3]]

## ------------------------------------------------------------------------
SList.XList = vector("list",2)
#  Fields:                     variable ncoef derivative factor
SList.XList[[1]] <- make.Xterm(1,       1,    0,         -1)
SList.XList[[2]] <- make.Xterm(2,       2,    0,          1)
SList.FList = NULL
SList = make.Variable("Speed", 1, SList.XList, SList.FList)
# List object for the control equation
CList.XList <- vector("list",2)
CList.XList[[1]] <- make.Xterm(1,       3,    0,         -1)
CList.XList[[2]] <- make.Xterm(2,       4,    0,          1)
CList.FList <- vector("list",1)
#  Fields:                 variable ncoef Ufd      factor
CList.FList[[1]] <- make.Fterm(3, SetPtfd, 1)
CList <- make.Variable("Control", 1, CList.XList, CList.FList)
#  List array for the whole system
cruiseVariableList <- vector("list",2)
cruiseVariableList[[1]] <- SList
cruiseVariableList[[2]] <- CList

## ------------------------------------------------------------------------
 cruiseModelList <- make.Model(cruiseBasisList, cruiseVariableList, cruiseCoefList)

## ------------------------------------------------------------------------
rhoMat = 0.5*matrix(1,1,2)
Data2LDList <- Data2LD(cruiseDataList, cruiseBasisList, 
                       cruiseModelList, cruiseCoefList, rhoMat)
MSE     <- Data2LDList$MSE
DpMSE   <- Data2LDList$DpMSE

## ------------------------------------------------------------------------
conv    <- 1e-4  
iterlim <- 20    
dbglev  <-  1   

Gvec    <- c(0:7)
nrho    <- length(Gvec)
rhoMat  <- matrix(exp(Gvec)/(1+exp(Gvec)),nrho,2)
dfesave <- matrix(0,nrho,1)
gcvsave <- matrix(0,nrho,1)
MSEsave <- matrix(0,nrho,2)
thesave <- matrix(0,nrho,ntheta)

cruiseCoefList.opt <- cruiseCoefList  
for (irho in 1:nrho) {   
  rhoVeci <- rhoMat[irho,]
  print(paste(" ------------------  rhoVeci <- ", round(rhoVeci[1],4), 
              " ------------------"))
  Data2LDOptList <- Data2LD.opt(cruiseDataList, cruiseBasisList, 
                                cruiseModelList, cruiseCoefList.opt, 
                                rhoVeci, conv, iterlim, dbglev)
  theta.opti     <- Data2LDOptList$theta  
  cruiseCoefList.opti  <- modelVec2List(theta.opti, cruiseCoefList)  
  Data2LDList    <- Data2LD(cruiseDataList, cruiseBasisList, 
                            cruiseModelList, cruiseCoefList.opti, 
                            rhoVeci)
  thesave[irho,]  <- theta.opti
  dfesave[irho]   <- Data2LDList$df
  gcvsave[irho]   <- Data2LDList$gcv
  x1fd            <- Data2LDList$XfdParList[[1]]$fd
  x2fd            <- Data2LDList$XfdParList[[2]]$fd
  x1vec           <- eval.fd(tobs,  x1fd)
  x2vec           <- eval.fd(tobs,  x2fd)
  MSEsave[irho,1] <- mean((x1vec - speedObs)^2)
  MSEsave[irho,2] <- mean((x2vec - controlObs)^2)
  cruiseCoefList.opt    <- cruiseCoefList.opti  #  update the optimizing coefficient values
}

## ------------------------------------------------------------------------
#   rho  df     gcv   RMSE
# 0.500  44.7   20.1  1.98
# 0.731  34.0   26.0  2.10
# 0.881  25.4   33.8  2.33
# 0.953  12.8   34.6  0.86
# 0.982   6.9   37.5  0.61
# 0.993   4.0   39.3  0.59
# 0.998   2.7   40.4  0.82
# 0.999   2.2   41.3  1.27

## ------------------------------------------------------------------------
rho.opt   <- rhoMat[nrho,]
theta.opt <- thesave[nrho,]
cruiseCoefList.opt <- modelVec2List(theta.opt, cruiseCoefList)
DataLDList <- Data2LD(cruiseDataList, cruiseBasisList, 
                      cruiseModelList, cruiseCoefList.opt, rho.opt)

## ------------------------------------------------------------------------
stddev.opt <- sqrt(diag(DataLDList$Var.theta))

#  True   Est.  Std. Err. Low CI  Upr CI
#  1.00  1.037    0.196    0.645   1.430
#  0.25  0.259    0.049    0.161   0.356
#  1.00  0.942    0.045    0.853   1.032

## ------------------------------------------------------------------------
XfdParList <- Data2LDList$XfdParList
Xfd1 <- XfdParList[[1]]$fd
Xfd2 <- XfdParList[[2]]$fd
Xvec1 <- eval.fd(tfine, Xfd1)
Xvec2 <- eval.fd(tfine, Xfd2)
Uvec  <- eval.fd(tfine, SetPtfd)
cruiseDataList1 <- cruiseDataList[[1]]
cruiseDataList2 <- cruiseDataList[[2]]
plot(tfine, Xvec1, type="l", xlim=c(0,80), ylim=c(0,100), 
     ylab="Speed",
     main=paste("RMSE =",round(sqrt(MSEsave[3,1]),4)))
points(cruiseDataList1$argvals, cruiseDataList1$y, pch="o")
lines(tfine, Uvec, lty=2)
plot(tfine, Xvec2, type="l", 
     xlab="Time (sec)", ylab="Control")
points(cruiseDataList2$argvals, cruiseDataList2$y, pch="o")

