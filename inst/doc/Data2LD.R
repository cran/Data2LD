## -----------------------------------------------------------------------------
library("Data2LD")
N <- 194
TimeData <- RefineryData[,1]
TrayData <- RefineryData[,2]
ValvData <- RefineryData[,3]

## -----------------------------------------------------------------------------
plot(TimeData, TrayData, type="p", xlab="Time (sec)", ylab="Speed (km/hr)") 
lines(c(67,67), c(0,4.0), type="l")
plot(TimeData, ValvData, type="p", xlab="Time (sec)", ylab="Control")
lines(c(67,67), c(0,0.5), type="l")

## -----------------------------------------------------------------------------
  TrayList <- list(argvals=RefineryData[,1], y=RefineryData[,2])
  TrayDataList  <- vector("list",1)
  TrayDataList[[1]] <- TrayList

## -----------------------------------------------------------------------------
Valvbreaks <- c(0,67,193)
Valvnbasis <- 2
Valvnorder <- 1
Valvbasis  <- create.bspline.basis(c(0,193), Valvnbasis, Valvnorder, Valvbreaks)
Valvfd     <- smooth.basis(TimeData, ValvData, Valvbasis)$fd

## -----------------------------------------------------------------------------
  Traynorder <- 5
  Traybreaks <- c(0, rep(67,3), seq(67, 193, len=15))
  Traynbasis <- 22
  TrayBasis  <- create.bspline.basis(c(0,193), Traynbasis, Traynorder, 
                                   Traybreaks)

## -----------------------------------------------------------------------------
  par(mfrow=c(1,1))
  plot(TrayBasis, xlab="", ylab="B-spline basis functions")
  TrayBasisList    <- vector("list",1)
  TrayBasisList[[1]] <- TrayBasis

## -----------------------------------------------------------------------------
conbasis   <- create.constant.basis(c(0,193))
confd      <- fd(0,conbasis)
betafdPar  <- fdPar(confd)
alphafdPar <- fdPar(confd)

## -----------------------------------------------------------------------------
  #  Xterm Fields:    funobj     parvec  estimate  variable deriv. factor
  XTerm <- make.Xterm(betafdPar, 0.04,   TRUE,     1,       0,     -1)
  XList <- vector("list", 1)
  XList[[1]] <- XTerm

## -----------------------------------------------------------------------------
  #  Fterm Fields:    funobj      parvec  estimate  Ufd      factor
  FTerm <- make.Fterm(alphafdPar, 1.0,    TRUE,     Valvfd,  1)
  FList <- vector("list", 1)
  FList[[1]] <- FTerm

## -----------------------------------------------------------------------------
  TrayVariable <- make.Variable(XList=XList, FList=FList, 
                                name="Tray level", order=1)
  TrayModelList <- vector("list",1)
  TrayModelList[[1]] <- TrayVariable

## -----------------------------------------------------------------------------
checkModelList <- checkModel(TrayBasisList, TrayModelList, summarywrd=TRUE)
TrayModelList  <- checkModelList$model
nparam         <- checkModelList$nparam
print(paste("Number of parameters =",nparam))

## -----------------------------------------------------------------------------
rhoVec <- 0.5
Data2LDList <- Data2LD(TrayDataList, TrayBasisList, TrayModelList, rhoVec)
MSE  <- Data2LDList$MSE    # Mean squared error for fit to data
DMSE <- Data2LDList$DpMSE  #  gradient with respect to parameter values

## -----------------------------------------------------------------------------
dbglev   <-  1    #  debugging level
iterlim  <- 50    #  maximum number of iterations
convrg   <- 1e-6  #  convergence criterion

gammavec <- 0:7
nrho     <- length(gammavec)
rhoMat   <- as.matrix(exp(gammavec)/(1+exp(gammavec)),nrho,1)
dfesave  <- matrix(0,nrho,1)
gcvsave  <- matrix(0,nrho,1)
MSEsave  <- matrix(0,nrho,1)
ISEsave  <- matrix(0,nrho,1)
thesave  <- matrix(0,nrho,nparam)

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
print("    rho      df         gcv,       MSE       ISE")
for (irho in 1:nrho) {
  print(round(c(rhoMat[irho,1], dfesave[irho], gcvsave[irho], 
                MSEsave[irho], ISEsave[irho]),5))
}

## -----------------------------------------------------------------------------
irho <- nrho  #  evaluate solution for (highest rho value
irho <- 8  #  evaluate solution for (highest rho value
TrayModelList <- modelVec2List(TrayModelList, thesave[irho,])
Data2LDList <- Data2LD(TrayDataList, TrayBasisList, TrayModelList, rhoMat[irho,])
df        <- Data2LDList$df
gcv       <- Data2LDList$gcv 
MSE       <- Data2LDList$MSE 
ISE       <- Data2LDList$ISE 
Var.theta <- Data2LDList$Var.theta
beta  <- thesave[irho,1]
alpha <- thesave[irho,2]
print(paste("beta  =",round(beta,3)))
print(paste("alpha =",round(alpha,3)))
print(paste("df  = ", round(df,1)))
print(paste("gcv = ", round(gcv,5)))
print(paste("MSE = ", round(MSE,5)))
print(paste("ISE = ", round(ISE,5)))

## -----------------------------------------------------------------------------
# plot fit of solution to data
Trayfd   <- Data2LDList$XfdParList[[1]]$fd
tfine    <- seq(0,193,len=101)
Trayfine <- eval.fd(tfine, Trayfd)
Ufine    <- eval.fd(tfine, Valvfd)
par(mfrow=c(1,1))
plot(tfine, Trayfine, type="l", lwd=2, 
     xlab="Time", ylab="Tray level", xlim=c(0,193), ylim=c(-0.5,4.5))
points(TimeData, TrayData)

## -----------------------------------------------------------------------------
theta   <- thesave[nrho,]
stderr  <- sqrt(diag(Var.theta))
thetaUP <- theta + 2*stderr
thetaDN <- theta - 2*stderr
print(round(c(theta[1], thetaDN[1], thetaUP[1], stderr[1]),4))
print(round(c(theta[2], thetaDN[2], thetaUP[2], stderr[2]),4))

## -----------------------------------------------------------------------------
# compute the derivative Dx(t) and the right side of the equation
theta <- thesave[irho,]
DTrayData <- eval.fd(TimeData, Trayfd, 1)
DTrayfit  <- -theta[1]*Trayfine + theta[2]*Ufine
# plot the left and right sides of the equation
plot(TimeData, DTrayData, type="b", lwd=2,
     xlab="Time", ylab="DTray level", xlim=c(0,193), ylim=c(-0.01,0.12))
lines(tfine, DTrayfit)

## -----------------------------------------------------------------------------
T     <- 80  #  seconds
rng   <- c(0,T)
n     <- 41
tobs  <- seq(0,T,len=n)
tfine <- seq(0,T,len=501)

## -----------------------------------------------------------------------------
steporder  <- 1  #  step function basis
stepnbasis <- 4  #  four basis functions
stepbreaks <- c(0,20,40,60,80)
stepbasis  <- create.bspline.basis(rng, stepnbasis, steporder, stepbreaks)
stepcoef   <- c(60,40,80,60)  # target speed for each step
SetPtfd    <- fd(stepcoef, stepbasis)  # define the set point function

## -----------------------------------------------------------------------------
cruise0 <- function(t, y, parms) {
  DSvec    <- matrix(0,2,1)
  Uvec     <- eval.fd(t, parms$SetPtfd)
  DSvec[1] <- -y[1] + y[2]/4
  DSvec[2] <-  Uvec - y[1]
  return(list(DSvec=DSvec))
}

## -----------------------------------------------------------------------------
y0 <- matrix(0,2,1)
parms = list(SetPtfd=SetPtfd)                   #  define the input to the cruise0 function above
ytrue = lsoda(y0, tfine[1:500], cruise0, parms) # solve the equation using an initial value solver
ytrue = rbind(ytrue,matrix(c(80,60,240),1,3))   # combine with the final values

## -----------------------------------------------------------------------------
speedObs    <- approx(tfine, ytrue[,2], seq(0,80,len=41))$y
controlObs  <- approx(tfine, ytrue[,3], seq(0,80,len=41))$y

## -----------------------------------------------------------------------------
sigerr <- 2
yobs     <- matrix(0,length(tobs),2)
yobs[,1] <- as.matrix(  speedObs + rnorm(41)*sigerr)
yobs[,2] <- as.matrix(controlObs + rnorm(41)*sigerr*4)

## -----------------------------------------------------------------------------
plot(tfine, ytrue[,2], type="l", xlab="Time (secs)", ylab="Speed")
lines(c(0,T), c(60,60), lty=3)
points(tobs, yobs[,1], pch="o")
plot(tfine, ytrue[,3], type="l", xlab="Time (secs)", ylab="Control level")
lines(c(0,T), c(240,240), lty=3)
points(tobs, yobs[,2], pch="o")

## -----------------------------------------------------------------------------
Sdata <- list(argvals=tobs, y=yobs[,1])
Cdata <- list(argvals=tobs, y=yobs[,2])
cruiseDataList <- vector("list",2)
cruiseDataList[[1]] <- Sdata
cruiseDataList[[2]] <- Cdata

## -----------------------------------------------------------------------------
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

## -----------------------------------------------------------------------------
conbasis    <- create.constant.basis(rng)
confdPar    <- fdPar(conbasis)

## -----------------------------------------------------------------------------
# List object for the speed equation: A term for speed and a term for control, but no forcing
SList.XList = vector("list",2)
#  Fields:                     funobj    parvec  estimate  variable deriv. factor
SList.XList[[1]] <- make.Xterm(confdPar, 1,      TRUE,     1,       0,     -1)
SList.XList[[2]] <- make.Xterm(confdPar, 1/4,    TRUE,     2,       0,      1)
SList.FList = NULL
SList = make.Variable("Speed", 1, SList.XList, SList.FList)
# List object for the control equation:  a term for speed, and a zero-multiplied term for control plus a term for the forcing function SetPtfd
CList.XList <- vector("list",2)
#  Fields:                     funobj    parvec  estimate  variable deriv. factor
CList.XList[[1]] <- make.Xterm(confdPar, 1,      TRUE,     1,       0,     -1)
CList.XList[[2]] <- make.Xterm(confdPar, 0,      FALSE,    2,       0,      1)
CList.FList <- vector("list",1)
#  Fields:                     funobj    parvec  estimate  Ufd      factor
CList.FList[[1]] <- make.Fterm(confdPar, 1,      TRUE,     SetPtfd, 1)
CList <- make.Variable("Control", 1, CList.XList, CList.FList)
#  List array for the whole system
cruiseVariableList <- vector("list",2)
cruiseVariableList[[1]] <- SList
cruiseVariableList[[2]] <- CList

## -----------------------------------------------------------------------------
cruiseModelList <- vector("list",2)
cruiseModelList[[1]] <- SList
cruiseModelList[[2]] <- CList

## -----------------------------------------------------------------------------
cruiseModelCheckList <- checkModel(cruiseBasisList, cruiseModelList)
cruiseModelList <- cruiseModelCheckList$modelList
nparam <- cruiseModelCheckList$nparam
print(paste("total number of parameters = ", nparam))

## -----------------------------------------------------------------------------
rhoVec <- 0.5*matrix(1,1,2)
Data2LDList <- Data2LD(cruiseDataList, cruiseBasisList, cruiseModelList, rhoVec)
print(Data2LDList$MSE)
print(Data2LDList$DpMSE)

## -----------------------------------------------------------------------------
conv    <- 1e-4  
iterlim <- 20    
dbglev  <-  1   

Gvec    <- c(0:7)
nrho    <- length(Gvec)
rhoMat  <- matrix(exp(Gvec)/(1+exp(Gvec)),nrho,2)
dfesave <- matrix(0,nrho,1)
gcvsave <- matrix(0,nrho,1)
MSEsave <- matrix(0,nrho,2)
thesave <- matrix(0,nrho,nparam)

## -----------------------------------------------------------------------------
cruiseModelList.opt <- cruiseModelList  #  initialize the optimizing coefficient values
for (irho in 1:nrho) {   
  rhoVeci <- rhoMat[irho,]
  print(paste(" ------------------  rhoVeci <- ", round(rhoVeci[1],4), 
              " ------------------"))
  Data2LDOptList <- Data2LD.opt(cruiseDataList, cruiseBasisList, cruiseModelList.opt, 
                                rhoVeci, conv, iterlim, dbglev)
  theta.opti <- Data2LDOptList$theta  # optimal parameter values
  cruiseModelList.opt <- modelVec2List(cruiseModelList, theta.opti) 
  #  evaluate fit at optimal values and store the results
  Data2LDList    <- Data2LD(cruiseDataList, cruiseBasisList, cruiseModelList.opt, rhoVeci,
                            TRUE)
  thesave[irho,]  <- theta.opti
  dfesave[irho]   <- Data2LDList$df
  gcvsave[irho]   <- Data2LDList$gcv
  x1fd            <- Data2LDList$XfdParList[[1]]$fd
  x1vec           <- eval.fd(tobs,  x1fd)
  msex1           <- mean((x1vec - speedObs)^2)
  x2fd            <- Data2LDList$XfdParList[[2]]$fd
  x2vec           <- eval.fd(tobs,  x2fd)
  msex2           <- mean((x2vec - controlObs)^2)
  MSEsave[irho,1] <- msex1
  MSEsave[irho,2] <- msex2
}

## -----------------------------------------------------------------------------
ind <- 1:nrho
#  print df, gcv and MSEs
print("    rho      df        gcv        RMSE:")
print(cbind(round(rhoMat,4),  round(dfesave,1), 
            round(gcvsave,1), round(sqrt(MSEsave),2)))

## -----------------------------------------------------------------------------
matplot(rhoMat, thesave, type="b", xlab="rho", ylab="theta(rho)")

## -----------------------------------------------------------------------------
rho.opt   <- rhoMat[nrho,]
theta.opt <- thesave[nrho,]
#  convert the optimal parameter values to optimal cruiseModelList
cruiseModelList.opt <- modelVec2List(cruiseModelList, theta.opt)
#  evaluate the solution at the optimal solution
DataLDList <- Data2LD(cruiseDataList, cruiseBasisList, cruiseModelList.opt, rho.opt,
                      TRUE)

## -----------------------------------------------------------------------------
XfdParList <- Data2LDList$XfdParList
Xfd1 <- XfdParList[[1]]$fd
Xfd2 <- XfdParList[[2]]$fd
Xvec1 <- eval.fd(tfine, Xfd1)
Xvec2 <- eval.fd(tfine, Xfd2)
Uvec  <- eval.fd(tfine, SetPtfd)
cruiseDataList1 <- cruiseDataList[[1]]
plot(tfine, Xvec1, type="l", xlim=c(0,80), ylim=c(0,100), 
     ylab="Speed",
     main=paste("RMSE =",round(sqrt(MSEsave[3,1]),4)))
lines(tfine, ytrue[,2], lty=2)
points(cruiseDataList1$argvals, cruiseDataList1$y, pch="o")
#lines(tfine, Uvec, lty=2)
cruiseDataList2 <- cruiseDataList[[2]]
plot(tfine, Xvec2, type="l", 
     xlim=c(0,80), ylim=c(0,400),
     xlab="Time (sec)", ylab="Control",
     main=paste("RMSE =",round(sqrt(MSEsave[3,2]),4)))
lines(tfine, ytrue[,3], lty=2)
points(cruiseDataList2$argvals, cruiseDataList2$y, pch="o")

## -----------------------------------------------------------------------------
stddev.opt <- sqrt(diag(DataLDList$Var.theta))
theta.tru <- c(1, 1/4,  1, 0,  1)
print("    True      Est.      Std. Err. Low CI    Upr CI:")
for (i in 1:nparam) {
  print(round(c(theta.tru[i], 
                theta.opt[i], 
                stddev.opt[i], 
                theta.opt[i]-2*stddev.opt[i], 
                theta.opt[i]+2*stddev.opt[i]), 4))
}

## -----------------------------------------------------------------------------
#  parMap = matrix(c(0, 0, 1, 0, -1),5,1)  
#  Data2LDOptList <- Data2LD.opt(cruiseDataList, cruiseBasisList, cruiseModelList.opt, 
#                                rhoVeci, conv, iterlim, dbglev, parMap)

