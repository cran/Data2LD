#  In this example, a coefficient is nonconstant and also
#  a nonlinear function of its parameters.  The reaction speed
#  with which average temperature in Montreal responds to
#  solar forcing is required to be positive.
#  set up five-day block averages with winter centering
#  Here we use 73 five-day block averages as the data with the block
#  centers as the time points in order to speed up computation
daytime73 <- seq(2.5,362.5,len=73)
dayrange  <- c(0,365)
#  set up block averages for Canadian temperature data,
#  available from the fda package
tempav <- CanadianWeather$dailyAv[,,"Temperature.C"]
tempav73 <- matrix(0,73,35)
m2 <- 0
for (  i in 1 : 73 ) {
  m1 <- m2 + 1
  m2 <- m2 + 5
  tempavi <- apply(tempav[m1:m2,1:35],2,mean)
  tempav73[i,] <- tempavi
} 
#  center the time of the year on winter
winterind73  <- c( 37:73,1: 36)
tempav73 <- tempav73[winterind73,]
#  select the data for Montreal
station <- 12 
#  Define the two forcing functions:
#  constant function Ufd for constant forcing
Uconbasis <- create.constant.basis(dayrange)
Uconfd    <- fd(1, Uconbasis)
Uconvec   <- matrix(1,73,1)
#  cosine function for solar radiation forcing
uvec      <- -cos((2*pi/365)*(daytime73+10+182))
Ucosbasis <- create.fourier.basis(dayrange, 3)
Ucosfd    <- smooth.basis(daytime73, uvec, Ucosbasis)$fd
#  A fourier basis for defining the positive homogeneous 
#  coefficient
nWbasis   <- 7
Wbasisobj <- create.fourier.basis(dayrange, nWbasis)
#  constant forcing coefficient for constant forcing
nAbasisC   <- 1
AbasisobjC <- create.constant.basis(dayrange)
#  fourier coefficint function for radiative forcing
nAbasisF   <- 1
AbasisobjF <- create.constant.basis(dayrange)
#  Set up the list object for the positive coefficient for 
#  the homogeneous term.  
linfun <- list(fd=fun.explinear, Dfd=fun.Dexplinear, more=Wbasisobj)
#  Define the coefficient
coefListW  <-list(fun=linfun, parvec=matrix(0,nWbasis,1), 
                  estimate=TRUE, coeftype="beta")
#  Coefficient for constant forcing
coefListA1 <- make.coef(fun=AbasisobjC, parvec=1,
                        estimate=TRUE, coeftype="alpha")
#  Coeffcient for cosine forcing
coefListA2 <- make.coef(fun=AbasisobjF, parvec=1, 
                        estimate=TRUE, coeftype="alpha")
#  coefList constructed
coefList <- vector("list",3)
coefList[[1]] <- coefListW
coefList[[2]] <- coefListA1
coefList[[3]] <- coefListA2
#  check the coefficient list
coefResult <- coefCheck(coefList)
coefList <- coefResult$coefList
ntheta   <- coefResult$ntheta
print(paste("ntheta = ",ntheta))
#  define homogeneous term and list container
Xterm <- make.Xterm(variable=1, derivative=0, ncoef=1, factor=-1)
XList <- vector("list", 1)
XList[[1]] <- Xterm
#  define two forcing terms
Fterm1 <- make.Fterm(Ufd=Uconfd, ncoef=2, factor=1)
Fterm2 <- make.Fterm(Ufd=Ucosfd, ncoef=3, factor=1)
#  forcing term container
FList <- vector("list", 2)
FList[[1]] <- Fterm1
FList[[2]] <- Fterm2
#  model list object
TempVar <- make.variable(name="temperature", order=1, XList=XList, FList=FList)
# set up AvTempList
AvTempList <- vector("list",1)
AvTempList[[1]] <- TempVar
AvTempList <- modelCheck(AvTempList, coefList)
#  basis for 5-day block averages
norder    <- 5
daybreaks <- seq(0,365,5)
nbreaks   <- length(daybreaks)
nbasis    <- norder + nbreaks - 2
daybasis  <- create.bspline.basis(dayrange, nbasis)
XbasisList <- vector("list",1)
XbasisList[[1]] <- daybasis
#  set up yList
yList1 <- list(argvals=as.matrix(daytime73), y=as.matrix(tempav73[,station]))
yList  <- vector("list",1)
yList[[1]]    <- yList1
# An evaluation of the criterion at the initial values
rhoVec <- 0.5
#  This command causes Data2LD to set up and save the tensors
Data2LDResult <- Data2LD(yList, XbasisList, AvTempList, coefList, rhoVec)
MSE        <- Data2LDResult$MSE        # Mean squared error for fit to data
DpMSE      <- Data2LDResult$DpMSE      #  gradient with respect to parameter values
D2ppMSE    <- Data2LDResult$D2ppMSE    #  Hessian matrix
XfdParList <- Data2LDResult$XfdParList #  List of fdPar objects for variable values 
df         <- Data2LDResult$df         #  Degrees of freedom for fit
gcv        <- Data2LDResult$gcv        #  Generalized cross-validation coefficient
ISE        <- Data2LDResult$ISE        #  Size of second term, integrated squared error
Var.theta  <- Data2LDResult$Var.theta  #  Estimate sampling variance for parameters
#  set constants for estimation algorithm Data2LD.Opt
dbglev  <-  1    
iterlim <- 50    
convrg  <- c(1e-8, 1e-4)  
#  define rhovec using the logit function
gammavec <- seq(0,7,1)
rhoVec   <- exp(gammavec)/(1+exp(gammavec))
nrho     <- length(rhoVec)
#  Matrices to hold results
dfesave <- matrix(0,nrho,1)
gcvsave <- matrix(0,nrho,1)
MSEsave <- matrix(0,nrho,1)
ISEsave <- matrix(0,nrho,1)
thesave <- matrix(0,nrho,ntheta)
#  Initialize coefficient list
coefList.opt <- coefList
#  Loop through rho values
for (  irho  in  1 : nrho ) {
  rhoVeci <- rhoVec[irho]
  print(paste('Rho = ',round(rhoVeci,5)))
  OptList <- Data2LD.Opt(yList, XbasisList, AvTempList, coefList.opt, 
                         rhoVeci, convrg, iterlim, dbglev)
  theta.opti <- OptList$thetastore
  coefList.opti <- modelVec2List(theta.opti, coefList)    
  Data2LDResult <- Data2LD(yList, XbasisList, AvTempList, coefList.opt, rhoVeci)
  thesave[irho,] <- theta.opti
  dfesave[irho]  <- Data2LDResult$df
  gcvsave[irho]  <- Data2LDResult$gcv
  MSEsave[irho]  <- Data2LDResult$MSE
  ISEsave[irho]  <- Data2LDResult$ISE
  coefList.opt   <- coefList.opti
} 
# display degrees of freedom and gcv values
print('    rho      df         gcv')
print(cbind(rhoVec, dfesave, gcvsave))
# Evaluate the fit for parameter values at highest rho value
irho <- nrho
theta <- thesave[irho,]
coefList <- modelVec2List(theta, coefList)
rhoi <- rhoVec[irho]
Data2LDList <- Data2LD(yList, XbasisList, modelList, coefList.opti, rhoi)
MSE       <- Data2LDList$MSE 
df        <- Data2LDList$df
gcv       <- Data2LDList$gcv 
ISE       <- Data2LDList$ISE 
Var.theta <- Data2LDList$Var.theta
print(paste("MSE = ", round(Data2LDResult$MSE,4)))
print(paste("df  = ", round(Data2LDResult$df, 4)))
print(paste("gcv = ", round(Data2LDResult$gcv,4)))
#  Set up functional data object and fine mesh for plotting variable
XfdParList <- Data2LDList$XfdParList
tempfdPar <- XfdParList[[1]]
tempfd    <- tempfdPar$fd
tfine     <- seq(0,365,len=101)
tempfine  <- eval.fd(tfine, tempfd)
#  Plot the fit to the data
plot(tfine, tempfine, type="l", lwd=2, xlim=c(0,365), ylim=c(-15,25),
     xlab="Time (days)", ylab="Temperature (deg C)")
points(daytime73, tempav73[,station], pch="o") 
lines(daytime73, uvec, lty=4)
# plot the optimal parameter values as functions of rho
indW <- 1:nWbasis
indA1 <- (nWbasis+1):(nWbasis + nAbasisC)
indA2 <- (nWbasis+nAbasisC+1):ntheta
#  plot flow of beta(t) parameters
matplot(rhoVec, thesave[,indW], type="b", lwd=2, 
        xlab="rho", ylab="log beta coefficients")
#  plot flow of alpha 1 and alpha 2
par(mfrow=c(2,1))
plot(rhoVec, thesave[,indA1], type="b", ylab="alpha.1")
plot(rhoVec, thesave[,indA2], type="b", xlab="rho", ylab="alpha.2")
#  plot \beta with confidence intervals
betafd     <- fd(thesave[irho,indW], Wbasisobj)
betavec    <- exp(eval.fd(daytime73,betafd))
Var.thetaW <- Var.theta[indW,indW]
basismatW  <- (betavec %*% matrix(1,1,nWbasis)) * eval.basis(daytime73, Wbasisobj)
CI.beta    <- 2*sqrt(diag(basismatW %*% Var.thetaW %*% t(basismatW)))
betamat <- cbind(betavec,betavec  + CI.beta,betavec  - CI.beta)
par(mfrow=c(1,1))
matplot(daytime73, betamat , type="l", lty=c(1,4,4), col=1, lwd=2, 
        xlab="Time (days)", ylab="beta", xlim=c(0,365), ylim=c(0,0.05))
#  display 95# confidence limits for \alpha1 and \alpha2
alpha1 <- thesave[irho,indA1]
alpha2 <- thesave[irho,indA2]
alpha1.stderr <- sqrt(Var.theta[indA1])
alpha2.stderr <- sqrt(Var.theta[indA2])
print('    alpha.1   lower CI upper CI')
print(c(alpha1, alpha1-2*alpha1.stderr, alpha1+2*alpha1.stderr))
print('    alpha.2   lower CI upper CI')
print(c(alpha2, alpha2-2*alpha2.stderr, alpha2+2*alpha2.stderr))



