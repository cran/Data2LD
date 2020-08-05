#  Analyses of head impact data.
#  The second order linear forced differential equation is:
#           D2x(t) <- -beta0 x(t) - beta1 Dx(t) + alpha u(t)

load("HeadImpact.RData")

#  Set up and plot the time and variable vectors

attach("../Data/HeadImpactData.rda")

ImpactTime <- HeadImpactData[,2]  #  time in milliseconds
ImpactData <- HeadImpactData[,3]*100  #  acceleration in millimeters/millisecond^2
ImpactRng  <- c(0,60) # Define range time for estimated model
# plot the data along with a unit pulse
plot(ImpactTime, ImpactData, type="p", xlim=c(0,60),  ylim=c(-100,150),
     xlab="Time (milliseconds)", ylab="Acceleration (mm/msec^2)")
lines(c( 0,60), c(0,0), lty=3)
lines(c(14,14), c(0,1), lty=2)
lines(c(15,15), c(0,1), lty=2)
lines(c(14,15), c(1,1), lty=2)
#  Define the List array containing the data
ImpactDataList1 <- list(argvals=ImpactTime, y=ImpactData)
ImpactDataList <- vector("list",1)
ImpactDataList[[1]] <- ImpactDataList1

# Define a unit pulse function located at times 14-15
Pulsebasis <- create.bspline.basis(ImpactRng, 3, 1, c(0,14,15,60))
Pulsefd    <- fd(matrix(c(0,1,0),3,1),Pulsebasis)

#  Construct and plot basis for output x(t), 
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
matplot(ImpactTimefine, ImpactBasisfine, type="l", 
        xlim=c(0,60), ylim=c(0,1),
        xlab="Time t", ylab="Basis functions phi(t)")
lines(c(14,14), c(0,1), lty=2) 
lines(c(15,15), c(0,1), lty=2)
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

MSE
DpMSE
D2ppMSE

#  Optimization of the criterion for a range of smoothing values rho
#  algorithm constants
dbglev   <-  1    #  debugging level
iterlim  <- 50    #  maximum number of iterations
convrg  <- c(1e-8, 1e-4)  #  convergence criterion
#  sequence of values of \rho
gammavec <- c(0:8)
rhoVec   <- exp(gammavec)/(1+exp(gammavec))
nrho     <- length(rhoVec)
#  Matrices to hold results
dfesave <- matrix(0,nrho,1)
gcvsave <- matrix(0,nrho,1)
MSEsave <- matrix(0,nrho,1)
ISEsave <- matrix(0,nrho,1)
thesave <- matrix(0,nrho,nparam)

#  Initialize coefficient list
ImpactModelList.opt <- ImpactModelList
#  Loop through rho values
for (irho in 1:nrho) {
    rhoi <- rhoVec[irho]
    print(paste(" ------------------  rhoVeci <- ", round(rhoi,4), 
                " ------------------"))
    Data2LDResult <- Data2LD.opt(ImpactDataList, ImpactBasisList, ImpactModelList.opt,  
                                 rhoi, convrg, iterlim, dbglev)
    theta.opti     <- Data2LDResult$theta
    ImpactModelList.opt  <- modelVec2List(ImpactModelList, theta.opti)
    DataListResult <- Data2LD(ImpactDataList, ImpactBasisList, 
                              ImpactModelList.opt, rhoi)
    thesave[irho,] <- theta.opti
    dfesave[irho]  <- DataListResult$df
    gcvsave[irho]  <- DataListResult$gcv
    MSEsave[irho]  <- DataListResult$MSE
    ISEsave[irho]  <- DataListResult$ISE
}

#  Display the results for all rho values
# display the optimal parameter values
print("  rho   Stiffness   Damping   Forcing")
print(round(cbind(rhoVec,thesave),4))
# display degrees of freedom and gcv values
print(" rho    df      gcv")
print(round(cbind(rhoVec, dfesave, gcvsave),4))

#  plot the evolution of the parameters over the values of rho
par(mfrow=c(1,1))
matplot(rhoVec, thesave, type="b", pch="o", 
        xlab="\rho", ylab="parameter")

#  Evaluate the fit for parameter values at highest rho value
#  and display the resultsd
irho <- 5
theta <- thesave[irho,]
ImpactModelList <- modelVec2List(ImpactModelList, theta)
rho <- rhoVec[irho]
Data2LDResult <- Data2LD(ImpactDataList, ImpactBasisList, ImpactModelList, rho)
print(paste("MSE = ", round(Data2LDResult$MSE,4)))
print(paste("df  = ", round(Data2LDResult$df, 4)))
print(paste("gcv = ", round(Data2LDResult$gcv,4)))
Impactfd     <- Data2LDResult$XfdParList[[1]]$fd
Impactfine   <- eval.fd(ImpactTimefine, Impactfd)

# plot fit to the data
delta  <- 1
impact <- 14
plot(ImpactTimefine, Impactfine, type="l", xlim=c(0,60), ylim=c(-100,150),
     xlab="Time (msec)", ylab="Acceleration (cm/msec^2)")
points(ImpactTime, ImpactData, pch="o") 
lines(c(impact,      impact),       c(0,1), lty=2)
lines(c(impact+delta,impact+delta), c(0,1), lty=2) 
lines(c(impact+delta,impact+delta), c(0,1), lty=2) 
lines(c(impact,      impact+delta), c(1,1), lty=2) 
lines(c(0,60), c(0,0), type="l", lty=3)

# compute standard error and confidence limits for forcing
stderr  <- sqrt(diag( Data2LDResult$Var.theta))
theta   <- thesave[irho,]
thetaUP <- theta + 2*stderr
thetaDN <- theta - 2*stderr
#  display parameters along with confidence limits and standard error
print("rate constant beta0, confidence limits and std. error") 
print(round(c(theta[1], thetaDN[1], thetaUP[1], stderr[1]),4))
print("rate constant beta1, confidence limits and std. error") 
print(round(c(theta[2], thetaDN[2], thetaUP[2], stderr[2]),4))
print("rate constant alpha, confidence limits and std. error") 
print(round(c(theta[3], thetaDN[3], thetaUP[3], stderr[3]),4))

