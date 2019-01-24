#  Analyses of head impact data.
#  The second order linear forced differential equation is:
#  D2x(t) <- -beta0 x(t) - beta1 Dx(t) + alpha u(t)

#  Set up and plot the time and variable vectors
ImpactTime <- HeadImpactData[,2]  #  time in milliseconds
ImpactData <- HeadImpactData[,3]  #  acceleration in millimeters/millisecond^2
ImpactRng  <- c(0,60) # Define range time for estimated model
# plot the data along with a unit pulse
plot(ImpactTime, ImpactData, type="p", xlim=c(0,60),  ylim=c(-1.0,1.5),
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

# Set up the constant basis, make three coefficients and check them
conbasis <- create.constant.basis(ImpactRng)
# Define the three constant coefficient functions
coef1 <- make.Coef(fun=conbasis, parvec=0.073, estimate=TRUE)
coef2 <- make.Coef(fun=conbasis, parvec=0.01,  estimate=TRUE)
coef3 <- make.Coef(fun=conbasis, parvec=0.25,  estimate=TRUE)
#  Set up the coefficient list
ImpactCoefList  <- vector("list",3)
ImpactCoefList[[1]] <- coef1
ImpactCoefList[[2]] <- coef2
ImpactCoefList[[3]] <- coef3
#  check coefficient list
coefCheckList <- coefCheck(ImpactCoefList)
ImpactCoefList <- coefCheckList$coefList
coefListPrint(ImpactCoefList)
ntheta        <- coefCheckList$ntheta
print(paste("ntheta = ",ntheta))

# Define the terms in the second order linear equation
# Define the two terms in the homogeneous part of the equation
Xterm0 <- make.Xterm(variable=1, derivative=0, ncoef=1, factor=-1)
Xterm1 <- make.Xterm(variable=1, derivative=1, ncoef=2, factor=-1)
# Set up the XList vector of length two
XList <- vector("list",2)
XList[[1]] <- Xterm0
XList[[2]] <- Xterm1
# Define the forcing term
Fterm      <- make.Fterm(ncoef=3, Ufd=Pulsefd, factor=1)
# Set up the forcing term list of length one
FList      <- vector("list",1) 
FList[[1]] <- Fterm

#  Define the single differential equation in the model
ImpactVariable <- make.Variable(order=2, XList=XList, FList=FList)
#  Define list of length one containing the equation definition
ImpactVariableList=vector("list",1)
ImpactVariableList[[1]] <- ImpactVariable

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

#  Check the object for internal consistency and 
#  set up the model object
ImpactModelList <- make.Model(ImpactBasisList, ImpactVariableList, 
                              ImpactCoefList)

# An evaluation of the criterion at the initial values
rhoVec <- 0.5  #  light smoothing
Data2LDResult <- Data2LD(ImpactDataList, ImpactBasisList, 
                         ImpactModelList, ImpactCoefList, rhoVec)
MSE        <- Data2LDResult$MSE        # Mean squared error for fit to data
DpMSE      <- Data2LDResult$DpMSE      #  gradient with respect to parameter values
D2ppMSE    <- Data2LDResult$D2ppMSE    #  Hessian matrix
XfdParList <- Data2LDResult$XfdParList #  List of fdPar objects for variable values 
df         <- Data2LDResult$df         #  Degrees of freedom for fit
gcv        <- Data2LDResult$gcv        #  Generalized cross-validation coefficient
ISE        <- Data2LDResult$ISE        #  Size of second term, integrated squared error
Var.theta  <- Data2LDResult$Var.theta  #  Estimate sampling variance for parameters

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
thesave <- matrix(0,nrho,ntheta)

#  Initialize coefficient list
ImpactCoefList.opt <- ImpactCoefList
#  Loop through rho values
for (irho in 1:nrho) {
    rhoi <- rhoVec[irho]
    print(paste(" ------------------  rhoVeci <- ", round(rhoi,4), 
                " ------------------"))
    Data2LDResult <- Data2LD.opt(ImpactDataList, ImpactBasisList, 
                                 ImpactModelList, ImpactCoefList, 
                                 rhoi, convrg, iterlim, dbglev)
    theta.opti     <- Data2LDResult$theta
    ImpactCoefList.opti  <- modelVec2List(theta.opti, ImpactCoefList)
    DataListResult <- Data2LD(ImpactDataList, ImpactBasisList, 
                              ImpactModelList, ImpactCoefList.opti, rhoi)
    thesave[irho,] <- theta.opti
    dfesave[irho]  <- DataListResult$df
    gcvsave[irho]  <- DataListResult$gcv
    MSEsave[irho]  <- DataListResult$MSE
    ISEsave[irho]  <- DataListResult$ISE
    ImpactCoefList.opt   <- ImpactCoefList.opti
}

#  Display the results for all rfho values
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
irho <- 9
theta <- thesave[irho,]
ImpactCoefList <- modelVec2List(theta, ImpactCoefList)
rho <- rhoVec[irho]
Data2LDResult <- Data2LD(ImpactDataList, ImpactBasisList, ImpactModelList, ImpactCoefList, rho)
print(paste("MSE = ", round(Data2LDResult$MSE,4)))
print(paste("df  = ", round(Data2LDResult$df, 4)))
print(paste("gcv = ", round(Data2LDResult$gcv,4)))
Impactfd     <- Data2LDResult$XfdParList[[1]]$fd
Impactfine   <- eval.fd(ImpactTimefine, Impactfd)

# plot fit to the data
delta  <- 1
impact <- 14
plot(ImpactTimefine, Impactfine, type="l", xlim=c(0,60), ylim=c(-1,1.5),
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

# Explore the basis system defined by Rmat
# We calculate the eigenvalues and eigenvectors of Rmat in ascending order.
# The matrix of eigenvectors multiplying the basis matrix produces a new
# orthogonal basis system, the R-basis functions, for representing the 
# variation within and betweenestimated functions.
# The first two eigenvectors are associated with near-zero eigenvalues, and
# show the kernel of the differential operator L.  
Rmat <- Data2LDResult$Rmat
Reigen <- eigen(Rmat)
D <- Reigen$values
V <- Reigen$vectors
neig <- length(D)
Xbasismat <- eval.basis(ImpactTime, ImpactBasis)
Rbasismat <- Xbasismat %*% V
Xbasismatfine <- eval.basis(ImpactTimefine, ImpactBasis)
Rbasismatfine <- Xbasismatfine %*% V
#  plot the log10 eigenvalues
plot(1:nbasis, log10(D), type="b", pch="o", 
     xlab="Eigenvalue number", ylab="log10 Eigenvalue")
#  plot the first two  basis functions
ind <- seq(neig,1,-1)
matplot(ImpactTimefine, Rbasismatfine[,ind[1:2]],  type="l", xlim=c(0,60), ylim=c(-0.5,0.5),
        xlab="Time (milliseconds)", ylab="R-function")
lines(c(0,60), c(0,0), lty=3)       
# plot the next five
matplot(ImpactTimefine, Rbasismatfine[,ind[3:7]],  type="l", xlim=c(0,60), ylim=c(-0.5,0.55),
        xlab="Time (milliseconds)", ylab="R-function")
lines(c(0,60), c(0,0), lty=3)       
# Plot the fits for successive bases
# As we add more and more and more R-basis functions we see how the fitted 
# functions are re-conListucted.  
# About 20 basis functions in this case produces a root-mean-squared-error 
# of about 2.0, which was what we used togenerate the data.
#
#  compute root mean squared errors for successive fits
nRbasis <- dim(Rmat)[1]
RMSEsave <- matrix(0,nRbasis,1)
for (i in 1:nRbasis) {
    Zmat <- as.matrix(Rbasismat[,1:i])
    Bvec <- lsfit(Zmat,ImpactData,intercept=FALSE)$coef
    yhat <- Zmat %*% Bvec
    RMSE <- sqrt(mean(mean((ImpactData - yhat)^2)))
    RMSEsave[i] <- RMSE
}
#  plot the RMSE values for fitting the data
index <- 1:nRbasis
plot(index, RMSEsave[ind[index]], type="b", pch="o", xlim=c(1,17), ylim=c(0,0.6),
        xlab="Number of R functions", ylab="Root-mean-squared-error")
# plot each basis function and the fit to the data up to that function
index=1:nRbasis
for (i in index) {
  par(mfrow=c(2,1),ask=TRUE)
    Zmat <- as.matrix(Rbasismat[,ind[i]:neig])
    Bvec <- lsfit(Zmat,ImpactData,intercept=FALSE)$coef
    yhat <- Zmat %*% Bvec
    plot(ImpactTimefine, Rbasismatfine[,ind[i]], type="l", 
         xlab="Time (milliseconds)", ylab="R-basis function", 
         main=paste("R-basis number", i, ",  Complexity ",round(D[ind[i]],5)))
    plot(ImpactTime, ImpactData, type="p", pch="o",     
         xlab="Time (milliseconds)", ylab="Acceleration (cm/msec^2)",
         main=paste("RMSE <- ",round(RMSEsave[i],3)))
         lines(ImpactTime, yhat, col=2)
}
