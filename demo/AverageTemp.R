###                 Data2LD analyses of weather data
#  Montreal's average temperature, centered on January 1, is fit by
#  a first order linear equation with a time-varying coefficient and
#  forced by a constant and by a cosine function representing solar
#  heating translated by 10 days:
#         DT(t) <- - \beta(t) + \alpha_1 + \alpha_2 U(t)
#  where U(t) <- -cos((2*pi/365)*(t+192)).
#  Rate function \beta(t) is constrained to be positive by expressing it
#  as the exponential of a B-spline function with seven basis functions.
#  This illustrates the use of a user-defined pair of functions for 
#  \beta and its derivative with respect to the coefficients.
#  The user-defined functions are:
#
# function bval <- fun_explinear(t, bvec, Bbasisobj)
# if isa_fdPar(Bbasisobj)
#     Bbasisobj <- getbasis(Bbasisobj);
# end
# basismat <- eval_basis(t,Bbasisobj);
# bval     <- exp(basismat*bvec);
# end
#
# function Dbval <- fun_Dexplinear(t, bvec, Bbasisobj)
# nbasis    <- length(bvec);
# if isa_fdPar(Bbasisobj)
#     Bbasisobj <- getbasis(Bbasisobj);
# end
# basismat  <- eval_basis(t, Bbasisobj);
# bval      <- exp(basismat*bvec);
# Dbval     <- basismat.*repmat(bval,1,nbasis);
# end

daytime   <- seq(0.5,364.5,len=365)  #  time in days
daytime   <- daytime*12/365  #  measure time in months
dayrange  <- c(0,12)
#  set up block averages for Canadian temperature data,
#  available from the fda package
tempav <- CanadianWeather$dailyAv[,,"Temperature.C"]
#  center the time of the year on winter
winterind  <- c(183:365,1:182)
tempData   <- tempav[winterind,]
#  select the data for Montreal
station    <- 12 
tempData12 <- matrix(tempData[,station],365,1)

#  set up TempDataList

TempDataList1 <- list(argvals=as.matrix(daytime), y=as.matrix(tempData12))
TempDataList  <- vector("list",1)
TempDataList[[1]]    <- TempDataList1

#  set up fourier basis for representing the fit curve

nbasis <- 51
daybasis <- create.fourier.basis(dayrange, nbasis)

TempBasisList <- vector("list",1)
TempBasisList[[1]] <- daybasis

#  set up a constant basis 

Cbasis <- create.constant.basis(dayrange)

#  List objects for each homogeneous term

#  Model: no forcing, harmonic only, fourier coefficients

nFbasis <- 3;
Fbasis  <- create.fourier.basis(dayrange, nFbasis)

Fparvec <- c(0.3, rep(1,nFbasis-1))
# Fields:            funobj   parvec    estimate  variable deriv. factor
Xterm1 <- make.Xterm(Fbasis,  Fparvec,  TRUE,     1,       1,     -1);

XList <- vector("list", 1)
XList[[1]] = Xterm1

#  struct objects for each forcing term

FList <- NULL

#  set up variable list object TempVariableList

TempVariableList <- make.Variable(name="Temperature", order=3, XList=XList, FList=FList)

# set up  model list object TempModelList

TempModelList <- vector("list",1)
TempModelList[[1]] <- TempVariableList

# check the model list object

TempModelcheck <- checkModel(TempBasisList, TempModelList)
TempModelList  <- TempModelcheck$modelList
nparam         <- TempModelcheck$nparam

print(paste("Number of parameters =",nparam))

# An evaluation of the criterion at the initial values

rhoVec <- 0.5

#  This command causes Data2LD to set up and save the tensors

Data2LDResult <- Data2LD(TempDataList, TempBasisList, TempModelList, rhoVec)

MSE        <- Data2LDResult$MSE        # Mean squared error for fit to data
DpMSE      <- Data2LDResult$DpMSE      #  gradient with respect to parameter values

round(MSE,3)
round(DpMSE,3)

#  set constants for estimation algorithm Data2LD.opt

dbglev  <-  1    
iterlim <- 50    
convrg  <- 1e-5  
#  define rhovec using the logit function
gammavec <- seq(0,7,1)
rhoVec   <- exp(gammavec)/(1+exp(gammavec))
nrho     <- length(rhoVec)
#  Matrices to hold results
dfesave <- matrix(0,nrho,1)
gcvsave <- matrix(0,nrho,1)
MSEsave <- matrix(0,nrho,1)
ISEsave <- matrix(0,nrho,1)
thesave <- matrix(0,nrho,nparam)

#  Initialize coefficient list
TempModelList.opt <- TempModelList
#  Loop through rho values
for (  irho  in  1 : nrho ) {
  rhoVeci <- rhoVec[irho]
  print(paste('Rho <- ',round(rhoVeci,5)))
  OptList <- Data2LD.opt(TempDataList, TempBasisList, TempModelList.opt, 
                         rhoVeci, convrg, iterlim, dbglev)
  theta.opti <- OptList$thetastore
  TempModelList.opt <- modelVec2List(TempModelList, theta.opti)    
  Data2LDResult <- Data2LD(TempDataList, TempBasisList, TempModelList.opt, rhoVeci)
  thesave[irho,] <- theta.opti
  dfesave[irho]  <- Data2LDResult$df
  gcvsave[irho]  <- Data2LDResult$gcv
  MSEsave[irho]  <- Data2LDResult$MSE
  ISEsave[irho]  <- Data2LDResult$ISE
} 

# display degrees of freedom and gcv values
print('    rho      df         gcv')
print(round(cbind(rhoVec, dfesave, gcvsave),3))

# Evaluate the fit for parameter values at highest rho value

irho <- 8
theta <- matrix(thesave[irho,],nparam,1)
TempModelList.opt <- modelVec2List(TempModelList.opt, theta)
rhoi <- rhoVec[irho]
Data2LDList <- Data2LD(TempDataList, TempBasisList, TempModelList.opt, rhoi)
MSE       <- Data2LDList$MSE 
df        <- Data2LDList$df
gcv       <- Data2LDList$gcv 
ISE       <- Data2LDList$ISE 
Var.theta <- Data2LDList$Var.theta
print(paste("MSE <- ", round(Data2LDResult$MSE,4)))
print(paste("df  <- ", round(Data2LDResult$df, 4)))
print(paste("gcv <- ", round(Data2LDResult$gcv,4)))

#  Set up functional data object and fine mesh for plotting variable

XfdParList <- Data2LDList$XfdParList
tempfdPar <- XfdParList[[1]]
tempfd    <- tempfdPar$fd
tfine     <- seq(0,12,len=101)
tempfine  <- eval.fd(tfine, tempfd)

#  Plot the fit to the data
par(mfrow=c(1,1),ask=FALSE)
plot(tfine, tempfine, type="l", lwd=2, xlim=c(0,12), ylim=c(-15,25),
     xlab="Time (days)", ylab="Temperature (deg C)")
points(daytime, tempData12, pch="o") 

# plot the optimal parameter values as functions of rho
#  plot flow of beta(t) parameters

matplot(rhoVec, thesave, type="b", lwd=2, 
        xlab="rho", ylab="beta coefficients")

#  plot \beta with confidence intervals

betafd   <- fd(theta, Fbasis)
betavec  <- eval.fd(daytime,betafd)
basismat <- (betavec %*% matrix(1,1,3)) * eval.basis(daytime, Fbasis)
CI.beta  <- 2*sqrt(diag(basismat %*% Var.theta %*% t(basismat)))
betamat  <- cbind(betavec, betavec+CI.beta, betavec-CI.beta)
par(mfrow=c(1,1))
matplot(daytime, betamat , type="l", lty=c(1,4,4), col=1, lwd=2, 
       xlab="Time (days)", ylab="beta", xlim=c(0,12))

#  plot fit to derivative, lower panel in Figure 12.10 in the book

D1dayvec = eval.fd(daytime, tempfd, 1)
D3dayvec = eval.fd(daytime, tempfd, 3)

D3dayfit = -betavec * D1dayvec

D3mat = cbind(D3dayfit,D3dayvec)
matplot(daytime, D3mat, type="l", lty=c(1,4), col=1,
     xlab="Time (days)", ylab="D3 Temperature", xlim=c(0,12))

