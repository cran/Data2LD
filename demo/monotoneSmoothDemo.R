###  Monotone smoothing of growth data

#  This demonstration is designed to highlight the fact that functions
#  with specific structural constraints are often best defined as 
#  solutions of differential equations.
#  The most elementary differential equation is 
#                       Dx(t) = beta x(t)
#  and this equation defines either exponential growth or exponential
#  decay.  In both cases the function x(t) is either always positive or 
#  always negative.  More generally, if parameter $beta$ is a function,
#  it still follows that the function's values are all on one side of zero.
#  This implies that the second ordr linear equation
#                      D2x(t) = beta Dx(t)
#  and its generalization define functions that are everwhere either
#  increasing (Dx >0) or decreasing (Dx < 0).  These functions are 
#  therefore strictly monotonic, and, beta is a function, define the
#  complete class of twice differentiable strictly monotonic functions.

#  There are many applications which require this property, and here we
#  illustrate one of them, the estimation of a smooth height function for 
#  growing children.

#  As a side issue, we also illustrate that fixing a parameter that is
#  poorly defined by the data will result in considerable improvement
#  in speed of computing with negligible loss of fitting power.

## plot the data

agerng = c(1,18)
age = growth$age

##  Define the cell array containing the data

igirl = 4
height = as.matrix(growth$hgtf[,igirl])

growthList$argvals = age
growthList$y       = height

growthDataList = vector("list", 1)
growthDataList[[1]] = growthList

## construct a basis for representing a growth curve

growthbreaks = age
growthnbasis = 33
growthbasis  = create.bspline.basis(agerng, growthnbasis, 4, growthbreaks)

growthtfine = seq(1,18,len=501)
growthyfine = eval.basis(growthtfine, growthbasis)

XbasisList = vector("list",1)
XbasisList[[1]] = growthbasis

## Set up the cell array growthModel as a second order linear equation

#  the coefficient beta is unconstrained and varies over time

betanbasis = 7
betabasis  = create.bspline.basis(agerng, betanbasis)
betafdPar  = fdPar(betabasis)

#  The damping coefficient is estimated but the reaction coefficient is 0.

# Fields:              funobj      parvec         est.  var. deriv. factor
XtermList = make.Xterm(betafdPar,  matrix(0,7,1), TRUE, 1,   1,     1)
XList    = vector("list",1)
XList[[1]] = XtermList

#  There is no forcing term for this model.

FList    = vector("list",1)
FList[[1]] = NULL

# List object for the single variable in the equation

# Fields:                 name      order XList 
growthList = make.Variable("Height", 2,   XList)

#  Set up the model list

growthModel = vector("list",1)
growthModel[[1]] = growthList

#  check and print the initial model

checkModelList = checkModel(XbasisList, growthModel, TRUE)
growthModel = checkModelList$modelList 
nparam      = checkModelList$nparam

print(paste("Number of parameters = ",nparam))

## print the initial parameter vector

# The parameters in the model are arranged in the vector in two groups:
# first:  all the parameters associated with the homogeneous terms
#         in the differential equations, in order of the variables and
#         the terms within the variables
# second: all the parameters associated with the forcing terms
#         in the differential equations, in order of the variables and
#         the terms within the variables

thetavec = modelList2Vec(growthModel, nparam)

# print the parameter vector

print(t(thetavec))

## An evaluation of the criterion at the initial values

rho = 0.5  #  light smoothing

#  this command causes Data2LD to set up and save the tensors

Data2LDList = Data2LD(growthDataList, XbasisList, growthModel, rho)

MSE   = Data2LDList$MSE
DMSE  = Data2LDList$DpMSE
D2MSE = Data2LDList$D2ppMSE

#  print the initial criterion, gradient and hessian matrix

MSE
DMSE
D2MSE

##  Set up the optimisation

dbglev   =  1    #  debugging level
iterlim  = 50    #  maximum number of iterations
conv     = c(1e-4, 1e-3)  #  convergence criterion

gammavec = 0:8
rhoVec   = exp(gammavec)/(1+exp(gammavec))
nrho   = length(rhoVec)
dfesave = matrix(0,nrho,1)
gcvsave = matrix(0,nrho,1)
MSEsave = matrix(0,nrho,1)
thesave = matrix(0,nrho,nparam)

## Optimization of the criterion

growthModel.opt = growthModel

for (irho in 1:nrho) {
    rhoi = rhoVec[irho]
    optList = Data2LD.opt(growthDataList, XbasisList, growthModel.opt, 
                             rhoi, conv, iterlim, dbglev)
    theta.opti <- optList$thetastore
    growthModel.opt = modelVec2List(growthModel, theta.opti)
    Data2LDList = Data2LD(growthDataList, XbasisList, growthModel.opt, rhoi)
    MSE     = Data2LDList$MSE   
    DpMSE   = Data2LDList$DpMSE
    D2ppMSE = Data2LDList$D2ppMSE
    XfdList = Data2LDList$XfdList
    df      = Data2LDList$df
    gcv     = Data2LDList$gcv
    thesave[irho,] = t(theta.opti)
    dfesave[irho]   = df
    gcvsave[irho]   = gcv
    MSEsave[irho]   = MSE
}

## print degrees of freedom and gcv values

ind <- 1:nrho
#  print df, gcv and MSEs
print("    rho      df        gcvE:")
print(cbind(round(rhoVec[ind],4), 
            round(dfesave[ind],3), 
            round(gcvsave[ind],3)))

## Evaluate the fit for parameter values at highest rho value

irho = 8  # the value with the lowest gcv index

#  set up the final model 

theta    = thesave[irho,]
growthModel.opt = modelVec2List(growthModel, theta)

##  print the model

printModel(growthModel.opt, paste("The final model for girl ",igirl))

## print the optimal parameter vector

# The parameters in the model are arranged in the vector in two groups:
# first:  all the parameters associated with the homogeneous terms
#         in the differential equations, in order of the variables and
#         the terms within the variables
# second: all the parameters associated with the forcing terms
#         in the differential equations, in order of the variables and
#         the terms within the variables

##  print MSE, degrees of freedom and gcv index for this model

rho      = rhoVec[irho]
Data2LDList    = Data2LD(growthDataList, XbasisList, growthModel.opt, rho)
    MSE        = Data2LDList$MSE   
    DpMSE      = Data2LDList$DpMSE
    D2ppMSE    = Data2LDList$D2ppMSE
    XfdParList = Data2LDList$XfdParList
    df         = Data2LDList$df
    gcv        = Data2LDList$gcv
    ISE        = Data2LDList$ISE
    Var.theta  = Data2LDList$Var.theta
print(paste("MSE = ", MSE))
print(paste("df  = ", df))
print(paste("gcv = ", gcv))

##  plot the estimated time-varying homogeneous parameter

betafd   = fd(theta,betabasis)
betafine = eval.fd(growthtfine, betafd)

plot(growthtfine, betafine, type="l", xlim=c(1,18), ylim=c(-1.6,3.3),
     xlab="Time (years)", ylab="beta(t)",
     main=paste("RMSE =",round(sqrt(MSE),1)))

##  print the data and the fit of the growth curve to the data

growthfd     = XfdParList[[1]]$fd
growthfine   = eval.fd(growthtfine, growthfd)

plot(growthtfine, growthfine, type="l", xlim=c(1,18), ylim=c(60,200), 
     xlab="Time (years)", ylab="height (cm)",
     main=paste("Girl",igirl,",  Std. Error = ",round(sqrt(gcv),2)," cm"))
points(age, height)

##  print the parameter values along with their 95# confidence limits

stderr = sqrt(diag(Var.theta))
plot(1:nparam, theta, type="b", xlim=c(1,7), ylim=c(-10,10),
     xlab="Parameter number", ylab="Parameter") 
lines(1:nparam, theta + 2*stderr, lty=2)
lines(1:nparam, theta - 2*stderr, lty=2)
 
## print the estimated linear differential operator and its 
# fit to the second derivative

Dgrowthfine  = eval.fd(growthtfine, growthfd, 1)
D2growthfine = eval.fd(growthtfine, growthfd, 2)
D2fitfine    = Dgrowthfine*betafine

plot(growthtfine, D2growthfine, type="l", xlim=c(1,18), ylim=c(-7,2),
     xlab="Time (years)", ylab = "Height Acceleration (cm/y/y)")
lines(growthtfine, D2fitfine, lty=2) 
lines(c(1,18), c(0,0), lty=4)

##  plot the evolution of the parameters over the values of rho

matplot(rhoVec, thesave, type="b", xlab="rho", ylab="parameter beta")

##  Fix the final coefficient at zero
#  The wide confidcence limits suggest that the data do not provide
#  strong information about the value of the final coefficient.
#  This is because the growth curve is flat over the interval
#  between the final pair of knots. 
#  Other girls will not print this, however.

#  fix the final coefficient value at 3

estim = rep(TRUE,7)
estim[7] = FALSE
parvec = matrix(0,7,1)
parvec[7] = 3

# Fields:              funobj      parvec      est.   var. deriv. factor
XtermList = make.Xterm(betafdPar,  parvec,     estim, 1,   1,     1)
XList    = vector("list",1)
XList[[1]] = XtermList

FList    = vector("list",1)
FList[[1]] = NULL

# List object for the single variable in the equation

# Fields:                 name      order XList 
growthList = make.Variable("Height", 2,    XList)

growthModel = vector("list",1)
growthModel[[1]] = growthList

checkModelList = checkModel(XbasisList, growthModel, TRUE)
growthModel = checkModelList$modelList 
nparam      = checkModelList$nparam

print(paste("Number of parameters = ",nparam))

#  Rerun the analysis using the code above
#  You will notice that the convergence is now very fast for each 
#  value of rho.
