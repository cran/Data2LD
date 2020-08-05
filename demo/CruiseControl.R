# Cruise control: A Two-variable Feedback Loop
#
# A driver starts his car and (1) accelerates to the 60 km/h speed limit
# holding on most main streets of Ottawa, (2) reduces speed to 40 km/h in a
# school zone, (3) enters a controlled-access boulevard with an 80 km/h
# limits and finally (4) against enters a typical principal Listeet.
# Operating on a snowy winter morning, the driver takes 8 seconds to reach
# 60 km/h, or 4/30 seconds per unit change in speed.
#
# The speed of the car under acceleration is modelled by a first order
# constant coefficient equation:
#
# The feedback equations: 
# The two differental equations specifying the feedback model are:
#
# $$DS(t) <- \beta_{11} S(t) + \beta_{12} C(t)$$
# $$DC(t) <- \beta_{21} S(t) + \beta_{22} C(t) + \alpha_[[2]] S.0(t)$$
#
# The right side term
# $\beta_{12} C(t)$ in the first equation is the contribution of the 
# control variable $C$ to the change in speed $DS$.
# In the second term on the right side the variable $S_0$ 
# is the set-point function specifying the target speed, and it is
# called a forcing function.
#
# We will use a special case of these equations to generate some 
# simulated data and then estimate the parameters defining the 
# equation using these simulated data.
#
# The specialized equations are
#
# $$DS(t) = -S(t)   + C(t)/4$$
# $$DC(t) =  S_0(t) - S(t)$$
#
# These equations correspond to:
# $\beta_{11}$ = -1
# $\beta_{12}$ =  1/4
# $\beta_{21}$ = -1
# $\beta_{22}$ =  0
# $\alpha_2$   =  1
#
# We see that when (1) the speed $S(t)$ is less than the set point value
#                      S_0(t), the control variable increases to a 
#                      positive value and forces the speed to increase, and
#                  (2) the speed $S(t)$ is greater than the set point value
#                      S_0(t), the control variable decreases to a 
#                      negative value and forces the speed to decrease.
# The value $\beta_{22}$ =  0 implies that the controller responds
# instantly to a change in the difference $S_0(t) - S(t)$.
#
# We will try out two models.  In the first we estimate the four 
# parameters $\beta_{11}$, $\beta_{12}$, $\beta_{21}$ and $\alpha_2$.
# In the second we coerce $\beta_{21} - \alpha_2$ to be 0.  In both
# models we will set up $\beta_{22}$ as a parameter but fix its value
# to be 0.

# Last modified 17 April 2020

#  Define the problem:
#  Set up the time span, a set of observation times, and a fine grid for
#  plotting:

T     <- 80  #  seconds
rng   <- c(0,T)
n     <- 41
nfine <- 501
tobs  <- seq(0,T,len=n)
tfine <- seq(0,T,len=nfine)

# Set up the set-point forcing function.
#  UfdList is a list array having the same two-level listListucture as 
#  AwtList, but the contents of the lists are functional data objects
#  specifying the external input to the system. If a list is empty, it is
#  assumed that there is no forcing of the corresponding equation.
#
# The set-point function uses an order 1 step function B-spline basis.  The 
# knots are placed at the points where the set-point changes values.

steporder  <- 1  #  step function basis
stepnbasis <- 4  #  four basis functions
stepbreaks <- c(0,20,40,60,80)
stepbasis  <- create.bspline.basis(rng, stepnbasis, steporder, stepbreaks)
stepcoef   <- c(60,40,80,60)  # target speed for each step
SetPtfd    <- fd(stepcoef, stepbasis)  # define the set point function

# Set up cruiseModelList
#  The total number of coefficients defining the estimated coefficient 
#  functions is three because each coefficient function has only one 
#  coefficient, and only three of them are estimated.

#  set up a constant basis over this range

conbasis  <- create.constant.basis(rng)
confdPar  <- fdPar(conbasis)

# Solving the equations for known parameter values and initial conditions.
# In order to simulate data, we need to know the true values of .S(t).
# and .C(t). at the time points where the process is observed.  We also 
# need to specify the initial state of the system at time 0, which we define 
# to be zero for both variables.  We get the solution by using an initial 
# value approximation algorithm, which is here the Runge-Kutta fourth order
# method coded in Matlab"s function |ode45|.
# The function |cruuise1| evaluates the right side of the equations at a 
# set of time values.  
#
# Here is a function that evaluates the right side of the equation for a 
# time value:

cruise0 <- function(t, y, parms) {
  DSvec    <- matrix(0,2,1)
  Uvec     <- eval.fd(t, parms$SetPtfd)
  DSvec[1] <- -y[1] + y[2]/4
  DSvec[2] <-  Uvec - y[1]
  return(list(DSvec=DSvec))
}

# We first have a look at the solution at a fine mesh of values by solving
# the equation for the points in |tfine| and plotting them

y0 <- matrix(0,2,1)
parms = list(SetPtfd=SetPtfd)
ytrue = lsoda(y0, tfine[1:500], cruise0, parms)
ytrue = rbind(ytrue,matrix(c(80,60,240),1,3))

#  Plot the true solution

par(mfrow=c(2,1))
#  speed panel
plot(tfine, ytrue[,2], type="l", lwd=2, ylim=c(0,100), ylab="Speed S (mph)")
lines(c( 0,20),c(60,60),type="l",lty=2)
lines(c(20,20),c(60,40),type="l",lty=2)
lines(c(20,40),c(40,40),type="l",lty=2)
lines(c(40,40),c(40,80),type="l",lty=2)
lines(c(40,60),c(80,80),type="l",lty=2)
lines(c(60,60),c(80,60),type="l",lty=2)
lines(c(60,80),c(60,60),type="l",lty=2)
#  control panel
plot(tfine, ytrue[,3], type="l", lwd=2, 
     xlab="Time (mins)", ylab="Control level C")

#  interpolate at 41 equally spaced time points the values of these variables

speedResult    <- approx(tfine, ytrue[,2], seq(0,80,len=41))
controlResult  <- approx(tfine, ytrue[,3], seq(0,80,len=41))

# Simulate noisy data at n observation points
# We simulate data by adding a random zero-mean Gaussian deviate to each 
# curve value.  The deviates for speed have a standard deviation 2 speed 
# units, and those for the control level have a standard deviation of 8.

sigerr <- 2
yobs     <- matrix(0,length(tobs),2)
yobs[,1] <- as.matrix(  speedResult$y + rnorm(41)*sigerr)
yobs[,2] <- as.matrix(controlResult$y + rnorm(41)*sigerr*4)

# plot the data along with the true solution:

par(mfrow=c(2,1))
plot(tfine, ytrue[,2], type="l", ylab="Speed")
lines(c(0,T), c(60,60), lty=3)
points(tobs, yobs[,1], pch="o")

plot(tfine, ytrue[,3], type="l", xlab="Time (mins)", ylab="Control level")
lines(c(0,T), c(240,240), lty=3)
points(tobs, yobs[,2], pch="o")

#  Define cruiseDataList, and insert these into structs into the corresponding
#  lists.

cruiseDataList1 <- list(argvals=tobs, y=yobs[,1])
cruiseDataList2 <- list(argvals=tobs, y=yobs[,2])

cruiseDataList <- vector("list",2)
cruiseDataList[[1]] <- cruiseDataList1
cruiseDataList[[2]] <- cruiseDataList2

# Define cruiseBasisList containing the basis system for each varible.
# We also have to provide a basis system for each variable that is large
# enough to allow for any required sharp curvature in the solution to the
# differential equation system.  
#
# First we set the order of the B-spline basis to 5 so that the first 
# derivative will be smooth when working with a second order derivative in
# the penalty term.  Then we position knots at 41 positions where we willl
# simulate noisy observations.  We use the same basis for both variables
# and load it into a list array of length 2.

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

#  Now set up the list objects for each of the two cells in 
#  list cruiseModelList

# List object for the speed equation: A term for speed and a term for control, but no forcing

SList.XList = vector("list",2)
#  Fields:                     funobj    parvec  estimate  variable deriv. factor
SList.XList[[1]] <- make.Xterm(confdPar, 1,      TRUE,     1,       0,     -1)
SList.XList[[2]] <- make.Xterm(confdPar, 1/4,    TRUE,     2,       0,      1)
SList.FList = NULL
SList = make.Variable("Speed", 1, SList.XList, SList.FList)

# List object for the control equation:  a term for speed, and a zero-multiplied term for control
# plus a term for the forcing function SetPtfd

CList.XList <- vector("list",2)
#  Fields:                     funobj    parvec  estimate  variable deriv. factor
CList.XList[[1]] <- make.Xterm(confdPar, 1,      TRUE,     1,       0,     -1)
CList.XList[[2]] <- make.Xterm(confdPar, 0,      FALSE,    2,       0,      1)
CList.FList <- vector("list",1)
#  Fields:                     funobj    parvec  estimate  Ufd      factor
CList.FList[[1]] <- make.Fterm(confdPar, 1,      TRUE,     SetPtfd, 1)
CList <- make.Variable("Control", 1, CList.XList, CList.FList)
#  Now set up the Listuct objects for each of the two lists in 
#  list array cruiseVariableList

#  List array for the whole system

cruiseModelList <- vector("list",2)
cruiseModelList[[1]] <- SList
cruiseModelList[[2]] <- CList

#  check the system specification for consistency

cruiseModelCheckList <- checkModel(cruiseBasisList, cruiseModelList)
cruiseModelList <- cruiseModelCheckList$modelList
nparam <- cruiseModelCheckList$nparam

print(paste("total number of parameters = ", nparam))

# A preliminary evaluation of the function and its derivatives

rhoVec <- 0.5*matrix(1,1,2)

Data2LDList <- Data2LD(cruiseDataList, cruiseBasisList, cruiseModelList, rhoVec)

print(Data2LDList$MSE)
print(Data2LDList$DpMSE)
print(Data2LDList$D2ppMSE)


#  Note that, for the parameter that is fixed, the gradient element is 0, and
#  the row and column of the hessian matrix are 0 except for the diagonal
#  value, which is 1.

#  Evaluate at the corrent parameter values, which in this case are
#  the right values since the data are without error

#  set up a loop through a series of values of rho
#  We know, because the signal is smooth and the data are rough, that the 
#  optimal value of rho will be rather close to one, here we set up a 
#  range of rho values using the logistic transform of equally spaced
#  values between 0 and 5.
#  For each value of rho we save the degrees of freedom, the gcv value,
#  the stop sum of squares for each equation, the mean squared stops for 
#  the parameters, and the parameter values.

#  set up a loop through a series of values of rho
#  We know, because the signal is smooth and the data are rough, that the 
#  optimal value of rho will be rather close to one, here we set up a 
#  set of three rho values using the logistic transform of equally spaced
#  values 2, 5, and 8, corresponding to rho values of 0.8808, 0.9933 and
#  0.9997.
#  For each value of rho we save the degrees of freedom, the gcv value,
#  the error sum of squares for each equation, the mean squared errors for 
#  the parameters, and the parameter values.

Gvec    <- c(0:7)
nrho    <- length(Gvec)
rhoMat  <- matrix(exp(Gvec)/(1+exp(Gvec)),nrho,2)
dfesave <- matrix(0,nrho,1)
gcvsave <- matrix(0,nrho,1)
MSEsave <- matrix(0,nrho,2)
thesave <- matrix(0,nrho,nparam)

conv    <- 1e-4  #  convergence criterion for mean square error values
iterlim <- 20    #  limit on number of iterations
dbglev  <-  1    #  displays one line of results per iteration

cruiseModelList.opt <- cruiseModelList  #  initialize the optimizing coefficient values

#  loop through rho values, with a pause after each value

for (irho in 1:nrho) {   
  rhoVeci <- rhoMat[irho,]
  print(paste(" ------------------  rhoVeci <- ", round(rhoVeci[1],4), 
              " ------------------"))
  Data2LDOptList <- Data2LD.opt(cruiseDataList, cruiseBasisList, cruiseModelList.opt, 
                                rhoVeci, conv, iterlim, dbglev)
  theta.opti <- Data2LDOptList$theta  # optimal parameter values
  cruiseModelList.opt <- modelVec2List(cruiseModelList, theta.opti) 
  #  evaluate fit at optimal values and store the results
  Data2LDList    <- Data2LD(cruiseDataList, cruiseBasisList, cruiseModelList.opt, rhoVeci)
  thesave[irho,]  <- theta.opti
  dfesave[irho]   <- Data2LDList$df
  gcvsave[irho]   <- Data2LDList$gcv
  x1fd            <- Data2LDList$XfdParList[[1]]$fd
  x1vec           <- eval.fd(tobs,  x1fd)
  msex1           <- mean((x1vec - speedResult$y)^2)
  x2fd            <- Data2LDList$XfdParList[[2]]$fd
  x2vec           <- eval.fd(tobs,  x2fd)
  msex2           <- mean((x2vec - controlResult$y)^2)
  MSEsave[irho,1] <- msex1
  MSEsave[irho,2] <- msex2
}

ind <- 1:nrho
#  print df, gcv and MSEs
print("    rho      df        gcv        RMSE:")
print(cbind(round(rhoMat[ind,1],4), round(dfesave[ind],1), 
            round(gcvsave[ind],1),  round(sqrt(MSEsave[ind,1]),2)))

#  plot parameters as a function of \rho

matplot(rhoMat[ind,1], thesave[ind,], type="b", xlab="rho", ylab="theta(rho)")

rho.opt   <- rhoMat[nrho,]
theta.opt <- thesave[nrho,]

#  convert the optimal parameter values to optimal cruiseModelList

cruiseModelList.opt <- modelVec2List(cruiseModelList, theta.opt)

#  evaluate the solution at the optimal solution

DataLDList <- Data2LD(cruiseDataList, cruiseBasisList, cruiseModelList.opt, rho.opt)

#  display parameters with 95% confidence limits

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

#  Plot the estimated solutions, as estimated from the data, rather
#  than from the initial value estimates in the above code

XfdParList <- Data2LDList$XfdParList
Xfd1 <- XfdParList[[1]]$fd
Xfd2 <- XfdParList[[2]]$fd

Xvec1 <- eval.fd(tfine, Xfd1)
Xvec2 <- eval.fd(tfine, Xfd2)

Uvec  <- eval.fd(tfine, SetPtfd)

par(mfrow=c(2,1))
cruiseDataList1 <- cruiseDataList[[1]]
plot(tfine, Xvec1, type="l", xlim=c(0,80), ylim=c(0,100), 
     ylab="Speed",
     main=paste("RMSE =",round(sqrt(MSEsave[3,1]),4)))
lines(tfine, ytrue[,2], lty=2)
points(cruiseDataList1$argvals, cruiseDataList1$y, pch="o")
lines(tfine, Uvec, lty=2)

cruiseDataList2 <- cruiseDataList[[2]]
plot(tfine, Xvec2, type="l", xlim=c(0,80), ylim=c(0,400),
     xlab="Time (sec)", ylab="Control",
     main=paste("RMSE =",round(sqrt(MSEsave[3,2]),4)))
lines(tfine, ytrue[,3], lty=2)
points(cruiseDataList2$argvals, cruiseDataList2$y, pch="o")


