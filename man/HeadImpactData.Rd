\name{HeadImpactData}
\alias{HeadImpactData}
\docType{data}
\title{
  Acceleration of brain tissue before and after a blow to the head.
}
\description{
  The data are the acceleration of brain tissue in mm/msec in a cadaver before 
  and after five blows to the cranium measured over 60 milliseconds.  
  The times of the blows within this interval were at 14 msec.
}
\usage{
  HeadImpactData
}
\format{
  A 133 by 3 matrix of real numbers.  The first column contains row numbers, 
  the second the times of observations, and third the tissue acceleration in millimeters
  per second per second.
}
\source{
  The data are widely reported in various books on data smoothing.  The first report 
  was in Wolfgang Hartle's 1990 book Appplied Nonparamemtric Regression.  
  The original experiment was conducted at the University of Heidleberg, 
  and as far as we know, seems to have been unpublished.
}
\examples{
HeadImpactTime <- HeadImpactData[,2]  #  time in milliseconds
HeadImpact     <- HeadImpactData[,3]  #  acceleration in millimeters/millisecond^2
HeadImpactRng  <- c(0,60) # Define range time for estimated model
# plot the data along with a unit pulse
plot(HeadImpactTime, HeadImpact, type="p", xlim=c(0,60),  ylim=c(-1.0,1.5),
     xlab="Time (milliseconds)", ylab="Acceleration (mm/msec^2)")
lines(c( 0,60), c(0,0), lty=3)
lines(c(14,14), c(0,1), lty=2)
lines(c(15,15), c(0,1), lty=2)
lines(c(14,15), c(1,1), lty=2)
}
\keyword{datasets}
