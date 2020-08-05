\name{fdascript}
\alias{fdascript}
\docType{data}
\title{
  Positions of the tip of a pen in centimeters while writing in cursive script the 
  letters "fda" twenty times.
}
\description{
  The writing was on a horizontal surface with the position of the pen recorded 400 
  times per second with a background noise level of about 0.05 centimeters.
}
\usage{
  fdascript
}
\format{
  An array object \code{fdascript} with dimensions 1401, 20, and 2.  
  The first dimension is the number of observations in each record, the second is 
  the records, and the third dimension is for the X and Y coordinate.
}
\details{
  The times are in seconds over an interval of 2.3 seconds, but in the analysis 
  seconds are converted to "centiseconds" by multiplication by 100 in order to 
  ensure numerical accuracy in the calculations.
}
\source{
  The data were published and analyzed in the 1997 and 2005 editions of the book by 
  J. O. Ramsay and B. W. Silverman, Functional Data Analysis, and made available 
  in the R package fda.}
\examples{
#  load the first record
fdascriptX <- fdascript[,1,1]
fdascriptY <- fdascript[,1,2]
#  Define the observation times in 100ths of a second
centisec <- seq(0,2.3,len=1401)*100
#  plot the script
plot(fdascriptX, fdascriptY, type="l")
}
\keyword{datasets}
