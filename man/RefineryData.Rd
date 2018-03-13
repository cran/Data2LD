\name{RefineryData}
\alias{RefineryData}
\docType{data}
\title{
  Reflux and tray level in a refinery
}
\description{
  194 observations on reflux and "tray 47 level" in a distallation
  column in an oil refinery.  
}
\usage{
  RefineryData
}
\format{
  A data.frame with the following components:  
  \itemize{
    \item{Time:}{   observation time 0:193}
    \item{Tray.level:}{ reflux flow centered on the mean of the first 60 observations}
    \item{Valve.setting:}{ tray 47 level centered on the mean of the first 60 observations}
  }
}
\source{
  Ramsay, James O., and Silverman, Bernard W. (2005), \emph{Functional
    Data Analysis, 2nd ed.}, Springer, New York, p. 4, Figure 1.4, and
  chapter 17.   
}
\examples{
#  input the data
TimeData <- RefineryData[,1]
TrayData <- RefineryData[,2]
ValvData <- RefineryData[,3]
#  plot the data
par(mfrow=c(2,1))
plot(TimeData, TrayData, type="p") 
lines(c(67,67), c(0,4.0), type="l")
plot(TimeData, ValvData, type="p")
lines(c(67,67), c(0,0.5), type="l")
}
\keyword{datasets}

