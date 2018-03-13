fdascript.force.plot <- function(Aknots, Acoef, Aname, vertrng) {
  nAknots <- length(Aknots)
  k <- 1
  plot( c(Aknots[k],  Aknots[k+1]), c(Acoef[k], Acoef[k]), type="l", xlim=c(0,230), ylim=c(-1.5,1.5))
  lines(c(Aknots[k+1],Aknots[k+1]), c(Acoef[k], Acoef[k+1]))
  for (k in 2:(nAknots-2)) {
    lines(c(Aknots[k],  Aknots[k+1]), c(Acoef[k],   Acoef[k]))
    lines(c(Aknots[k+1],Aknots[k+1]), c(Acoef[k],   Acoef[k+1]))
  }
  lines(c(Aknots[nAknots-1],Aknots[nAknots]), c(Acoef[nAknots-1], Acoef[nAknots-1]))
  lines(c(Aknots[1],        Aknots[nAknots]), c(0,0), lty=2)
}
