DASarrayFn <- function(nXbasisw, nWbasisw, nUbasisj, nAbasisj, 
                      nrep, Bvecw, Ucoefj, BAtens) {
  BAtens <- as.vector(BAtens)
  DSvecC <- rep(0,nXbasisw*nrep*nAbasisj)
  inprodDAList <- .C("DASarrayFn", 
                     as.integer(nXbasisw), 
                     as.integer(nWbasisw), 
                     as.integer(nUbasisj), 
                     as.integer(nAbasisj), 
                     as.integer(nrep),
                     as.double(Bvecw),
                     as.double(Ucoefj),
                     as.double(BAtens), 
                     as.double(DSvecC),
                     PACKAGE = "Data2LD")
  DASarrayC <- inprodDAList[[9]]
  return(DASarrayC)
}

