DBSarrayFn <- function(nXbasisw, nWbasisw, 
                      nUbasisj, nAbasisj, 
                      nrep, Avecj, Ucoefj, BAtens) {
  BAtens <- as.vector(BAtens)
  DSvecC <- rep(0,nXbasisw*nrep*nWbasisw)
  inprodDBList <- .C("DBSarrayFn", as.integer(nXbasisw), 
                     as.integer(nWbasisw), 
                     as.integer(nUbasisj), 
                     as.integer(nAbasisj), 
                     as.integer(nrep),
                     as.double(Avecj),
                     as.double(Ucoefj),
                     as.double(BAtens), 
                     as.double(DSvecC),
                     PACKAGE = "Data2LD")
  DBSarrayC <- inprodDBList[[9]]
  return(DBSarrayC)
}
