SmatFn <- function(nXbasisw, nWbasisw, nUbasisj, nAbasisj, nrep, 
                  Bvecw, Avecj, Ucoefj, BAtens) {
  BAtens <- as.vector(BAtens)
  SvecC <- rep(0,nXbasisw*nrep)
  inprodSList <- .C("SmatFn", as.integer(nXbasisw), 
                     as.integer(nWbasisw), 
                     as.integer(nUbasisj), 
                     as.integer(nAbasisj), 
                     as.integer(nrep),
                     as.double(Bvecw), 
                     as.double(Avecj),
                     as.double(Ucoefj),
                     as.double(BAtens), 
                     as.double(SvecC),
                     PACKAGE = "Data2LD")
  SmatC <- matrix(inprodSList[[10]],nXbasisw,nrep)
  return(SmatC)
}
