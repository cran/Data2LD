DRarrayFn <- function(nXbasisw, nWbasisw, nXbasisx, nWbasisx, 
                     Bvecx, Btens) {
  DRvec <- rep(0,nXbasisw*nXbasisx*nWbasisw)
  Btens <- as.vector(Btens)
  inprodDList <- .C("DRarrayFn", as.integer(nXbasisw), 
                     as.integer(nWbasisw), 
                     as.integer(nXbasisx), 
                     as.integer(nWbasisx), 
                     as.double(Bvecx), 
                     as.double(Btens), 
                     as.double(DRvec),
                     PACKAGE="Data2LD")
  DRarrayC <- array(inprodDList[[7]],c(nXbasisw,nXbasisx,nWbasisw))
  return(DRarrayC)
}
