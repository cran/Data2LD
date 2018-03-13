RmatFn  <- function(nXbasisw, nWbasisw, nXbasisx, nWbasisx, 
                   Bvecw, Bvecx, Btens) {
  Btens <- as.vector(Btens)
  Rvec <- rep(0,nXbasisw*nXbasisx)
  inprodList <- .C("RmatFn", 
                    as.integer(nXbasisw), 
                    as.integer(nWbasisw), 
                    as.integer(nXbasisx), 
                    as.integer(nWbasisx), 
                    as.double(Bvecw), 
                    as.double(Bvecx), 
                    as.double(Btens), 
                    as.double(Rvec),
                    PACKAGE = "Data2LD")
  RmatC <- matrix(inprodList[[8]],nXbasisw,nXbasisx)
  return(RmatC)
}