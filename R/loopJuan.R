loopJuan <- function( bmat1, bmat2, bmat3, bmat4 , wvec) {
  nrow    = nrow(bmat1)
  nbasis1 = ncol(bmat1)
  nbasis2 = ncol(bmat2)
  nbasis3 = ncol(bmat3)
  nbasis4 = ncol(bmat4)
  bmat1 = matrix(bmat1,nbasis1*nrow,1)
  bmat2 = matrix(bmat2,nbasis2*nrow,1)
  bmat3 = matrix(bmat3,nbasis3*nrow,1)
  bmat4 = matrix(bmat4,nbasis4*nrow,1)
  Btensor = .Call("loopJuanC", as.double(bmat1), as.double(bmat2), 
                               as.double(bmat3), as.double(bmat4), 
                               as.double(wvec))
  return(Btensor)
  
}

