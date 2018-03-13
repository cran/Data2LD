tensorprodfn <- function(bmat1,bmat2,bmat3,bmat4,wtvec) {
  n1 <- ncol(bmat1)
  n2 <- ncol(bmat2)
  n3 <- ncol(bmat3)
  n4 <- ncol(bmat4)
  nr <- nrow(bmat1)
  # does this call to loopJim allocate storage to tensorprod pointer?
  .C("loopJim", 
     as.double(bmat1), as.integer(n1), 
     as.double(bmat2), as.integer(n2),
     as.double(bmat3), as.integer(n3),
     as.double(bmat4), as.integer(n4),
     as.double(wtvec), as.integer(nr),
     PACKAGE="Data2LD",
     tensorprod <- double(n1*n2*n3*n4))$tensorprod  
}
