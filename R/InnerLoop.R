InnerLoop <- function( in1, in2, in3, in4 , wvec){

output_args <- .Call("loopJuan", in1, in2, in3, in4, wvec, PACKAGE="Data2LD")

return(output_args)

}
 
