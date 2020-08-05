\name{Atensorfn}
\alias{Atensorfn}
\title{
  Compute the four-way tensors corresponding to pairs of forcing terms in systems 
  of linear differential equations.
}
\description{
  Forcing terms consist of the product of a coefficient function that must be estimated 
  and a known forcing function that does not.  When both of these functions in a forcing 
  term are defined by B-splines, the product involves an inner product of two B-spline 
  basis systems.  When a product of two forcing terms are required, as is usual in the
  the use of the Data2LD package, a great improvement in efficiency of computation can be 
  acheived by an initial computation of the four-way array or tensor resulting by 
  taking the inner products of all possible quadruples of the B-spline basis functions.  
  Memoization is the process of storing these tensors in memory so that they do
  not need to be re-computed each time the Data2LD.R function is called.  
  Memoization is taken care of automatically in the code using the R.cache package, 
  and is activated the first time a new \code{modelList} object is encountered.
  Normally the user does not have to worry about the memorization procedure.  It is
  possible, however, to manually re-activate the memoization.  
  However, users may also want to construct these four-way tensors
  manually for debugging and other purposes, and this function is made available for 
  this reason.
}
\usage{
  Atensorfn(modelList)
}
\arguments{
  \item{modelList}{A list object containing the specification of a Data2LD model. 
  Each member of this list contains a list object that defines a single linear 
  differential equation.}
}
\details{
  Variable specifications can be set by an invocation of function
  \code{make.Variable} for each linear differential equation in the system.
  Functions \code{Atensorfn}, \code{BAtensorfn} and \code{Btensorfn} will not
  normally be required by users since they are invoked automatically with the
  function \code{checkModel}, which is required to be invoked before the
  analysis of the data.
}

\value{
  A list object of length equal to the number of variables in the system.  
  Each of the members of this list is a two-dimensional list object, and the members 
  of this list are the four-way tensors set up as vectors for each of the possible pairs 
  of forcing terms.  All levels of the this list structure are designed to be accessed 
  numerically by a call like \code{myAtensor[[ivar]][[ntermj]][[ntermk]]}.
}
\references{
    J. O. Ramsay and G. Hooker (2017) \emph{Dynamic Data Analysis}. Springer.
}
