\name{printModel}
\alias{printModel}
\title{
Display the definition of each variable in a system of linear differential equations.
}
\description{
The definitions of the linear differential equations in a system are containing in an
un-named list object of the same length as the number of equations.  Each member of the list
is displayed in turn.  For each equation, each of the homogenous terms in turn are
displayed, and then each forcing term displayed next.  
}
\usage{
printModel(modelList, titlestr=NULL)
}
\arguments{
  \item{modelList}{This is a list object defining each equation in a system of
  linear differential equations.}
  \item{titlestr}{A character string used as a title for the equation system.}
}
\details{
For each equation in turn, the details of the homogenous terms are displayed first,
and then the details of the forcing terms.

The homogeneous terms are those involving a product of a coefficient function 
and a derivative of a varible in the system.  

The coefficient function is defined by a set of parameters, some or all of which can be either fixed or 
estimated.  The display of the homogenous coefficient function is of: (1) the positions in the
master parameter vector of these parameters, (2) the values of these parameters, (3) TRUE/FALSE value(s)
indicating which parameters are fixed and which are estimated, and (4) the type of the coefficient function.
The coefficient function may be a functional data basis object, a functional data object, or a functional
parameter object, in which case the parmeters define a linear combination of functional basis objects.
Or, alternatively, the coefficient function may be named list whose members contain (1) \code{fd:} the definition of an R 
function, (2) \code{Dfd:} the definition of its derivative with respect to the parameters and (3) \code{more:} a 
list object containing any other information that the coefficient needs to define the values of itself and its
derivative.

The homogeneous term is also defined by the number of the variable that the coefficient function multiplies,
the order of the derivative of that function, and a fixed constant factor multiplying the product, which defaults
to one; but, for example, is often -1.

A forcing term is the product of a coefficient function and a forcing function. Each forcing function is
a fixed fixed functional data object of class \code{fd}.  The function may be single, or involve a number of
replications.  The forcing function also involves a fixed multiplying constant.
}
\references{
    J. O. Ramsay and G. Hooker (2017) \emph{Dynamic Data Analysis}. Springer.
}
