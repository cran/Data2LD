\name{checkModel}
\alias{checkModel}
\title{
  Check linear differential equation specifications.
}
\description{
  A system of linear differential equations is defined by a model list of length equal 
  to the number of variables in the system.  Each variable is defined an invocation
  of function \code{make.Variable}, and then assigned to the corresponding position
  in the model list. so that each member of this list defines a single 
  linear differential equation.  
  
  Within each equation there are typically one or more 
  terms involving a coefficient function multiplying either a derivative of one of the 
  variables, or multiplying a forcing function.  
  
  This function checks that all linear 
  differential equation specifications conform to what is required.  A list of
  errors is displayed if there are problems.
  Default values are supplied where needed.  New named members are also supplied. 
  
  In addition, function \code{checkModel} also assigns position(s) in the master 
  parameter vector corresponding to coefficient functions in the system.  This
  master parameter vector is optimized by function \code{Data2LD.opt}.
  
  Finally, \code{checkModel} also displays a detailed summary of the linear
  differential system if the argument \code{summarywrd} is TRUE.
  It is easy to make mistakes in setting up lists of variable specifications,
  and it is ESSENTIAL to invoke this function prior to invoking any function using 
  this list.
}
\usage{
checkModel(XbasisList, modelList, summarywrd=TRUE)
}
\arguments{
  \item{XbasisList}{A list object containing the functional basis objects for 
  representing the variables as functional data objects.}
  \item{modelList}{A list object containing the specification of a Data2LD model. 
  Each member of this list contains a list object that defines a single linear 
  differential equation.}
  \item{summarywrd}{If TRUE, display a report on the structure of the dynamic system.}
}
\details{
  A list object containing:
  \itemize{
    \item{model list object:} {A list object with length equal to the number of variables
    in the system.  It contains variable specifications that can either be made manually, 
    or can set up by an invocation of function \code{make.Variable}.}
    \item{nparam:} {The total number of parameters in the master parameter vector.
    This master vector contains all of the parameters in the system, including any that
    are fixed as well as those that are to be optimized.}
  }
  
  If certain values are not supplied, default values are assigned.  
  
  Two new named 
  members are also supplied for each variable.  Member \code{nallXterm} is the number 
  of terms in the right side homogeneous portion of the equation, which
  involve a derivative of a variable in the system multiplied by a coefficient function.  
  Member \code{nallFterm} is the number of forcing terms in the equation.  
  These two new members are essential for function \code{Data2LD}, and,
  although \code{checkModel} is also invoked inside \code{Data2LD}, an inovking of 
  \code{checkModel} is essential prior to beginning the data analysis.
}
\value{
  A named list object with the same structure as the argument along with added 
  default values and members \code{nallXterm} and \code{nallFterm}.
}
\references{
    J. O. Ramsay and G. Hooker (2017) \emph{Dynamic Data Analysis}. Springer.
}
\seealso{
  \code{\link{make.Fterm}},
  \code{\link{make.Xterm}},
  \code{\link{make.Variable}}
}
\examples{
#  For examples of the use of this function, 
#  see these examples in the description of the function \code{make.variable}.
#
#  Example 1:  The refinery data
#
#
#  Example 2:  The X coordinate of the fda script data
#
#
#  Example 3:  The X coordinate of the fda script data
#
#
#  Example 4  Average temperature for Montreal
#
}
