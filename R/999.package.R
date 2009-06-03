#########################################################################/**
# @RdocPackage sfit
#
# \description{
#   @eval "packageDescription('sfit')$Description".
# }
#
# \section{Requirements}{
#   This package requires the package \pkg{R.oo}.
# }
#
# \section{Installation and updates}{
#
#   To install this package, see \url{http://www.braju.com/R/}.
# }  
#
# \section{To get started}{
#   To get started, see:
#   \enumerate{
#     \item @see "cfit.matrix" - To fit a K-dimensional simplex in an
#        N-dimensional space.
#   }
# }
#
# \section{How to cite this package}{
#   Please site [1] and [2] below.
# }
#
# \section{Wishlist}{
#   The current interface from R is very ad hoc; it dumps all data to
#   file(s), calls 'cfit' with the correct parameters via \code{pipe()}
#   (see @see "base::connections") and parses the output files from 'cfit'.
#   Ideally, we would link the cfit code to R via \code{.Call()} (see 
#   @see "base::Foreign"), but that is for the future.  The current
#   solutions has been verified to work on both Windows XP, Linux and OSX.
# }
#
# \author{
#   The algorithm and its C source code implementation is work
#   of Pratyaksha Wirapati (now at ISREC, Lousanne).
#   The R wrapper is work of Henrik Bengtsson (\url{http://www.braju.com/R/}).
# }
#
# \references{
#  [1] P. Wirapati, & T. Speed, \emph{Fitting polyhedrial cones and
#     simplices to multivariate data points}, Walter and Eliza Hall Institute 
#     of Medical Research, December 30, 2001.\cr
#  [2] P. Wirapati and T. Speed, \emph{An algorithm to fit a simplex 
#     to a set of multidimensional points}, Walter and Eliza Hall Institute 
#     of Medical Research, January 15, 2002.\cr
# }
#
# \section{License}{
#   The releases of this package is licensed under 
#   LGPL version 2.1 or newer.
#
#   The development code of the packages is under a private licence 
#   (where applicable) and patches sent to the author fall under the
#   latter license, but will be, if incorporated, released under the
#   "release" license above.
# }
#
#*/#########################################################################  
