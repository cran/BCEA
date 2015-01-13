\name{evppi}
\alias{evppi}
\alias{evppi.default}
\title{
Expected Value of Perfect Partial Information (EVPPI) for selected parameters
}
\description{
Computes the Expected Value of Perfect Partial Information (EVPPI) with respect 
to a given parameter. 
}
\usage{
evppi(parameters,inputs,he,n.blocks=NULL,n.seps=NULL,method=c("so","sal"))

\method{evppi}{default}(parameters,inputs,he,n.blocks=NULL,n.seps=NULL,method=c("so","sal"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parameters}{
A string (or vector of strings) containing the names of the parameters for which the EVPPI 
should be computed. If more than one parameter is specified, the EVPPI will be computed for 
each separately.
}
  \item{inputs}{
A matrix containing the simulations for all the parameters monitored by the call to 
JAGS or BUGS. The matrix should have column names matching the names of the parameters 
and the values in the vector parameter should match at least one of those values.
}
  \item{he}{
A \code{bcea} object (the result of the call to the function \code{\link{bcea}}).
}
  \item{n.blocks}{
The number of blocks in which the matrix of inputs should be decomposed to perform the 
calculations. The number of blocks must be such that the ratio of the number of MCMC 
simulations to the number of blocks is an integer number. The choice of this parameter 
is potentially quite sensitive, since too large values imply an upward bias in the EVPPI. 
If n.blocks = the number of simulations, then the computed EVPPI actually coincides with 
the overall EVPI. If n.blocks = 1, then the estimated EVPPI is identically 0. This 
parameter is used only if \code{method = "so"}.}
  \item{n.seps}{
The number of points defining the separation in the matrix of inputs. The available 
number of separation points is 1 (default) or 2. This parameter is only used if 
\code{method = "sal"}.}
  \item{method}{
The method used to perform the EVPPI calculation. Options are Strong and Oakley (\code{"so"}) or 
Sadatsafavi et al (\code{"sal"}).
}
}
\details{
Computes the Expected Value of (Partial) Partial Information with respect to a given 
parameter, using the algorithm for fast computation of Strong and Oakley (2013) or of 
Sadatsafavi et al. (2013). 
}
\value{
\item{evppi}{The computed values of evppi for all values of the parameter of 
willingness to pay}
\item{parameters}{The string vector with the names of the parameters for which the 
EVPPI is computed}
\item{k}{The vector of values for the willingness to pay}
\item{evi}{The vector of values for the overall EVPI}
\item{method}{The string reporting the method selected for calculations}
}
\references{
Strong, M. and Oakley, J.E. (2013). An efficient method for computing partial 
         expected value of perfect information for correlated inputs

Sadatsafavi et al. (2013). Need for Speed: An efficient algorithm for calculation
               of single-parameter expected value of partial perfect information

Baio G. (2012). Bayesian Methods in Health Economics. CRC/Chapman Hall, London
}
\author{
Gianluca Baio, Andrea Berardi
}
\seealso{
\code{\link{plot.evppi}}, \code{\link{bcea}}
}
\keyword{Health economic evaluation}
\keyword{Expected value of partial information}
