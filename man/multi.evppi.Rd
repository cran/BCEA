\name{multi.evppi}
\alias{multi.evppi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expected Value of Perfect Partial Information for multi-parameters
}
\description{
Computes the Expected Value of Perfect Partial Information (EVPPI) with respect 
to a set of parameters using non-parametric GAM regression 
}
\usage{
multi.evppi(regr.mod,inputs,he,mcsim=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{regr.mod}{a string specifying the form of the non parametric regression in terms 
of the parameters involved. "Main" effects are constructed with the format \code{s(par)}, where 
par is one of the parameters for the model being considered. Also, it is possible to 
consider "interaction" terms, to account for correlation among parameters. This is done 
using the notation \code{te(par1,par2)}, which defines the tensor function between \code{par1} and \code{par2}
}
  \item{inputs}{
A matrix containing the simulations for all the parameters monitored by the call to 
JAGS or BUGS. The matrix should have column names matching the names of the parameters 
and the values in the vector parameter should match at least one of those values
}
  \item{he}{
A bcea object (the result of the call to the function "bcea")
}
  \item{mcsim}{the number of Monte Carlo simulations to be performed in order to compute
an estimation of the standard error and the bias of the approximation to the EVPPI, due 
to the use of the non-parametric GAM regression}
}
\details{
Computes the Expected Value of (Partial) Partial Information with respect to a set of 
parameters, using the algorithm of Strong and Oakley (2013), based on non-parametric
GAM regression. 
}
\value{
\item{evppi}{The computed values of evppi for all values of the parameter of 
willingness to pay}
\item{parameters}{The string vector with the names of the parameters for which the 
EVPPI is computed. If many parameters are considered at once and independently,
then the string is in the form \code{par1+par2}. Interaction between parameters is
indicated with the form \code{par1 : par2}}
\item{k}{The vector of values for the willingness to pay}
\item{evi}{The vector of values for the overall EVPI}
\item{SE}{The standard error of the approximation to the EVPPI}
\item{upward.bias}{The upward bias of the approximation to the EVPPI}
\item{method}{The string reporting the method selected for calculations}
}
\references{
Strong, M. and Oakley, J.E. (2013). Estimating multi-parameter partial Expected 
Value of Perfect Information from a probabilistic sensitivity analysis sample: 
a non-parametric regression approach

Baio G. (2012). Bayesian Methods in Health Economics. CRC/Chapman Hall, London
}
\author{
Gianluca Baio
}
\seealso{
\code{\link{plot.evppi}}, \code{\link{bcea}}
}
\keyword{Health economic evaluation}
\keyword{Expected value of partial information}
