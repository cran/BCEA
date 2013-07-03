\name{evppi.func}
\alias{evppi.func}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expected Value of Perfect Partial Information
}
\description{
Computes the Expected Value of Perfect Partial Information (EVPPI) with respect to a given parameter. 
}
\usage{
evppi.func(param, nbs, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{
The parameter of interest. Must be in the format of a \code{n.sim} x 1 vector of simulations from its posterior.
}
  \item{nbs}{
A matrix of net benefits for the interventions being compared.
  }
  \item{\dots}{
Externally set segmentation points for the algorithm.
}
}
\details{
Computes the Expected Value of (Partial) Partial Information with respect to a given parameter, using the algorithm for fast computation of Sadatsafavi et al. (2013). Source: \url{http://webservices.core.ubc.ca/research/other/1evppi}
}
\value{
\item{evppi}{The computed values of evppi.}
\item{segPoints}{The segmentation points used in the algorithm.}
}
\references{
\url{http://webservices.core.ubc.ca/research/other/1evppi}

Baio G. (2012). Bayesian Methods in Health Economics. CRC/Chapman Hall, London
}
\author{
Gianluca Baio
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{evppi}}, \code{\link{evppi.plot}}
}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%
%}

\keyword{Health economic evaluation}
\keyword{Expected value of partial information}