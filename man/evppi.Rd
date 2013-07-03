\name{evppi}
\alias{evppi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expected Value of Perfect Partial Information
}
\description{
Computes the Expected Value of Perfect Partial Information (EVPPI) with respect to a given parameter. 
}
\usage{
evppi(param, he)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{
The parameter of interest. Must be in the format of a \code{n.sim} x 1 vector of simulations from its posterior.
}
  \item{he}{
A "bcea" object containing the results of the Bayesian modelling and the economic evaluation.
}
}
\details{
Computes the Expected Value of (Partial) Partial Information with respect to a given parameter, using the algorithm for fast computation of Sadatsafavi et al. (2013). Source: \url{http://webservices.core.ubc.ca/research/other/1evppi}
}
\value{
\item{evppi}{The computed values of evppi.}
}
\references{
\url{http://webservices.core.ubc.ca/research/other/1evppi}

Baio G. (2012). Bayesian Methods in Health Economics. CRC/Chapman Hall, London
}
\author{
Gianluca Baio, Andrea Berardi
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{evppi.func}}, \code{\link{evppi.plot}}
}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function (param, he) 
%{
%    f <- numeric()
%    nb <- he$U
%    evppi.graphics = FALSE
%    for (k in 1:length(he$k)) {
%        nbs <- nb[, k, ]
%        f[k] <- evppi.func(param, nbs)$evppi
%    }
%    list(evppi = f)
%  }
%}

\keyword{Health economic evaluation}
\keyword{Expected value of partial information}