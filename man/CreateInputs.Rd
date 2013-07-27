\name{CreateInputs}
\alias{CreateInputs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CreateInputs
}
\description{
Creates an object containing the matrix with the parameters simulated using the MCMC
procedure (using JAGS or BUGS) and a vector of parameters (strings) that can be 
used to perform the expected value of partial information analysis
}
\usage{
CreateInputs(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A rjags or bugs object, containing the results of a call to either "jags" (under R2jags) 
or bugs (under "R2WinBUGS" or "R2OpenBUGS")
}
}
\value{
\item{mat}{A matrix contaning all the simulations for all the monitored parameters}
\item{parameters}{A character vectors listing the names of all the monitored parameters}
}
\author{
Gianluca Baio
}
\seealso{
\code{\link{bcea}},\code{\link{evppi}},
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{R2jags}
\keyword{R2WinBUGS}
\keyword{R2OpenBUGS}% __ONLY ONE__ keyword per line
