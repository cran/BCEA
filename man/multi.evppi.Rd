\name{multi.evppi}
\alias{multi.evppi}

\title{
multi.evppi
}
\description{
Performs calculations for multiparameter EVPPI
}
\usage{
multi.evppi(regr.mod, inputs, he, mcsim = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{regr.mod}{
A  string specifying the form of the non parametric regression in terms of the parameters involved. "Main" effects are constructed with the format s(par), where par is one of the parameters for the model being considered. Also, it is possible to consider "interaction" terms, to account for correlation among parameters. This is done using the notation te(par1,par2), which defines the tensor function between par1 and par2
}
  \item{inputs}{
A matrix containing the simulations for all the relevant parameters. It can be constructed by calling the function CreateInputs
}
  \item{he}{
A bcea object including the health economic model
}
  \item{mcsim}{
The number of Monte Carlo simulations to be performed in order to compute an estimation of the standard error and the bias of the approximation to the EVPPI, due to the use of the non-parametric GAM regression
}
}
\references{
Baio G. (2012). Bayesian Methods in Health Economics. CRC/Chapman Hall, London
}
\author{
Gianluca Baio
}
\seealso{
\code{\link{bcea}}
}