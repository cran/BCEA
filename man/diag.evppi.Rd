\name{diag.evppi}
\alias{diag.evppi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
diag.evppi
}
\description{
Performs diagnostic plots for the results of the multiparameter EVPPI
}
\usage{
diag.evppi(x, k = NULL, t = 1, diag = c("residuals", "qqplot", "bias", "se"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An \code{evppi} object containing the results of the multiparameter computations
}
  \item{k}{
The value of the willingness to pay selected for the analysis. If the diagnostics selected is either residual or qqplot, then it needs to be specified. If it is either bias or se plot, it is not relevant
}
  \item{t}{
The intervention selected for the analysis. If it is left unspecified, the default value 1 for the first intervention will be used
}
  \item{diag}{
A text string specifying the diagnostics to be performed. The default value is 'residuals'. Other possibilities are 'qqplot', 'bias' and 'se'
}
}
\references{
Baio G. (2012). Bayesian Methods in Health Economics. CRC/Chapman Hall, London
}
\author{
Gianluca Baio
}
\seealso{
\code{\link{evppi}}
}