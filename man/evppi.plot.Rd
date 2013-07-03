\name{evppi.plot}
\alias{evppi.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expected value of Perfect Information Plot
}
\description{
Plots the Expected Value of Perfect Partial Information (EVPPI)
}
\usage{
evppi.plot(he, evppi=NULL, param = NULL, graph = c("base", "ggplot2"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{he}{
A "bcea" object containing the results of the Bayesian modelling and the economic evaluation.
}
  \item{evppi}{
The values of evppi. Must be in the format of a vector of the same length of \code{he$k}, that is the willingness to pay grid of values. If both \code{evppi} and \code{param} are supplied, the second will be ignored with a warning.
  }
  \item{param}{
The parameter of interest. Must be in the format of a \code{n.sim} x 1 vector of simulations from its posterior. Will be used to compute the values (via \code{\link{evppi}}) if they are not contained in \code{evppi}.
}
  \item{graph}{
A string used to select the graphical engine to use for plotting. Should (partial-)match the two options \code{"base"} or \code{"ggplot2"}. Default value is \code{"base"}.
}
}
%\details{
%%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\item{evppi}{
A ggplot object containing the requested plot. Returned only if \code{graph="ggplot2"}.
}
}
\references{
Baio G. (2012). Bayesian Methods in Health Economics. CRC/Chapman Hall, London
}
\author{
Andrea Berardi
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{evppi}}, \code{\link{evppi.func}}
}
%\examples{
%}

\keyword{Health economic evaluation}
\keyword{Expected value of partial information}