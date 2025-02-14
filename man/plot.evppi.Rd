% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.evppi.R
\name{plot.evppi}
\alias{plot.evppi}
\title{Plot Expected Value of Partial Information With Respect to a
Set of Parameters}
\usage{
\method{plot}{evppi}(x, pos = c(0, 0.8), graph = c("base", "ggplot2"), col = c(1, 1), ...)
}
\arguments{
\item{x}{An object in the class \code{evppi},
obtained by the call to the function \code{\link[=evppi]{evppi()}}.}

\item{pos}{Parameter to set the position of the legend (only relevant for
multiple interventions, ie more than 2 interventions being compared).
Can be given in form
of a string \code{(bottom|top)(right|left)} for base graphics and
\code{bottom|top|left|right} for ggplot2. It can be a two-elements vector,
which specifies the relative position on the x and y axis respectively, or
alternatively it can be in form of a logical variable, with \code{FALSE}
indicating to use the default position and \code{TRUE} to place it on the
bottom of the plot.}

\item{graph}{A string used to select the graphical engine to use for
plotting. Should (partial-) match the two options \code{"base"} or
\code{"ggplot2"}. Default value is \code{"base"}.}

\item{col}{Sets the colour for the lines depicted in the graph.}

\item{...}{Arguments to be passed to methods, such as graphical parameters
(see \code{\link[=par]{par()}}).}
}
\value{
Plot with base R or \pkg{ggplot2}.
}
\description{
Plot Expected Value of Partial Information With Respect to a
Set of Parameters
}
\examples{

\dontrun{
data(Vaccine, package = "BCEA")
treats <- c("Status quo", "Vaccination")

# Run the health economic evaluation using BCEA
m <- bcea(e.pts, c.pts, ref = 2, interventions = treats)

# Compute the EVPPI for a bunch of parameters
inp <- createInputs(vaccine_mat)

# Compute the EVPPI using INLA/SPDE
if (require("INLA")) {
  x0 <- BCEA::evppi(m, c("beta.1." , "beta.2."), input = inp$mat)
  
  plot(x0, pos = c(0,1))

  x1 <- BCEA::evppi(m, c(32,48,49), input = inp$mat)
  plot(x1, pos = "topright")

  plot(x0, col = c("black", "red"), pos = "topright")
  plot(x0, col = c(2,3), pos = "bottomright")

  plot(x0, pos = c(0,1), graph = "ggplot2")
  plot(x1, pos = "top", graph = "ggplot2")

  plot(x0, col = c("black", "red"), pos = "right", graph = "ggplot2")
  plot(x0, col = c(2,3), size = c(1,2), pos = "bottom", graph = "ggplot2")

  plot(x0, graph = "ggplot2", theme = ggplot2::theme_linedraw())
}

if (FALSE)
 plot(x0, col = 3, pos = "topright")
# The vector 'col' must have the number of elements for an EVPI
# colour and each of the EVPPI parameters. Forced to black
}

}
\references{
\insertRef{Baio2013}{BCEA}
}
\seealso{
\code{\link[=bcea]{bcea()}}, \code{\link[=evppi]{evppi()}}
}
\author{
Gianluca Baio, Andrea Berardi
}
\keyword{hplot}
