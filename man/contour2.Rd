\name{contour2}
\alias{contour2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Specialised contour plot for objects in the class "bcea"
}
\description{
Produces a scatterplot of the cost-effectiveness plane, with a contour-plot of the bivariate density of the differentials of cost (y-axis) and effectiveness (x-axis). Also adds the sustainability area (ie below the selected value of the willingness-to-pay threshold).
}
\usage{
contour2(he, wtp=25000, xl=NULL, yl=NULL, comparison=NULL, graph=c("base","ggplot2"),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{he}{
A "bcea" object containing the results of the Bayesian modelling and the economic evaluation
}
  \item{wtp}{
The selected value of the willingness-to-pay. Default is 25000
}
  \item{xl}{
Limits on the x-axis (default=\code{NULL}, so that R will select appropriate limits)
}
  \item{yl}{
Limits on the y-axis (default=\code{NULL}, so that R will select appropriate limits)
}
  \item{comparison}{
The comparison being plotted. Default to \code{NULL} chooses the first comparison with \code{graph="base"}, while if \code{graph="ggplot2"} it will choose all the possible comparisons. Any subset of the possible comparisons can be selected (e.g., \code{comparison=c(1,3)}). Not relevant if head-to-head comparison.
}
  \item{graph}{
   A string used to select the graphical engine to use for plotting. Should (partial-)match the two options \code{"base"} or \code{"ggplot2"}. Default value is \code{"base"}.
  }
  \item{...}{
Arguments to be passed to \code{\link{ceplane.plot}}. See the relative manual page for more details.
  }
}
\value{
\item{contour}{
A ggplot item containing the requested plot. Returned only if \code{graph="ggplot2"}.
}
Plots the cost-effectiveness plane with a scatterplot of all the simulated values from the (posterior) bivariate distribution of (Delta_e,Delta_c), the differentials of effectiveness and costs; superimposes a contour of the distribution and prints the value of the ICER, together with the sustainability area 
}
\references{
Baio, G., Dawid, A. P. (2011). Probabilistic Sensitivity Analysis in Health Economics. Statistical Methods in Medical Research doi:10.1177/0962280211419832.

Baio G. (2012). Bayesian Methods in Health Economics. CRC/Chapman Hall, London
}
\author{
Gianluca Baio, Andrea Berardi
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bcea}},
\code{\link{ceplane.plot}},
\code{\link{contour.bcea}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Health economic evaluation}
\keyword{ Bayesian model }% __ONLY ONE__ keyword per line
