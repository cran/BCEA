\name{ceplane.plot}
\alias{ceplane.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cost-effectiveness plane plot
}
\description{
Produces a scatter plot of the cost-effectiveness plane, together with the sustainability area, as a function of the selected willingness to pay threshold
}
\usage{
ceplane.plot(he, comparison = NULL, wtp = 25000, pos=c(1,1), 
size=NULL, graph=c("base","ggplot2"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{he}{
A "bcea" object containing the results of the Bayesian modelling and the economic evaluation
}
  \item{comparison}{
Selects the comparator, in case of more than two interventions being analysed. Default as NULL plots all the comparisons together. With \code{graph="ggplot2"}, any subset of the possible comparisons can be selected (e.g., \code{comparison=c(1,3)}), while only one can be chosen using base graphics (the default).
}
  \item{wtp}{
The value of the willingness to pay parameter
}
  \item{pos}{
Parameter to set the position of the legend; for a single comparison plot, the ICER legend position. Can be given in form of a string \code{(bottom|top)(right|left)} for base graphics and \code{bottom|top|left|right} for ggplot2. It can be a two-elements vector, which specifies the relative position on the x and y axis respectively, or alternatively it can be in form of a logical variable, with \code{FALSE} indicating to use the default position and \code{TRUE} to place it on the bottom of the plot. Default value is \code{c(1,1)}, that is the topright corner inside the plot area.
  }
  \item{size}{
Value (in millimetres) of the size of the willingness to pay label. Used only if \code{graph="ggplot2"}, otherwise is ignored with a message.
  }
  \item{graph}{
A string used to select the graphical engine to use for plotting. Should (partial-)match the two options \code{"base"} or \code{"ggplot2"}. Default value is \code{"base"}.
  }
  \item{...}{
If \code{graph="ggplot2"} and a named theme object is supplied, it will be added to the ggplot object. Also, there are some other options: \code{ICER.size=} can be used to resize the red dot representing the ICER, while \code{label.pos=FALSE} uses a different position for the willingness to pay label.
  }
}
\value{
\item{ceplane}{
A ggplot object containing the plot. Returned only if \code{graph="ggplot2"}.
}
The function produces a plot of the cost-effectiveness plane. Grey dots show the simulated values for the joint distribution of the effectiveness and cost differentials. The red larger dot shows the ICER and the grey area identifies the sustainability area, ie the part of the plan for which the simulated values are below the willingness to pay threshold. The proportion of points in the sustainability area effectively represents the CEAC for a given value of the willingness to pay. If the comparators are more than 2 and no pairwise comparison is specified, all scatterplots are graphed using different colors.
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
\code{\link{bcea}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Health economic evaluation}
\keyword{Cost Effectiveness Plane}% __ONLY ONE__ keyword per line
