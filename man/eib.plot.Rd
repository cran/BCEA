\name{eib.plot}
\alias{eib.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Expected Incremental Benefit plot
}
\description{
Plots the Expected Incremental Benefit as a function of the willingness to pay
}
\usage{
eib.plot(he, comparison=NULL, pos=c(1,0), size=NULL, graph=c("base","ggplot2"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{he}{
A "bcea" object containing the results of the Bayesian modelling and the economic evaluation
}
  \item{comparison}{
Will be ignored with a message using \code{graph="base"}. If \code{graph="ggplot2"} selects the comparators, in case of more than two interventions being analysed. Default as \code{NULL} chooses all the possible comparisons. Any subset of the possible comparisons is accepted.
  }
  \item{pos}{
Parameter to set the position of the legend; for a single comparison plot, the ICER legend position. Can be given in form of a string \code{(bottom|top)(right|left)} for base graphics and \code{bottom|top|left|right} for ggplot2. It can be a two-elements vector, which specifies the relative position on the x and y axis respectively, or alternatively it can be in form of a logical variable, with \code{FALSE} indicating to use the default position and \code{TRUE} to place it on the bottom of the plot. Default value is \code{c(1,0)}, that is the bottomright corner inside the plot area.
}
  \item{size}{
Value (in millimetres) of the size of the willingness to pay label. Used only if \code{graph="ggplot2"}, otherwise is ignored with a message.
  }
  \item{graph}{
A string used to select the graphical engine to use for plotting. Should (partial-)match the two options \code{"base"} or \code{"ggplot2"}. Default value is \code{"base"}.
  }
  \item{...}{
If \code{graph="ggplot2"} and a named theme object is supplied, it will be added to the ggplot object.
  }
}
\value{
\item{eib}{
A ggplot object containing the requested plot. Returned only if \code{graph="ggplot2"}.
}
The function produces a plot of the Expected Incremental Benefit as a function of the discrete grid approximation of the willingness to pay parameter. The break even point (ie the point in which the EIB=0, ie when the optimal decision changes from one intervention to another) is also showed. The value k* is the discrete grid approximation of the ICER.
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
\code{\link{ib.plot}},
\code{\link{ceplane.plot}}
}
\examples{
# See Baio G., Dawid A.P. (2011) for a detailed description of the 
# Bayesian model and economic problem
#
# Load the processed results of the MCMC simulation model
data(Vaccine)
# 
# Runs the health economic evaluation using BCEA
m <- bcea(e=e,c=c,          # defines the variables of 
                            #  effectiveness and cost
      ref=2,                # selects the 2nd row of (e,c) 
                            #  as containing the reference intervention
      interventions=treats, # defines the labels to be associated 
                            #  with each intervention
      Kmax=50000,           # maximum value possible for the willingness 
                            #  to pay threshold; implies that k is chosen 
                            #  in a grid from the interval (0,Kmax)
      plot=FALSE            # inhibits graphical output
)

### Plots the Expected Incremental Benefit for the "bcea" object m
eib.plot(m)

### Using ggplot2:
eib.plot(m,graph="ggplot2")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Health economic evaluation}
\keyword{Expected Incremental Benefit}
