\name{contour.bcea}
\alias{contour.bcea}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Contour method for objects in the class "bcea"
}
\description{
Produces a scatterplot of the cost-effectiveness plane, with a contour-plot of the bivariate density of the differentials of cost (y-axis) and effectiveness (x-axis)
}
\usage{
\method{contour}{bcea}(x, comparison = 1, scale = 0.5, nlevels = 4, levels = NULL, 
pos = c(1,0), graph=c("base","ggplot2"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A "bcea" object containing the results of the Bayesian modelling and the economic evaluation
}
  \item{comparison}{
In case of more than 2 interventions being analysed, selects which plot should be made. By default the first comparison among the possible ones will be plotted. If \code{graph="ggplot2"} any subset of the possible comparisons can be selected, and \code{comparison=NULL} will yield a plot of all the possible comparisons together.
}
  \item{scale}{
Scales the plot as a function of the observed standard deviation
}
  \item{levels}{
Numeric vector of levels at which to draw contour lines. Will be ignored using \code{graph="ggplot2"}.
}
  \item{nlevels}{
Number of levels to be plotted in the contour
}
  \item{pos}{
Parameter to set the position of the legend. Can be given in form of a string \code{(bottom|top)(right|left)} for base graphics and \code{bottom|top|left|right} for ggplot2. It can be a two-elements vector, which specifies the relative position on the x and y axis respectively, or alternatively it can be in form of a logical variable, with \code{FALSE} indicating to use the default position and \code{TRUE} to place the legend on the bottom of the plot. Default value is \code{c(1,0)}, that is the bottomright corner inside the plot area.
  }
  \item{graph}{
  A string used to select the graphical engine to use for plotting. Should (partial-)match the two options \code{"base"} or \code{"ggplot2"}. Default value is \code{"base"}.
  }
\item{...}{
Additional arguments to 'plot.window', 'title', 'Axis' and 'box', typically graphical parameters such as 'cex.axis'. Will be ignored if \code{graph="ggplot2"}.
}
}
\value{
\item{ceplane}{
A ggplot object containing the plot. Returned only if \code{graph="ggplot2"}.
}
Plots the cost-effectiveness plane with a scatterplot of all the simulated values from the (posterior) bivariate distribution of (Delta_e,Delta_c), the differentials of effectiveness and costs; superimposes a contour of the distribution and prints the estimated value of the probability of each quadrant (combination of positive/negative values for both Delta_e and Delta_c) 
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
\code{\link{contour2}}
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
#
# Plots the contour and scatterplot of the bivariate 
# distribution of (Delta_e,Delta_c)
contour(m,          # uses the results of the economic evalaution 
                    #  (a "bcea" object)
      comparison=1, # if more than 2 interventions, selects the 
                    #  pairwise comparison 
      nlevels=4,    # selects the number of levels to be 
                    #  plotted (default=4)
      levels=NULL,  # specifies the actual levels to be plotted 
                    #  (default=NULL, so that R will decide)
      scale=0.5,    # scales the bandwiths for both x- and 
                    #  y-axis (default=0.5)
      graph="base"  # uses base graphics to produce the plot
)

# Using gggplot2:
contour(m,            # uses the results of the economic evalaution 
                      #  (a "bcea" object)
      comparison=1,   # if more than 2 interventions, selects the 
                      #  pairwise comparison 
      nlevels=4,      # selects the number of levels to be 
                      #  plotted (default=4)
      scale=0.5,      # scales the bandwiths for both x- and 
                      #  y-axis (default=0.5)
      graph="ggplot2" # uses ggplot2 to produce the plot
)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Health economic evaluation}
\keyword{ Bayesian model }% __ONLY ONE__ keyword per line
