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
contour2(he, wtp=25000, xl=NULL, yl=NULL)
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
Limits on the x-axis (default=NULL, so that R will select appropriate limits)
}
  \item{yl}{
Limits on the y-axis (default=NULL, so that R will select appropriate limits)
}
}
\value{
Plots the cost-effectiveness plane with a scatterplot of all the simulated values from the (posterior) bivariate distribution of (Delta_e,Delta_c), the differentials of effectiveness and costs; superimposes a contour of the distribution and prints the value of the ICER, together with the sustainability area 
}
\references{
Baio, G., Dawid, A. P. (2011). Probabilistic Sensitivity Analysis in Health Economics. Statistical Methods in Medical Research doi:10.1177/0962280211419832.

Baio G. (2012). Bayesian Methods in Health Economics. CRC/Chapman Hall, London
}
\author{
Gianluca Baio
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bcea}}
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
      Kmax=50000            # maximum value possible for the willingness 
                            #  to pay threshold; implies that k is chosen 
                            #  in a grid from the interval (0,Kmax)
)
#
# Plots the contour and scatterplot of the bivariate 
#   distribution of (Delta_e,Delta_c)
contour2(m,       # uses the results of the economic evalaution 
                  #  (a "bcea" object)
      wtp=25000,  # selects the willingness-to-pay threshold
      xl=NULL,    # assumes default values
      yl=NULL     # assumes default values
)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Health economic evaluation}
\keyword{ Bayesian model }% __ONLY ONE__ keyword per line
