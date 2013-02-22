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
eib.plot(he,pos="topleft")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{he}{
A "bcea" object containing the results of the Bayesian modelling and the economic evaluation
}
  \item{pos}{
The position of the legend to be included in the plot. Default value is "topleft", but can be changed by the user
}
}
\value{
The function produces a plot of the Expected Incremental Benefit as a function of the discrete grid approximation of the willingness to pay parameter. The break even point (ie the point in which the EIB=0, ie when the optimal decision changes from one intervention to another) is also showed. The value k* is the discrete grid approximation of the ICER
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
      Kmax=50000            # maximum value possible for the willingness 
                            #  to pay threshold; implies that k is chosen 
                            #  in a grid from the interval (0,Kmax)
)
#
# Plots the Expected Incremental Benefit for the "bcea" object m
eib.plot(m)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Health economic evaluation}
\keyword{Expected Incremental Benefit}
