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
ceplane.plot(he, comparison = NULL, wtp = 25000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{he}{
A "bcea" object containing the results of the Bayesian modelling and the economic evaluation
}
  \item{comparison}{
Selects the comparator, in case of more than two interventions being analysed. Default as NULL chooses the first possible. 
}
  \item{wtp}{
The value of the willingness to pay parameter
}
}
\value{
The function produces a plot of the cost-effectiveness plane. Grey dots show the simulated values for the joint distribution of the effectiveness and cost differentials. The red larger dot shows the ICER and the grey area identifies the sustainability area, ie the part of the plan for which the simulated values are below the willingness to pay threshold. The proportion of points in the sustainability area effectively represents the CEAC for a given value of the willingness to pay. If the comparators are more than 2 and no pairwise comparison is specified, all scatterplots are graphed using different colors
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
ceplane.plot(m,      # plots the Cost-Effectiveness plane
      comparison=1,  # if more than 2 interventions, selects the
                     #  pairwise comparison 
      wtp=25000      # selects the relevant willingness to pay 
                     #  (default: 25,000) 
)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Health economic evaluation}
\keyword{Cost Effectiveness Plane}% __ONLY ONE__ keyword per line
