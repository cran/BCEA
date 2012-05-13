\name{ceac.plot}
\alias{ceac.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cost Effectiveness Acceptability Curve (CEAC) plot
}
\description{
Produces a plot of the CEAC against the willingness to pay threshold
}
\usage{
ceac.plot(he)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{he}{
A "bcea" object containing the results of the Bayesian modelling and the economic evaluation
}
}
\value{
The function produces a plot of the cost-effectiveness acceptability curve against the discrete grid of possible values for the willingness to pay parameter. Values of the CEAC closer to 1 indicate that uncertainty in the cost-effectiveness of the reference intervention is very low. Similarly, values of the CEAC closer to 0 indicate that uncertainty in the cost-effectiveness of the comparator is very low
}
\references{
Baio, G., Dawid, A. P. (2011). Probabilistic Sensitivity Analysis in Health Economics. Statistical Methods in Medical Research doi:10.1177/0962280211419832.

Baio G. (2012). Bayesian Methods in Health Economics. CRC/Chapman Hall, London
}
\author{
Gianluca Baio
}

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
# Produces a plot of the CEAC against a grid of values for the 
# willingness to pay threshold
ceac.plot(m)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Health economic evaluation}
\keyword{Cost Effectiveness Acceptability Curve}% __ONLY ONE__ keyword per line
