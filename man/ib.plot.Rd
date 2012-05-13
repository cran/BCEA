\name{ib.plot}
\alias{ib.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Incremental Benefit distribution plot
}
\description{
Plots the distribution of the Incremental Benefit for a given value of the willingness to pay threshold, k
}
\usage{
ib.plot(he, comparison = NULL, wtp = 25000, bw = nbw, n = 512, xlim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{he}{
A "bcea" object containing the results of the Bayesian modelling and the economic evaluation
}
  \item{comparison}{
In the case of multiple interventions, specifies the one to be used in comparison with the reference. Default value of NULL forces R to consider the first non-reference intervention as the comparator 
}
  \item{wtp}{
The value of the willingness to pay threshold. Default value at 25000
}
  \item{bw}{
Identifies the smoothing bandwith used to construct the kernel estimation of the IB density.
}
  \item{n}{
The number of equally spaced points at which the density is to be estimated
}
  \item{xlim}{
The limits of the plot on the x-axis 
}
}
\value{
The function produces a plot of the distribution of the Incremental Benefit for a given value of the willingness to pay parameter. The dashed area indicates the positive part of the distribution (ie when the reference is more cost-effective than the comparator) 
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
# Plots the distribution of the Incremental Benefit
ib.plot(m,        # uses the results of the economic evalaution 
                  #  (a "bcea" object)
    comparison=1, # if more than 2 interventions, selects the 
                  #  pairwise comparison 
    wtp=25000     # selects the relevant willingness 
                  #  to pay (default: 25,000) 
)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Health economic evaluation}
