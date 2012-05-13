\name{summary.bcea}
\alias{summary.bcea}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary method for objects in the class "bcea"
}
\description{
Produces a table printout with some summary results of the health economic evaluation
}
\usage{
\method{summary}{bcea}(object, wtp = 25000, ...)
%%summary(object, wtp = 25000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A "bcea" object containing the results of the Bayesian modelling and the economic evaluation
}
  \item{wtp}{
The value of the willingness to pay threshold to be used in the summary table
}
\item{...}{
Additional arguments affecting the summary produced
}
}

\value{
Prints a summary table with some information on the health economic output and synthetic information on the economic measures (EIB, CEAC, EVPI)
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
summary(m,      # uses the results of the economic evalaution 
                #  (a "bcea" object)
      wtp=25000	# selects the particular value for k 
)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Health economic evaluation}% __ONLY ONE__ keyword per line
