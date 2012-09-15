\name{ceaf.plot}
\alias{ceaf.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots the Cost-effectiveness acceptability frontier 
}
\description{
Plots the cost-effectiveness acceptability frontier
}
\usage{
ceaf.plot(mce)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mce}{
mce is the output of the call to the function multi.ce}
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
mce <- multi.ce(m           # uses the results of the economic analysis 
)
ceaf.plot(mce)              # plots the CEAF
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Health economic evaluation}
\keyword{Multiple comparison}
