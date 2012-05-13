\name{plot.bcea}
\alias{plot.bcea}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary plot of the health economic analysis
}
\description{
Plots in a single graph the Cost-Effectiveness plane, the Expected Incremental Benefit, the CEAC and the EVPI
}
\usage{
\method{plot}{bcea}(x,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{
A "bcea" object containing the results of the Bayesian modelling and the economic evaluation
}
\item{...}{
Arguments to be passed to methods, such as graphical parameters (see 'par')
}
}

\value{
The function produces a plot with four graphical summaries of the health economic evaluation
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
\code{\link{ceplane.plot}},
\code{\link{ceac.plot}},
\code{\link{evi.plot}}
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
# Plots the summary plots for the "bcea" object m
plot(m)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Health economic evaluation}% __ONLY ONE__ keyword per line
