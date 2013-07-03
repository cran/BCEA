\name{bcea}
\alias{bcea}
\alias{bcea.default}
\alias{CEanalysis}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian Cost-Effectiveness Analysis
}
\description{
Cost-effectiveness analysis based on the results of a simulation model for a variable of clinical benefits (e) and of costs (c). Produces results to be post-processed to give the health economic analysis. The output is stored in an object of the class "bcea"
}
\usage{
bcea(e, c, ref = 1, interventions = NULL, Kmax = 50000, plot = FALSE)

\method{bcea}{default}(e, c, ref = 1, interventions = NULL, Kmax = 50000, plot = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{e}{
An object containing nsim simulations for the variable of clinical effectiveness for each intervention being considered. In general it is a matrix with nsim rows and nint columns
}
  \item{c}{
An object containing nsim simulations for the variable of cost for each intervention being considered. In general it is a matrix with nsim rows and nint columns
}
  \item{ref}{
Defines which interventions (columns of e or c) is considered to be the reference intervention. The default value is that the intervention associated with the first column of e or c is the reference and the one(s) associated with the other column(s) is(are) the comparators
}
  \item{interventions}{
Defines the labels to be associated with each intervention. By default and if NULL, assigns labels in the form "Intervention1", ... , "Intervention T"
}
  \item{Kmax}{
Maximum value of the willingness to pay to be considered. Default value is k=50000. The willingness to pay is then approximated on a discrete grid in the interval (0,Kmax), using steps of 100, ie [0,100,200,...,Kmax-100,Kmax].
}
  \item{plot}{
A logical value indicating whether the function should produce the summary plot or not.
  }
}
\value{
An object of the class "bcea" containing the following elements
  \item{n.sim}{Number of simulations produced by the Bayesian model}
  \item{n.comparators}{Number of interventions being analysed}
  \item{n.comparisons}{Number of possible pairwise comparisons}
  \item{delta.e}{For each possible comparison, the differential in the effectiveness measure}
  \item{delta.c}{For each possible comparison, the differential in the cost measure}
  \item{ICER}{The value of the Incremental Cost Effectiveness Ratio}
  \item{Kmax}{The maximum value assumed for the willingness to pay threshold}
  \item{k}{The vector of values for the grid approximation of the willingness to pay}
  \item{ceac}{The value for the Cost Effectiveness Acceptability Curve, as a function of the willingness to pay}
  \item{ib}{The distribution of the Incremental Benefit, for a given willingness to pay}
  \item{eib}{The value for the Expected Incremental Benefit, as a function of the willingness to pay}
  \item{kstar}{The grid approximation of the break even point (ICER)}
  \item{best}{A vector containing the numeric label of the intervention that is the most cost-effective for each value of the willingness to pay in the selected grid approximation}
  \item{U}{An array including the value of the expected utility for each simulation from the Bayesian model, for each value of the grid approximation of the willingness to pay and for each intervention being considered}
  \item{vi}{An array including the value of information for each simulation from the Bayesian model and for each value of the grid approximation of the willingness to pay}
  \item{Ustar}{An array including the maximum "known-distribution" utility for each simulation from the Bayesian model and for each value of the grid approximation of the willingness to pay}
  \item{ol}{An array including the opportunity loss for each simulation from the Bayesian model and for each value of the grid approximation of the willingness to pay}
  \item{evi}{The vector of values for the Expected Value of Information, as a function of the willingness to pay}
  \item{interventions}{A vector of labels for all the interventions considered}
  \item{ref}{The numeric index associated with the intervention used as reference in the analysis}
  \item{comp}{The numeric index(es) associated with the intervention(s) used as comparator(s) in the analysis}
  \item{step}{The step used to form the grid approximation to the willingness to pay}
}
\references{
Baio, G., Dawid, A. P. (2011). Probabilistic Sensitivity Analysis in Health Economics. Statistical Methods in Medical Research doi:10.1177/0962280211419832.

Baio G. (2012). Bayesian Methods in Health Economics. CRC/Chapman Hall, London
}
\author{
Gianluca Baio
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
      plot=TRUE          # plots the results
)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Health economic evaluation}
