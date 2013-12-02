\name{struct.psa}
\alias{struct.psa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
struct.psa
}
\description{
Computes the weights to be associated with a set of competing models in order to 
perform structural PSA model is a list containing the output from either R2jags 
or R2OpenBUGS/R2WinBUGS for all the models that need to be combined in the model 
average effect is a list containing the measure of effectiveness computed from the 
various models (one matrix with n.sim x n.ints simulations for each model) cost 
is a list containing the measure of costs computed from the various models (one 
matrix with n.sim x n.ints simulations for each model)
}
\usage{
struct.psa(models,effect,cost,ref=1,interventions=NULL,Kmax=50000,plot=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{A list containing the output from either R2jags or 
R2OpenBUGS/R2WinBUGS for all the models that need to be combined in the 
model average
}
  \item{effect}{A list containing the measure of effectiveness computed 
from the various models (one matrix with n.sim x n.ints simulations for 
each model)
}
  \item{cost}{A list containing the measure of costs computed from the 
various models (one matrix with n.sim x n.ints simulations for each model)
}
  \item{ref}{Defines which interventions (columns of \code{e} or \code{c}) 
is considered to be the reference intervention. The default value is that the 
intervention associated with the first column of \code{e} or \code{c} is the 
reference and the one(s) associated with the other column(s) is(are) the 
comparators
}
  \item{interventions}{
Defines the labels to be associated with each intervention. By default and 
if NULL, assigns labels in the form "Intervention1", ... , "Intervention T"
}
  \item{Kmax}{
Maximum value of the willingness to pay to be considered. Default value is 
k=50000. The willingness to pay is then approximated on a discrete grid in 
the interval (0,Kmax), using steps of 100, ie [0,100,200,...,Kmax-100,Kmax].
}
\item{plot}{
A logical value indicating whether the function should produce the summary 
plot or not.
}
}
\value{
\item{he}{A "bcea" object containing the economic analysis for the 
average model, obtained by suitably weighting the suggested models}
\item{w}{A vector of weights, each one associated with one of the 
model being combined}
\item{DIC}{A vector of values for the models DIC - these are used 
to compute the weights by suitably re-scaling (cfr Baio, 2012)}
}
\references{
Baio G. (2012). Bayesian Methods in Health Economics. CRC/Chapman Hall, London

Jackson, C. et al (2011). A framework for addressing structural uncertainty in 
decision models. Medical Decision Making, 31, 662-674

Jackson, C. et al (2010). Accounting for uncertainty in health economic decision 
models by using model averaging. JRSS A, 172(2), 383-404 
}
\author{
Gianluca Baio
}
\seealso{
\code{\link{bcea}},
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{R2jags}
\keyword{R2WinBUGS}
\keyword{R2OpenBUGS}% __ONLY ONE__ keyword per line
