\name{diag.evppi}
\alias{diag.evppi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diagnostic functions for the expected Value of Perfect Partial Information 
for multi-parameters
}
\description{
Computes a set of diagnostics to evaluate the performance of the non-parametric 
GAM regression used to compute the Expected Value of Perfect Partial Information 
(EVPPI) with respect to a set of parameters
}
\usage{
diag.evppi(x,k=NULL,t=1,diag=c("residuals","qqplot","bias","se"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{x}{an evppi objects obtained as a result of the call to the function 
\code{multi.evppi}}
\item{k}{is the value of the willingness to pay selected for the analysis if the 
diagnostics selected is either \code{residual} or \code{qqplot}, then it needs to be specified. 
If it is either \code{bias} or \code{se} plot, it is not relevant}
\item{t}{is the intervention selected for the analysis. If it is left unspecified, 
the default value 1 for the first intervention will be used}
\item{diag}{is a text string specifying the diagnostics to be performed. The default 
value is '\code{residuals}'. Other possibilities are '\code{qqplot}', '\code{bias}' and '\code{se}'}
}
\details{
Produces some graphical diagnostics for the Expected Value of (Partial) Partial 
Information with respect to a set of parameters, obtained using the algorithm of 
Strong and Oakley (2013), based on non-parametric GAM regression. 
}
\references{
Strong, M. and Oakley, J.E. (2013). Estimating multi-parameter partial Expected 
Value of Perfect Information from a probabilistic sensitivity analysis sample: 
a non-parametric regression approach

Baio G. (2012). Bayesian Methods in Health Economics. CRC/Chapman Hall, London
}
\author{
Gianluca Baio
}
\seealso{
\code{\link{multi.evppi}}, \code{\link{bcea}}
}
\keyword{Health economic evaluation}
\keyword{Expected value of partial information}
