\name{evppi}
\alias{evppi}
\alias{evppi.default}
\title{
Expected Value of Perfect Partial Information (EVPPI) for selected parameters
}
\description{
Calculates the Expected Value of Perfect Partial Information (EVPPI) for subsets 
of parameters. Uses GAM non-parameteric regression for single parameter EVPPI and 
the SPDE-INLA method for larger parameter subsets.
}
\usage{
evppi(parameter,input,he,N=NA,plot=F,residuals=F,h.value=0.00005)

\method{evppi}{default}(parameter,input,he,N=NA,plot=F,residuals=F,h.value=0.00005)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{parameter}{
A vector of parameters for which the EVPPI should be calculated. This can be given as a 
string (or vector of strings) of names or a numeric vector, corresponding to the column 
numbers of important parameters.
}
\item{input}{
A matrix containing the simulations for all the parameters monitored by the call to 
JAGS or BUGS. The matrix should have column names matching the names of the parameters 
and the values in the vector parameter should match at least one of those values.
}
\item{he}{
A \code{bcea} object (the result of the call to the function \code{\link{bcea}}).
}
\item{N}{
The number of PSA simulations used to calculate the EVPPI. The default uses all the 
available samples.
}
\item{plot}{
A logical value indicating whether the triangular mesh for SPDE-INLA should be plotted.
Default set to F.
}
\item{residuals}{
A logical value indicating whether the fitted values for the SPDE-INLA method should be 
outputted. Default set to F.
}
\item{h.value}{
A numerical value that indicates the step size for the INLA algorithm. If there is an 
error in the INLA algorithm, reduce step size.
}
}
\details{
The single parameter EVPPI has been calculated using the non-parametric GAM regression 
developed by Strong et al. (2014). The multi-parameter EVPPI is calculated using the 
SPDE-INLA regression method for Gaussian Process regression developed by Heath et al. 
(2015) 
}
\value{ 
\item{evppi}{The computed values of evppi for all values of the parameter of 
willingness to pay}
\item{index}{A numerical vector with the index associated with the parameters for which 
the EVPPI was calculated}
\item{k}{The vector of values for the willingness to pay}
\item{evi}{The vector of values for the overall EVPPI}
\item{fitted.costs}{The fitted values for the costs}
\item{fitted.effects}{The fitted values for the effects}
\item{parameters}{A single string containing the names of the parameters for which the EVPPI was calculated, used for plotting the EVPPI}
\item{time}{Computational time (in seconds)}
}
\references{
Strong M., Oakley J. and Brennan A. (2014), Estimating multi-parameter partial 
Expected Value of Perfect Information from a probabilistic sensitivity analysis sample: 
a non-parametric regression approach, Medical Decision Making.

Baio G. (2012). Bayesian Methods in Health Economics. CRC/Chapman Hall, London

Heath A., Manolopoulou I., Baio G. (2015). Efficient High-Dimensional Gaussian Process 
Regression to calculate the Expected Value of Partial Perfect Information in Health 
Economic Evaluations. http://arxiv.org/pdf/1504.05436v1.pdf
}
\author{
Anna Heath, Gianluca Baio
}
\examples{
# See Baio G., Dawid A.P. (2011) for a detailed description of the 
# Bayesian model and economic problem
#
# Load the processed results of the MCMC simulation model
data(Vaccine)
# 
# Runs the health economic evaluation using BCEA
m <- bcea(e,c,ref=2,interventions=treats)
#
\donttest{
# Computes the EVPPI for a bunch of parameters
inp <- CreateInputs(vaccine)
ev <- evppi(parameter=c(38:44),input=inp$mat,he=m)
}
}
\seealso{
\code{\link{plot.evppi}}, \code{\link{bcea}}
}
\keyword{Health economic evaluation}
\keyword{Expected value of partial information}
