% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_xxx.R
\name{compute_EIB}
\alias{compute_EIB}
\title{Compute Expected Incremental Benefit}
\usage{
compute_EIB(ib)
}
\arguments{
\item{ib}{Incremental benefit}
}
\value{
Array with dimensions (interv x k)
}
\description{
A summary measure useful to assess the potential changes in the decision
under different scenarios.
}
\details{
When considering a pairwise comparison
(e.g. in the simple case of a reference intervention \eqn{t = 1} and a comparator,
such as the status quo, \eqn{t = 0}), it is defined as the difference between the
expected utilities of the two alternatives:

\deqn{eib := \mbox{E}[u(e,c;1)] - \mbox{E}[u(e,c;0)] = \mathcal{U}^1 - \mathcal{U}^0.}

Analysis of the expected incremental benefit describes how the decision changes
for different values of the threshold. The EIB marginalises out the uncertainty,
and does not incorporate and describe explicitly the uncertainty in the outcomes.
To overcome this problem the tool of choice is the CEAC.
}
\seealso{
\code{\link[=ceac.plot]{ceac.plot()}}, \code{\link[=compute_CEAC]{compute_CEAC()}}, \code{\link[=compute_IB]{compute_IB()}}
}
