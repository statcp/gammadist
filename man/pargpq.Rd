% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pargpq.R
\name{pargpq}
\alias{pargpq}
\title{Generate the realizations of the GPQs for the parameters of the gamma distribution.}
\usage{
pargpq(x, B = 2000)
}
\arguments{
\item{x}{observations of a gamma distribution.}

\item{B}{number of realizations of the GPQs. If B is not specified, the default value is 2000.}
}
\value{
pargpq returns a dataframe which contains a list of the variables shape, scale and rate. Each variable has B realizations of the corresponding GPQ.
}
\description{
Generate the realizations of the GPQs for the parameters of the gamma distribution.
}
\details{
Assume that the data x_1, ..., x_n from a gamma distribution are available. 
The GPQs for the gamma parameters have been developed in Wang and Wu (2018). 
Because the exact distributions of the GPQs are difficult to derive, 
this pargpq function gives realizations of them using the Monte Carlo procedures.
}
\examples{
x <- rGamma(100, shape = 3, rate = 1)
pargpq(x)
pargpq(x, B=10000)
}
\references{
Wang, B. X. and Wu, F. (2018), “Inference on the Gamma Distribution,” Technometrics, 60(2), 235–244.
}
