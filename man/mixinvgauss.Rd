\name{Finite mixture of inverse Gaussian Distributions}
\alias{Finite mixture of inverse Gaussian Distributions}
  
\alias{pmixinvgauss}
\alias{dmixinvgauss}
\alias{mixinvgaussHazard}

\title{
Finite mixture of inverse Gaussian Distribution
}
\description{
Density (PDF), distribution function (CDF), and hazard function for Finite mixture of inverse Gaussian Distributions.
}
\usage{
dmixinvgauss(x, theta = .2, lambda = .1, gamma = .05, forceExpectation = F)
pmixinvgauss(q, theta = .2, lambda = .1, gamma = .05, forceExpectation = F)
mixinvgaussHazard(x, theta = .2, lambda = .1, gamma = .05, forceExpectation = F)
}
\arguments{
  \item{x, q}{
vector of quantiles.
}
  \item{theta, lambda, gamma}{
parameters, see 'Details'.
}
  \item{forceExpectation}{
logical; if \code{TRUE}, the expectation of the distribution is forced to be 1..
}
}
\details{
The finite mixture of inverse Gaussian distributions was used by Gomes-Deniz and Perez-Rodrigues (201X) for ACD-models. Its PDF is:


\deqn{f(x) = \frac{\gamma + x}{\gamma + \theta} \sqrt{\frac{\lambda}{2 \pi x^3}} \exp \left[ - \frac{\lambda(x-\theta)^2}{2 x \theta^2}\right].}


If \code{forceExpectation = TRUE} the distribution is transformed by dividing the random variable with its expectation and using the change of variable function.

}

\references{

Gomez-Deniz Perez-Rodriguez (201X) 
\emph{Non-exponential mixtures, non-monotonic financial hazard functions and the autoregressive conditional duration model.} Working paper. Retrieved June 16, 2015, from http://dea.uib.es/digitalAssets/254/254084_perez.pdf.

}