\name{acdFit-methods}
\alias{acdFit-methods}
\alias{print.acdFit}
\alias{coef.acdFit}
\alias{residuals.acdFit}

\title{
Methods for class acdFit
}
\description{
\code{residuals.acdFit()} returns the residuals and \code{coef.acdFit()} returns the coefficients of a fitted ACD model of class 'acdFit', while \code{print.acdFit()} prints the essential information.
}
\usage{
\method{residuals}{acdFit}(object, \dots)
\method{coef}{acdFit}(object, returnCoef = "all", \dots)
\method{print}{acdFit}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
the fitted ACD model of class 'acdFit' (as returned by the function \code{acdFit}).
}
  \item{x}{
same as \code{object}, ie. an object of class 'acdFit'.
}
 \item{returnCoef}{
on of "all", "distribution", or "model". Specifies whether all estimated parameters should be returned or only the distribution parameters or the model (for the conditional mean duration) parameters.
}
  \item{\dots}{
  additional arguments to \code{\link{print}}.
}
}