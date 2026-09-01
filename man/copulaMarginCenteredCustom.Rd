\name{copulaMarginCenteredCustom}
\alias{copulaMarginCenteredCustom}
\title{Construct a Centered User-Defined Eta Margin}
\usage{
copulaMarginCenteredCustom(
  name, parameters, lower, upper, log_density, cdf, quantile,
  center, scale, support_fixed, random = NULL, set_scale = NULL,
  free = rep(TRUE, length(parameters)), roles = NULL, metadata = list()
)
}
\arguments{
\item{name}{Distribution name used in output and diagnostics.}
\item{parameters}{Named finite initial parameter vector.}
\item{lower,upper}{Lower and upper bounds for the distribution parameters.}
\item{log_density,cdf,quantile}{Functions of a raw variate and the named
parameter vector.}
\item{center}{Function returning the raw distribution mean. The constructed
eta is the raw variate minus this value.}
\item{scale}{Function returning the actual marginal standard deviation.}
\item{support_fixed}{Whether the support of the centered eta is independent of
the free parameters. This declaration is required.}
\item{random}{Optional raw-distribution random generator with arguments
\code{n} and the parameter vector.}
\item{set_scale}{Optional function replacing the marginal standard deviation
while retaining other native parameters.}
\item{free}{Logical mask selecting estimated parameters.}
\item{roles}{Parameter roles: \code{"scale"} or \code{"shape"}; location
parameters are not allowed for a centered eta margin.}
\item{metadata}{Additional margin metadata.}
}
\description{
Adapts an arbitrary proper continuous distribution to the centered additive
eta representation used by flexible-margin Gaussian-copula FREM. Population
location remains in the saemix fixed effects, while the custom distribution
controls eta shape and scale.
}
\details{
Exact distributional FFEM translation requires a valid continuous CDF and
quantile function. Likelihood-score estimation additionally requires a proper
log density, identifiable finite-dimensional parameters, finite variance, and
the smoothness and moment conditions stated for the estimator. A fixed support
must be declared with \code{support_fixed = TRUE}; parameter-dependent support
uses the fixed-reference quantile route and its additional regularity
conditions.

This constructor does not make unrestricted discrete distributions,
infinite-variance distributions, or nonidentifiable mixtures compatible with
the current eta covariance representation.
}
\value{A centered continuous \code{saemix_copula_margin}.}
\seealso{\code{\link{copulaPopulation}}, \code{\link{copulaFremToFfem}}}
\keyword{models}
