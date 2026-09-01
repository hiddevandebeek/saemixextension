\name{copulaMargins}
\alias{copulaMargin}
\alias{copulaMarginDistribution}
\alias{copulaMarginNormal}
\alias{copulaMarginStudent}
\alias{copulaMarginLaplace}
\alias{copulaMarginCovariateNormal}
\alias{copulaMarginCovariateLognormal}
\alias{copulaMarginCovariateGamma}
\alias{copulaMarginCovariateWeibull}
\alias{copulaFitCovariateMargin}
\title{Marginal Distributions for Gaussian-Copula FREM}
\usage{
copulaMargin(name, parameters, lower, upper,
  free = rep(TRUE, length(parameters)), log_density, cdf, quantile,
  random = NULL, type = c("continuous", "discrete"), cdf_left = NULL,
  scale = NULL, set_scale = NULL, centered = TRUE, roles = NULL,
  scale_is_sd = FALSE, metadata = list(), support_fixed = NULL)

copulaMarginDistribution(distr, parameters, lower, upper,
  free = rep(TRUE, length(parameters)),
  type = c("continuous", "discrete"), centered = FALSE,
  scale = NULL, set_scale = NULL, cdf_left = NULL, roles = NULL,
  scale_is_sd = FALSE, support_fixed = NULL)

copulaMarginNormal(sd = 1)
copulaMarginStudent(sd = 1, df = 6)
copulaMarginLaplace(sd = 1)
copulaMarginCovariateNormal(mean, sd)
copulaMarginCovariateLognormal(meanlog, sdlog)
copulaMarginCovariateGamma(shape, scale)
copulaMarginCovariateWeibull(shape, scale)
copulaFitCovariateMargin(x,
  candidates = c("auto", "normal", "lognormal", "gamma", "weibull"))
}
\arguments{
\item{name,distr}{Distribution name. The distribution adapter resolves
\code{d*}, \code{p*}, and \code{q*} functions by this name.}
\item{parameters}{Named finite initial parameter vector.}
\item{lower,upper}{Bounds for the distribution parameters.}
\item{free}{Logical mask selecting estimated parameters.}
\item{log_density,cdf,quantile}{Distribution callbacks accepting values and
the named parameter vector.}
\item{random}{Optional random generator.}
\item{type,cdf_left}{Continuous or discrete measure declaration and, for a
discrete margin, its left-limit CDF.}
\item{scale,set_scale,scale_is_sd}{Callbacks and declaration for the marginal
scale. Eta margins used by FREM must expose their actual finite standard
deviation.}
\item{centered}{Whether the margin is centred at zero for an additive eta
coordinate.}
\item{roles}{Native parameter roles: location, scale, or shape.}
\item{metadata}{Additional margin metadata.}
\item{support_fixed}{Whether native support is parameter independent.}
\item{sd,df,mean,meanlog,sdlog,shape}{Built-in marginal parameters.}
\item{x}{Observed continuous covariate values.}
\item{candidates}{Candidate covariate marginal families.}
}
\description{
Defines the univariate marginal transformations used by Gaussian-copula FREM.
Built-in constructors cover common eta and covariate families. Arbitrary
continuous families can be supplied through callbacks with
\code{copulaMargin()} or adapted from compatible distribution functions with
\code{copulaMarginDistribution()}.
}
\details{
Copula likelihood evaluation requires a proper density or probability mass,
CDF, and quantile function. Additive eta margins must be centred so the
population location is not duplicated. Standard FREM covariance and UPV
reporting additionally require a finite marginal variance.

For a user-defined eta family, \code{copulaMarginCenteredCustom()} is usually
safer because it centres the raw distribution automatically and requires an
explicit support contract.
}
\value{A \code{saemix_copula_margin}, or for
\code{copulaFitCovariateMargin}, a fitted margin.}
\seealso{\code{\link{copulaMarginCenteredCustom}},
\code{\link{copulaPopulation}}}
\keyword{models}
