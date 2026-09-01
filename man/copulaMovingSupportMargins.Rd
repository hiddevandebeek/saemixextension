\name{copulaMovingSupportMargins}
\alias{copulaMarginUniformCentered}
\alias{copulaMarginCenteredGamma}
\alias{copulaMarginCovariateUniform}
\title{Moving-Support Margins for Gaussian-Copula FREM}
\usage{
copulaMarginUniformCentered(sd = 1)

copulaMarginCenteredGamma(shape = 4, sd = 1)

copulaMarginCovariateUniform(lower, upper)
}
\arguments{
\item{sd}{Positive marginal standard deviation.}
\item{shape}{Positive Gamma shape parameter.}
\item{lower,upper}{Finite increasing initial Uniform bounds.}
}
\description{
Construct continuous margins whose native support depends on estimated
parameters. The score estimator represents latent or missing values by their
fixed percentile and propagates the corresponding candidate quantile through
the structural model.
}
\details{
The centered constructors have mean zero and are therefore suitable for
additive eta coordinates without duplicating the saemix population location.
The natural-scale Uniform constructor is intended for a FREM covariate.

These margins declare \code{support_fixed = FALSE}. The likelihood-score
theorem applies while observed values remain in the interior of the support and
the quantile and response paths satisfy its differentiability and domination
conditions. It does not supply nonregular endpoint asymptotics when an observed
likelihood is maximized at a support boundary.
}
\value{A \code{saemix_copula_margin} object.}
\seealso{\code{\link{gaussianCopulaFrem}}}
\keyword{models}
