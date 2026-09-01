\name{copulaFremToFfem}
\alias{copulaFremToFfem}
\alias{copulaFfemTransform}
\alias{copulaFfemLocation}
\alias{copulaFfemQuantile}
\alias{copulaFfemSimulate}
\title{Translate Gaussian-copula FREM to a fixed-effects covariate model}
\usage{
copulaFremToFfem(object, covariates = NULL,
  representation = c("auto", "linear", "distributional"))

copulaFfemTransform(object, newdata)

copulaFfemLocation(object, newdata)

copulaFfemQuantile(object, newdata,
  probabilities = c(0.05, 0.5, 0.95))

copulaFfemSimulate(object, newdata, n = 1L, seed = NULL)
}
\arguments{
  \item{object}{For \code{copulaFremToFfem}, a fitted \code{SaemixObject}
  or fitted copula state. For the other functions, an object returned by
  \code{copulaFremToFfem}.}
  \item{covariates}{Names, indices, or a logical selector for the covariate
  subset. The default selects all fitted covariates.}
  \item{representation}{\code{"auto"} chooses a conventional linear FFEM
  when the parameter margins are Normal and selected covariates are
  continuous; otherwise it chooses the exact distributional representation.
  \code{"linear"} fails if an exact conventional FFEM is unavailable.}
  \item{newdata}{A data frame or matrix containing the selected covariates.}
  \item{probabilities}{Conditional quantile probabilities strictly between
  zero and one.}
  \item{n}{Number of conditional draws per row of \code{newdata}.}
  \item{seed}{Optional random-number seed.}
}
\description{
Calculates an exact conditional fixed-effects representation for any selected
covariate subset of a fitted Gaussian-copula FREM model.
}
\details{
When the parameter random-effect margins are Normal and the selected covariates
are continuous, the returned coefficient matrix acts on covariate Gaussian
scores \eqn{z_c=\Phi^{-1}\{F_c(c)\}} and the returned \code{omega} is the
corresponding unexplained parameter covariance. If the selected covariates are
also Normal, \code{native} contains the ordinary raw-scale FFEM coefficients
and covariate centres.

For non-Normal parameter margins or ordered categorical covariates, an exact
ordinary additive-Normal FFEM does not generally exist. The distributional
representation retains exact conditional quantiles and simulation through the
fitted marginal quantile functions. It never silently substitutes a
moment-matched Normal approximation.

The translation describes the conditional population model given the selected
covariates. It does not retain the covariate-density contribution to the joint
FREM likelihood.
}
\value{
\code{copulaFremToFfem} returns an object of class
\code{"saemixCopulaFfem"}. For an exact linear FFEM it contains
\code{coefficient}, \code{omega}, and, where available, raw-scale
\code{native} coefficients. The helper functions return transformed
covariates, conditional locations, conditional quantiles, or conditional
random-effect draws.
}
\examples{
\dontrun{
ffem <- copulaFremToFfem(fit, c("WT", "eGFR"))
z <- copulaFfemTransform(ffem, newSubjects)
conditionalMean <- copulaFfemLocation(ffem, newSubjects)
eta <- copulaFfemSimulate(ffem, newSubjects, n = 100)
}
}
