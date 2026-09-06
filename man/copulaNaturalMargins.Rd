\name{copulaNaturalMargins}
\alias{copulaNaturalMargin}
\alias{copulaNaturalMarginNormal}
\alias{copulaNaturalMarginStudent}
\alias{copulaNaturalMarginLaplace}
\alias{copulaNaturalMarginLognormal}
\alias{copulaNaturalMarginGamma}
\alias{copulaNaturalMarginWeibull}
\alias{copulaNaturalMarginGeneralizedGamma}
\alias{copulaNaturalMarginLognormalMixture}
\alias{copulaNaturalMarginSNP}
\alias{copulaNaturalMarginLogitNormal}
\alias{copulaNaturalMarginBeta}
\alias{copulaNaturalMarginFamilies}
\alias{copulaRegisterNaturalMargin}
\alias{copulaNaturalMarginTilted}
\alias{copulaNaturalMarginTiltedFit}
\alias{copulaSelectParameterMargins}
\alias{copulaSelectTiltedMargins}
\title{Natural-scale marginal distributions for individual parameters}
\description{
  Marginal laws declared directly on the natural individual-parameter scale.
  The distribution is independent of the working coordinate used by the
  estimation software: that transformation is retained only as a reversible
  coordinate map with its exact Jacobian, so a Gamma clearance stays Gamma
  whether the internal coordinate is the identity or the logarithm.

  Each margin is anchored to the typical value, so it contributes shape while
  the population location remains identified. \code{copulaNaturalMargin}
  builds a margin from user-supplied functions;
  \code{copulaRegisterNaturalMargin} adds one to the screening registry so it
  can be ranked by \code{copulaSelectParameterMargins}.
}
\usage{
copulaNaturalMargin(name, support, anchor, parameters, lower, upper,
  log_density, cdf, quantile, random = NULL, score = NULL,
  inverse_score = NULL, evaluate = NULL,
  free = rep(TRUE, length(parameters)), roles = NULL, metadata = list())

copulaNaturalMarginNormal(sd = 1)
copulaNaturalMarginStudent(sd = 1, df = 6)
copulaNaturalMarginLaplace(sd = 1)
copulaNaturalMarginLognormal(sdlog = .3,
  anchor = c("median", "mean", "geometric-mean"))
copulaNaturalMarginGamma(shape = 4,
  anchor = c("geometric-mean", "median", "mean"))
copulaNaturalMarginWeibull(shape = 3,
  anchor = c("geometric-mean", "median", "mean"))
copulaNaturalMarginLogitNormal(sdlogit = .5, anchor = c("median", "custom"))
copulaNaturalMarginBeta(precision = 10, anchor = c("mean", "median"))
copulaNaturalMarginGeneralizedGamma(sigma = .3, Q = 1,
  anchor = c("geometric-mean", "median"))
copulaNaturalMarginLognormalMixture(weight = .3, separation = .8,
  sdlog = .25, anchor = c("geometric-mean"))
copulaNaturalMarginSNP(degree = 2L, scale = .3,
  coefficients = NULL, support = c("positive", "real"))
copulaNaturalMarginTilted(degree = 2L, scale = .3, beta = NULL,
  support = c("positive", "real"))
copulaNaturalMarginTiltedFit(x, typical, degree = 2L,
  support = c("positive", "real"), scalePoints = 9L)
copulaSelectTiltedMargins(object, supports, maxDegree = 4L,
  posteriorDraws = 500L, max.iter = NULL, seed = 935001L, maxSweeps = 3L,
  minimumEssFraction = .01, minimumPosteriorEssFraction = .05,
  maximumMcse = .5)

copulaNaturalMarginFamilies(support = c("real", "positive", "unit"))
copulaRegisterNaturalMargin(name, support, start, builtin = FALSE)
}
\arguments{
  \item{name}{Family name, lower case.}
  \item{support}{One of \code{"real"}, \code{"positive"} or \code{"unit"}.}
  \item{anchor}{How the typical value pins the location.}
  \item{parameters}{Named numeric vector of margin parameters.}
  \item{lower, upper}{Bounds for the parameters.}
  \item{log_density, cdf, quantile, random}{Functions of
    \code{(x, typical, par)} giving the natural-scale log density, CDF,
    quantile function and sampler.}
  \item{score}{Optional closed form for the Gaussian score
    \eqn{\Phi^{-1}\{F(x)\}}. Supplying it avoids forming the uniform, which is
    faster and does not lose the far tails to underflow. Defaults to the
    probability integral transform written out.}
  \item{inverse_score}{Optional closed form for the inverse of \code{score}.}
  \item{evaluate}{Optional function returning \code{log_density} and
    \code{score} together, so quantities shared by both are formed once.}
  \item{free, roles, metadata}{Which parameters are estimated, their roles,
    and any extra information carried with the margin.}
  \item{sd, df, sdlog, shape, sdlogit, precision, Q, weight, separation, degree, coefficients}{Family parameters.}
  \item{beta}{Tilt coefficients; \code{NULL} gives the base family exactly.}
  \item{x, typical}{Values and their typical values, for a direct fit.}
  \item{scalePoints}{Number of points in the sweep over the base scale.}
  \item{object}{A fitted \code{SaemixObject} to select margins for.}
  \item{supports, maxDegree, posteriorDraws, max.iter, seed, maxSweeps}{Support
    per eta coordinate, the highest tilt degree on the ladder, draws per
    posterior pool, MCMC iterations behind them, the seed, and the maximum
    number of coordinate sweeps.}
  \item{minimumEssFraction, minimumPosteriorEssFraction, maximumMcse}{Overlap
    and mixing diagnostics a degree must satisfy to be eligible.}
  \item{start}{Function of \code{(values, typical)} returning a starting
    margin, used when the family is screened.}
  \item{builtin}{Internal flag marking the shipped families.}
}
\details{
  The shipped registry contains \code{normal}, \code{student} and
  \code{laplace} on the real line, \code{lognormal}, \code{gamma} and
  \code{weibull} on the positive line, and \code{logitnormal} and \code{beta}
  on the unit interval. Further proper, absolutely continuous families can be
  registered; they must be identifiable given the anchor and satisfy the
  differentiability and moment conditions the estimator assumes. Families
  whose maximum likelihood estimate sits on a support boundary are out of
  scope.

  A supplied \code{score} or \code{evaluate} is checked against the CDF it
  replaces when the margin is constructed, so an inconsistent closed form is
  rejected rather than used.

  \code{copulaNaturalMarginTilted} is a different kind of family: instead of a
  named shape it carries a base family multiplied by
  \eqn{\exp\{\sum_k \beta_k L_k\}}, with \eqn{L_k} orthonormal Legendre
  polynomials in the base probability integral transform. Degree 0 is exactly
  lognormal (or Normal), degrees are nested, and the statistics are bounded, so
  the family is a regular exponential family: at a fixed location the
  log-likelihood is concave in \eqn{\beta}, the score is
  \eqn{L_k - E_\beta[L_k]} and the Hessian is \eqn{-\mathrm{Var}_\beta(L)}.
  Anchoring ties the location to \eqn{\beta} and adds a rank-one term to both.
  \code{copulaNaturalMarginTiltedFit} maximises the anchored likelihood from
  \eqn{\beta = 0} using that analytic gradient, and sweeps only the
  one-dimensional base scale around it.
  \code{copulaSelectTiltedMargins} walks the degree ladder per eta coordinate,
  scoring each degree by the same exact posterior likelihood-ratio identity the
  family screen uses, penalised by BIC. Because the ladder is nested, there is
  no candidate list to supply and no \eqn{K^d} grid to enumerate.
}
\value{
  An object of class \code{saemix_natural_parameter_margin} for the
  constructors, a character vector of registered family names for
  \code{copulaNaturalMarginFamilies}.
}
\seealso{\code{\link{copulaMargins}}, \code{\link{gaussianCopulaFrem}}}
\examples{
\dontrun{
clearance <- copulaNaturalMarginGamma(2.5)
volume <- copulaNaturalMarginLognormal(0.3)
copulaNaturalMarginFamilies("positive")

## a bounded parameter, for example bioavailability
bioavailability <- copulaNaturalMarginLogitNormal(0.6)
}
}
\keyword{models}
