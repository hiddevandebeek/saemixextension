\name{copulaFisherInformation}
\alias{copulaStandardErrors}
\alias{copulaReportedParameters}
\alias{copulaScoreFisherInformation}
\alias{copulaScoreSandwich}
\title{Uncertainty for Gaussian-copula FREM fits}
\description{
  Standard errors, confidence intervals and a covariance matrix for a fitted
  copula population model, estimated from the score function alone.

  For a flexible-margin fit this is the only uncertainty available:
  \code{\link{fim.saemix}} linearises the structural model around the
  conditional means and assumes Gaussian random effects, so it is refused for
  copula fits.

  \code{copulaStandardErrors} is the ordinary entry point.
  \code{copulaReportedParameters} returns the same information as a covariance
  matrix. \code{copulaScoreFisherInformation} recomputes the information after
  the fact from fresh posterior draws, which is slower and needs more of them,
  but does not require the fit to have carried the metric.
}
\usage{
copulaStandardErrors(object, level = 0.95)
copulaReportedParameters(object)
copulaScoreFisherInformation(object, draws = 1000L, max.iter = NULL,
  seed = 935001L, step = 1e-4, check = TRUE)
copulaScoreSandwich(object, draws = 200L, max.iter = NULL, seed = 935001L,
  step = 1e-4, hessianStep = 1e-3)
}
\arguments{
  \item{object}{A fitted \code{SaemixObject} with a Gaussian-copula
    population.}
  \item{level}{Confidence level for the reported interval.}
  \item{draws, max.iter, seed}{Posterior draws per subject, MCMC iterations
    behind them, and the seed, for the post-hoc estimator.}
  \item{step, hessianStep}{Relative steps for the numerical first and second
    derivatives of the per-subject log density.}
  \item{check}{Whether to repeat the estimate on a halved chain and warn when
    the answer has not settled. See the note on posterior mixing below.}
}
\details{
  The estimator is that of Delattre and Kuhn (2023): with \eqn{f_i} the
  complete-data density of subject \eqn{i},
  \deqn{I = \frac{1}{N}\sum_i E[\nabla \log f_i \mid y_i]\,
    E[\nabla \log f_i \mid y_i]^{T},}
  which uses only first derivatives and the posterior expectations that
  Fisher's identity already needs for the score. Its inverse is the asymptotic
  covariance of the estimate.

  During a fit the per-subject scores are accumulated along the recursion,
  \eqn{\Delta_i^k = (1-\gamma_k)\Delta_i^{k-1} + \gamma_k s_i^k}, in the manner
  of Fisher-SGD (Baey et al. 2023). That matters for accuracy rather than
  convenience: for a population location the complete-data score is large while
  the observed-data score is small, because the random effect absorbs a shift
  in the location, so the whole cancellation depends on the posterior average.
  The posterior draws mix slowly -- an effective sample fraction of a few
  percent that does not improve when the chain is lengthened -- so a single
  batch of draws understates a location standard error badly, while a fit
  visits a new MCMC state every iteration.

  The same matrix serves as the preconditioner. It is used by default on the
  parameter scale and not on the transformed-additive scale, which is an
  empirical choice: a metric cannot move the stationary points of the marginal
  likelihood, so it only affects how quickly they are reached, and the
  measured answer differs by problem.

  The reported coordinates are not the ones the recursion steps in. Typical
  values are returned on their own scale rather than the working scale,
  dependence as pairwise correlations rather than Cholesky angles, and each
  marginal parameter under the name of the coordinate it belongs to; the
  covariance is carried across by the delta method.

  Being the outer-product form, this estimates the covariance of the score,
  which coincides with the Fisher information only when the model is correctly
  specified. \code{copulaScoreSandwich} drops that requirement: it keeps the
  same matrix as the meat and takes the bread from Louis' identity,
  \deqn{-\nabla^2 \log g_i = E[-\nabla^2 \log f_i \mid y_i] -
    \mathrm{Var}[\nabla \log f_i \mid y_i],}
  both terms being posterior expectations over the draws the score already
  needs, and returns \eqn{H^{-1} J H^{-1} / N}.

  Measured against the spread of estimates across replicates: with a correct
  model the two agree, as the information equality requires; with a
  misspecified one -- a Gamma clearance fitted as Normal -- the residual-error
  interval went from 0.56 of the spread it should have covered to 0.95, while
  the other coordinates were comparable under both. So the sandwich repairs
  the case it is for and costs little in the case it is not.
}
\value{
  For \code{copulaStandardErrors}, a data frame with one row per reported
  parameter and columns \code{parameter}, \code{estimate}, \code{se},
  \code{lower} and \code{upper}. For \code{copulaReportedParameters}, a list
  with the estimate vector and its covariance matrix.
  \code{copulaScoreFisherInformation} additionally returns the information
  matrix, the per-subject scores and their Monte Carlo standard errors, and a
  \code{stability} ratio comparing the answer against a halved chain.
}
\references{
  Delattre, M. and Kuhn, E. (2023). Estimating Fisher information matrix in
  latent variable models based on the score function. arXiv:1909.06094.

  Baey, C., Delattre, M., Kuhn, E., Leger, J.-B. and Lemler, S. (2023).
  Efficient preconditioned stochastic gradient descent for estimation in latent
  variable models. ICML, PMLR v202.
}
\seealso{\code{\link{gaussianCopulaFrem}}, \code{\link{copulaNaturalMargins}}}
\examples{
\dontrun{
fit <- gaussianCopulaFrem(model, data, control)
copulaStandardErrors(fit)
}
}
\keyword{models}
