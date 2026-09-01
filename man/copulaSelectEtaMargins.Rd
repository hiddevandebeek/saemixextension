\name{copulaSelectEtaMargins}
\alias{copulaSelectEtaMargins}
\title{Screen Eta-Margin Families Using Full Posterior Draws}
\usage{
copulaSelectEtaMargins(object,
  candidates = c("normal", "student", "laplace"),
  posteriorDraws = 500L, max.iter = NULL, seed = 930001L,
  maxModels = 81L, optimizerMaxit = 100L,
  minimumEssFraction = 0.01, minimumPosteriorEssFraction = 0.05,
  maximumMcse = 0.5)
}
\arguments{
\item{object}{A fitted fixed Gaussian-copula FREM \code{SaemixObject}.}
\item{candidates}{Candidate family names, or one character vector per eta
coordinate. Supported automatic candidates are \code{"normal"},
\code{"student"}, and \code{"laplace"}.}
\item{posteriorDraws}{Number of full posterior draws per subject in each of
two independently generated pools.}
\item{max.iter}{Optional MCMC iteration control passed to the individual
conditional sampler.}
\item{seed}{Seed for the fitting posterior pool. The validation pool uses
\code{seed + 1}.}
\item{maxModels}{Maximum allowed Cartesian candidate combinations.}
\item{optimizerMaxit}{Maximum iterations for each fixed-dependence posterior
bridge optimization.}
\item{minimumEssFraction}{Minimum subject-level validation importance-sampling
ESS fraction required for an eligible candidate.}
\item{minimumPosteriorEssFraction}{Minimum MCMC ESS fraction across subjects
and eta coordinates required for a resolved automatic ranking.}
\item{maximumMcse}{Maximum estimated validation log-likelihood Monte Carlo
standard error required for eligibility.}
}
\description{
Uses full subject-specific posterior draws from a fitted Gaussian incumbent to
screen eta-margin families by an exact likelihood-ratio identity. Family
selection occurs outside score stochastic approximation. The returned
population declaration has a fixed selected family fingerprint and is intended
for a fresh \code{saemix()} fit.
}
\details{
For incumbent population density \eqn{p_0} and candidate \eqn{p_m}, the
subject-specific observed-likelihood ratio is
\deqn{g_m(y_i,c_i)/g_0(y_i,c_i)=
E_0[p_m(\eta_i,c_i)/p_0(\eta_i,c_i)\mid y_i,c_i].}
Candidate parameters are optimized with one posterior pool and ranked with an
independent pool. The Gaussian dependence and all non-eta margins are held at
their incumbent values during screening. The final fresh fit re-estimates all
continuous parameters within the selected fixed model.

Heavy-tailed candidates can have unstable importance weights under a Gaussian
incumbent. Candidates failing the ESS or MCSE gates are ineligible. Increasing
\code{posteriorDraws}, improving MCMC mixing, or using a future defensive
proposal is required when overlap is insufficient.

The returned \code{diagnostics$selectionResolved} is true only when every
candidate is eligible and the leading validation-BIC advantage is separated
from its runner-up relative to Monte Carlo error. An unresolved result should
not be used for an automatic final fit; increase \code{posteriorDraws} or reduce
the scientifically prespecified candidate set.

The validation BIC advantage is relative to the incumbent and penalizes the
difference in the number of free eta-margin parameters. Screening does not
replace an independently evaluated observed likelihood after fresh refitting.
}
\value{
An object of class \code{saemixEtaMarginSelection} containing the ranking,
diagnostics, selected margins, \code{population} declaration for the winner,
and \code{populations} for all candidate combinations.
}
\examples{
\dontrun{
incumbent <- saemix(model, data, control,
  population = gaussianCopulaFrem(
    etaSd = c(.25, .30), covariates = subjectCovariates))
selection <- copulaSelectEtaMargins(incumbent)
final <- saemix(model, data, control, population = selection$population)
}
}
\seealso{\code{\link{gaussianCopulaFrem}}, \code{\link{saemix}}}
\keyword{models}
