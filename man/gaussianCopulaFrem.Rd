\name{gaussianCopulaFrem}
\alias{gaussianCopulaFrem}
\alias{copulaGaussianRvineFromCor}
\title{Declare a Gaussian-Copula FREM Population for saemix}
\usage{
gaussianCopulaFrem(etaSd = NULL, etaMargins = NULL,
  covariates = NULL, covariateMargins = "auto",
  correlation = NULL, structure = NULL, covariateNames = NULL,
  scoreScale = 0.05, scoreBurn = 50L,
  gainScale = 0.2, gainPower = 0.8,
  finiteDifference = 1e-4, projection = 24)

copulaGaussianRvineFromCor(correlation, structure = NULL)
}
\arguments{
\item{etaSd}{Positive initial standard deviation for each Normal random-effect
margin. Supply this or \code{etaMargins}.}
\item{etaMargins}{One declared fixed-support margin per random effect.}
\item{covariates}{Optional numeric matrix with one subject per row. Missing
values may be \code{NA}.}
\item{covariateMargins}{One declared margin per covariate, or \code{"auto"}
for marginal-AIC selection among Normal, lognormal, Gamma, and Weibull using
the observed non-missing values. Explicit declarations may contain multiple
dependent binary, ordinal, or ordered nominal categorical margins; categorical
meaning is not inferred automatically from numeric codes.}
\item{correlation}{Initial joint Gaussian-copula correlation matrix. The
default is independence.}
\item{structure}{Optional full R-vine structure used only to parameterize the
Gaussian correlation matrix.}
\item{covariateNames}{Optional covariate labels.}
\item{scoreScale,scoreBurn,gainScale,gainPower,finiteDifference,projection}{
Controls for the likelihood-score stochastic-approximation estimator. The
default gain power is theory-compatible and values at or below 0.75 are
rejected.}
}
\description{
Creates the minimal theorem-backed Gaussian-copula FREM population object for
the ordinary \code{saemix()} interface. Stock saemix behavior is unchanged
when its \code{population} argument is absent or \code{NULL}.
}
\details{
The returned model has a fixed full Gaussian copula, fixed declared marginal
families, joint observed-data likelihood, exact numerical policy, and the
likelihood-score estimator. It does not expose the legacy retained-particle,
general R-vine, adaptive-selection, or conditional-likelihood controls.

Automatic family selection applies to continuous covariates. The selected
family is fixed before stochastic-approximation estimation begins. Eta margins
are Normal when declared through \code{etaSd}; non-Normal eta families are
supplied explicitly through \code{etaMargins}.
After fitting this incumbent, \code{copulaSelectEtaMargins()} can screen
Normal, Student, and Laplace eta-family combinations from full posterior draws
and return a fixed declaration for a fresh fit.

Every custom margin must declare its support. Fixed-support margins use
\code{support_fixed = TRUE}. A proper continuous margin may use
\code{support_fixed = FALSE} when its CDF and quantile define a valid
fixed-percentile representation; the estimator then propagates the candidate
quantile through the structural model.

Multiple categorical coordinates use exact truncated-Gaussian conditional
sampling and fixed-support latent-uniform augmentation. Their score recursion
does not differentiate a numerical multivariate-Normal rectangle probability.

Moving-support likelihoods whose optimum lies on an observed support boundary
are nonregular. The convergence theorem applies to smooth interior trajectories,
not to endpoint asymptotics such as an unrestricted estimated Uniform bound.
}
\value{A \code{saemixPopulation} object to pass as
\code{population =} in \code{saemix()}.}
\examples{
\dontrun{
population <- gaussianCopulaFrem(
  etaSd = c(0.25, 0.30),
  covariates = subjectData[, c("WT", "eGFR")])
fit <- saemix(model, data, control, population = population)
}
}
\seealso{\code{\link{copulaPopulation}},
\code{\link{copulaSelectEtaMargins}}, \code{\link{saemix}}}
\keyword{models}
