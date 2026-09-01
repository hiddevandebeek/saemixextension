\name{copulaPopulation}
\alias{copulaPopulation}
\alias{copulaClear}
\alias{copulaGet}
\title{Declare a Fixed Gaussian-Copula Population Model}
\usage{
copulaPopulation(vine, margins = NULL, sd = NULL,
  scale = c("auto", "transformed-additive"), conditioning = NULL,
  populationAlgorithm = "score-sa", scoreScale = 0.05,
  scoreFiniteDifference = 1e-4, scoreProjection = 24,
  scoreGainScale = 0.2, scoreGainPower = 0.8,
  scoreGainOffset = 30, scoreBurn = 50L, ...)

copulaClear()

copulaGet(object = NULL)
}
\arguments{
\item{vine}{A full Gaussian \code{rvinecopulib::vinecop_dist}.}
\item{margins}{One declared marginal distribution per vine coordinate.}
\item{sd}{Positive Normal-margin scales used when \code{margins} is omitted.}
\item{scale}{Population scale. Only transformed-additive margins are supported.}
\item{conditioning}{Optional covariate matrix, or a list containing
\code{values} and \code{variableName}. Missing values may be \code{NA}.}
\item{populationAlgorithm}{The fixed value \code{"score-sa"}.}
\item{scoreScale,scoreFiniteDifference,scoreProjection,scoreGainScale,
scoreGainPower,scoreGainOffset,scoreBurn}{Controls for the likelihood-score
stochastic-approximation recursion.}
\item{...}{Reserved. Unsupported controls are rejected when the model is used.}
\item{object}{An optional fitted \code{SaemixObject}.}
}
\description{
Creates the population declaration consumed by the \code{population} argument
of \code{saemix()}. The fitted model has fixed marginal families, a fixed full
Gaussian copula, and the joint observed-data likelihood.
}
\details{
The population estimator is the controlled-MCMC likelihood-score recursion.
Family or structure selection, conditional-likelihood fitting, direct
natural-parameter margins, and general non-Gaussian R-vines are not part of
this interface. Use \code{gaussianCopulaFrem()} for the usual constructor.

\code{copulaGet(object)} returns the immutable fitted population snapshot.
The no-argument form and \code{copulaClear()} support low-level tests; ordinary
user code should pass the population object directly to \code{saemix()}.
}
\value{
\code{copulaPopulation()} returns a \code{saemixPopulation} declaration.
\code{copulaGet()} returns a population state list. \code{copulaClear()} is
called for its side effect.
}
\seealso{\code{\link{gaussianCopulaFrem}}, \code{\link{saemix}}}
\keyword{models}
