\name{copulaConvergence}
\alias{copulaConvergence}
\title{Did the score-SA recursion converge?}
\description{
  Reads back the per-iteration trace a Gaussian-copula fit already records and
  reports whether the recursion settled, rather than leaving the question to be
  answered by inspecting the estimates.
}
\usage{
copulaConvergence(object, tail = 0.2)
}
\arguments{
  \item{object}{A fitted \code{SaemixObject} whose population was estimated by
    the likelihood-score recursion.}
  \item{tail}{Fraction of the run treated as the tail, compared against the
    stretch immediately before it.}
}
\details{
  Four conditions are reported, each able to fail on its own.

  \emph{The score should be zero.} It is a Monte Carlo estimate, so it never is
  exactly, and comparing it against its own earlier average cannot distinguish
  a recursion that is still travelling from one sitting at its noise floor.
  When the fit carried per-subject scores the question is settled properly:
  they are independent across subjects, so their average has covariance
  \eqn{V/N}, and under the hypothesis that the observed score is zero
  \deqn{T = N\,\bar g' I^{-1} \bar g \sim \chi^2_p,}
  which is a test with a p-value rather than a ratio with a threshold. Both
  ingredients are already accumulated for the Fisher metric, so it costs
  nothing.

  \emph{The steps should shrink.} Stochastic approximation converges because
  the gain decays, so the distance moved per iteration should decay with it.
  The reported figure is the displacement ratio divided by the gain ratio, so
  a value near or below one means the steps are shrinking at least as fast as
  the schedule driving them.

  \emph{The truncation should be idle.} The Andrieu et al. framework tolerates
  excursions from the compact set, but the convergence result describes a tail
  that stays inside one, so reprojections after the adaptation freeze mean the
  realised trajectory is not that tail.

  \emph{The metric should have settled}, since the theory is written for a
  preconditioner that is eventually fixed.
}
\value{
  An object of class \code{saemixScoreConvergence} with the chi-square
  statistic and p-value where available, the score and displacement ratios,
  post-freeze counts of reprojections, step halvings and refused steps, the
  preconditioner conditioning, and a \code{settled} flag. Printing it lays the
  conditions out with the failing ones marked.
}
\seealso{\code{\link{copulaStandardErrors}}, \code{\link{gaussianCopulaFrem}}}
\examples{
\dontrun{
fit <- gaussianCopulaFrem(model, data, control)
copulaConvergence(fit)
}
}
\keyword{models}
