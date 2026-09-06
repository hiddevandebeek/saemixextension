suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})

## Stable coordinate transforms are exact inverses across every bound type.
native <- c(.4, 2, -3, .25)
lower <- c(-Inf, 0, -Inf, 0)
upper <- c( Inf, Inf, 0, 1)
internal <- copulaScoreToInternal(native, lower, upper)
stopifnot(max(abs(copulaScoreFromInternal(internal, lower, upper) - native)) < 1e-12)

set.seed(280834L)
n <- 100L; rho <- .5
R <- matrix(c(1, rho, rho, 1), 2, 2)
structure <- rvinecopulib::cvine_structure(c(2, 1))
truthVine <- copulaGaussianRvineFromCor(R, structure)
truthMargins <- list(copulaMarginNormal(.4),
  copulaMarginCovariateGamma(1.4, 48))
E <- copulaMarginsQuantile(rvinecopulib::rvinecop(n, truthVine), truthMargins)
E[, 1L] <- E[, 1L] + .12
conditioning <- matrix(E[, 2L], ncol = 1L)
beta <- mean(E[, 1L])
startMargins <- list(copulaMarginNormal(stats::sd(E[, 1L])),
  copulaFitCovariateMargin(E[, 2L], candidates = "gamma"))
initialResidual <- cbind(E[, 1L] - beta, E[, 2L])
initialScores <- copulaGaussianFremEvaluateMargins(initialResidual,
  startMargins)$z
startVine <- copulaGaussianRvineFromCor(stats::cor(initialScores), structure)
X <- matrix(1, n, 1L); locMap <- matrix(c(1, 0), 1L, 2L)
w <- rep(1 / n, n)
reference <- copulaMaximiseGaussianFrem(E, w, startMargins, startVine,
  2L, 1L, maxit = 400L, X = X, locMap = locMap,
  beta0 = beta, betaFree = 1L)

copulaSet(startVine, margins = startMargins, conditioning = conditioning,
  warmStartOnActivate = FALSE, guard = FALSE,
  populationAlgorithm = "score-sa", scoreScale = .07)
on.exit(copulaClear(), add = TRUE)
for (k in seq_len(350L)) {
  gain <- if (k <= 40L) 1 else (k - 40L)^-.62
  copulaScoreBatchUpdate(matrix(E[, 1L] - beta, ncol = 1L), 1L,
    phi = matrix(E[, 1L], ncol = 1L), design = X,
    locationMap = matrix(1, 1L, 1L), beta = beta, betaFree = 1L,
    subject = seq_len(n))
  copulaMstep(k, gamma = gain)
  updated <- copulaTakeBeta()
  if (!is.null(updated)) beta <- updated
}
state <- copulaGet()
residual <- cbind(E[, 1L] - beta, E[, 2L])
value <- sum(w * copulaGaussianFremLogPrior(residual,
  state$vine, state$margins, 1L, "joint"))
stopifnot(identical(state$lastJoint$backend,
    "gaussian-copula-frem-score-sa"),
  state$scoreState$projectionCount == 0L,
  isTRUE(state$lastJoint$scoreTheory$metricFrozen),
  state$lastJoint$postFreezeProjectionCount == 0L,
  state$lastJoint$postFreezeBacktrackCount == 0L,
  state$lastJoint$postFreezeNoMoveCount == 0L,
  reference$value - value < 2e-4,
  max(abs(copulaGaussianRvineCor(state$vine) -
    copulaGaussianRvineCor(reference$vine))) < 3e-3)

bad <- try(copulaSet(startVine, margins = startMargins,
  conditioning = conditioning, warmStartOnActivate = FALSE,
  populationAlgorithm = "score-sa", refitEvery = 2L), silent = TRUE)
stopifnot(inherits(bad, "try-error"))

cat(sprintf("score-SA backend check passed; common-Q gap %.3g\n",
  reference$value - value))
