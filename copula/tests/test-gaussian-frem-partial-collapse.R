## Exact one-arbitrary-margin profile: objective identity and speed/optimum
## comparison against the pre-existing literal all-parameter optimizer.
suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})
set.seed(91024L)
n <- 400L; d <- 4L; dEta <- 2L
R0 <- matrix(c(1, .28, .52, -.22, .28, 1, -.18, .36,
  .52, -.18, 1, .20, -.22, .36, .20, 1), d, d, byrow = TRUE)
R0 <- stats::cov2cor(copulaGaussianFremEnsurePd(R0, 1e-6))
structure <- rvinecopulib::cvine_structure(c(3, 1, 2, 4))
vine <- copulaGaussianRvineFromCor(R0, structure)
truth <- list(copulaMarginStudent(.48, 3.2), copulaMarginNormal(.30),
  copulaMarginCovariateLognormal(log(72), .32),
  copulaMarginCovariateNormal(90, 14))
E <- copulaMarginsQuantile(rvinecopulib::rvinecop(n, vine), truth)
E[, 1L] <- E[, 1L] + .12
E[, 2L] <- E[, 2L] - .08
w <- rep(1 / n, n)
startMargins <- list(copulaMarginStudent(.42, 5), copulaMarginNormal(.35),
  copulaMarginCovariateLognormal(log(70), .38),
  copulaMarginCovariateNormal(88, 16))
startVine <- copulaGaussianRvineFromCor(diag(d), structure)

profileTime <- system.time(profile <- copulaMaximiseGaussianFrem(E, w,
  startMargins, startVine, d, dEta, maxit = 250L))[["elapsed"]]
stopifnot(identical(profile$backend,
  "gaussian-copula-frem-one-arbitrary-profile"))
profileLiteral <- sum(w * copulaGaussianFremLogPrior(
  sweep(E, 2L, profile$delta, "-"), profile$vine, profile$margins,
  dEta, "joint"))
incomingLiteral <- sum(w * copulaGaussianFremLogPrior(
  E, startVine, startMargins, dEta, "joint"))
stopifnot(abs(profile$value - profileLiteral) < 1e-8,
  profileLiteral + copulaObjectiveTolerance(incomingLiteral) >= incomingLiteral)

gaussianMargins <- list(copulaMarginNormal(.42), copulaMarginNormal(.35),
  copulaMarginCovariateNormal(mean(E[, 3L]), sd(E[, 3L])),
  copulaMarginCovariateNormal(mean(E[, 4L]), sd(E[, 4L])))
gaussianTime <- system.time(gaussian <- copulaMaximiseGaussianFrem(E, w,
  gaussianMargins, startVine, d, dEta, maxit = 250L))[["elapsed"]]
stopifnot(identical(gaussian$backend, "gaussian-copula-frem-mvnormal-em"))

## An intercept design represents the same two eta locations while forcing the
## legacy literal path (the first partial implementation rejects any design).
Xliteral <- matrix(1, n, dEta)
locLiteral <- matrix(0, dEta, d); diag(locLiteral[, seq_len(dEta)]) <- 1
designed <- copulaMaximiseGaussianFrem(E, w, startMargins, startVine,
  d, dEta, maxit = 250L, X = Xliteral, locMap = locLiteral,
  beta0 = numeric(dEta))
stopifnot(identical(designed$backend,
  "gaussian-copula-frem-one-arbitrary-profile"),
  max(abs(designed$beta - profile$delta[seq_len(dEta)])) < 1e-8,
  abs(designed$value - profile$value) < 1e-8)
oldDisable <- getOption("saemix.copula.disablePartialCollapse")
options(saemix.copula.disablePartialCollapse = TRUE)
literalTime <- system.time(literal <- copulaMaximiseGaussianFrem(E, w,
  startMargins, startVine, d, dEta, maxit = 1000L,
  X = Xliteral, locMap = locLiteral, beta0 = numeric(dEta)))[["elapsed"]]
options(saemix.copula.disablePartialCollapse = oldDisable)
cat(sprintf(paste0("Gaussian %.3fs; profile %.3fs Q %.10f; literal ",
  "%.3fs Q %.10f; deltaQ %.3g\n"), gaussianTime,
  profileTime, profile$value, literalTime, literal$value,
  profile$value - literal$value))
## The exact profile should reach at least the finite-iteration literal optimum.
stopifnot(profile$value + 2e-6 >= literal$value)

## A genuinely non-affine observed covariate admits an exact recursive grouped
## statistic: subject mass, Gaussian-block first moments and cross-moments.
set.seed(91025L)
ns <- 40L; nchains <- 2L; batches <- 8L
gammaCov <- copulaMarginCovariateGamma(1.7, 40)
normalCov <- copulaMarginCovariateNormal(90, 18)
groupMargins <- list(copulaMarginNormal(.35), copulaMarginNormal(.32),
  gammaCov, normalCov)
conditioning <- cbind(
  gammaCov$quantile(runif(ns), gammaCov$parameters),
  normalCov$quantile(runif(ns), normalCov$parameters))
subject <- rep(seq_len(ns), each = nchains)
Xbatch <- matrix(1, ns * nchains, dEta)
locMap <- matrix(0, dEta, d); diag(locMap[, seq_len(dEta)]) <- 1
beta0 <- numeric(dEta)
groupStats <- NULL
poolE <- matrix(numeric(), 0L, d); poolX <- matrix(numeric(), 0L, dEta)
poolW <- numeric()
for (iteration in seq_len(batches)) {
  gain <- if (iteration <= 3L) 1 else 1 / (iteration - 2L)
  eta <- matrix(rnorm(ns * nchains * dEta), ncol = dEta) %*%
    diag(c(.35, .32))
  batch <- cbind(eta, conditioning[subject, , drop = FALSE])
  poolE <- rbind(poolE, batch); poolX <- rbind(poolX, Xbatch)
  poolW <- c((1 - gain) * poolW, rep(gain / nchains, nrow(batch)))
  keep <- poolW > 0
  poolE <- poolE[keep, , drop = FALSE]
  poolX <- poolX[keep, , drop = FALSE]
  poolW <- poolW[keep]
  groupStats <- copulaGaussianFremOneArbitraryStatsUpdate(
    groupStats, batch, Xbatch, subject, gain, nchains, groupMargins,
    dEta, locMap, beta0, seq_len(dEta))
  stopifnot(!is.null(groupStats))
}
finite <- copulaGaussianFremOneArbitraryStatsMaterialize(
  groupStats, groupMargins)
stopifnot(!is.null(finite), finite$rows == ns,
  finite$rows < nrow(poolE))
groupVine <- copulaGaussianRvineFromCor(diag(d), structure)
explicitGroupTime <- system.time(explicitGroup <- copulaMaximiseGaussianFremOneArbitrary(
  poolE, poolW, groupMargins, groupVine, d, dEta, maxit = 150L,
  X = poolX, locMap = locMap, beta0 = beta0, betaFree = seq_len(dEta)))[["elapsed"]]
finiteGroupTime <- system.time(finiteGroup <- copulaMaximiseGaussianFremOneArbitrary(
  finite$E, finite$w, groupMargins, groupVine, d, dEta, maxit = 150L,
  X = finite$X, locMap = locMap, beta0 = beta0, betaFree = seq_len(dEta),
  groupedStats = finite$groupedStats))[["elapsed"]]
stopifnot(!is.null(explicitGroup), !is.null(finiteGroup),
  abs(explicitGroup$value - finiteGroup$value) < 2e-8,
  max(abs(explicitGroup$beta - finiteGroup$beta)) < 2e-8,
  max(abs(unlist(lapply(explicitGroup$margins, `[[`, "parameters")) -
    unlist(lapply(finiteGroup$margins, `[[`, "parameters")))) < 2e-8,
  max(abs(copulaGaussianRvineCor(explicitGroup$vine) -
    copulaGaussianRvineCor(finiteGroup$vine))) < 2e-8)
cat(sprintf(paste0("Grouped non-affine statistic: %d pool rows -> %d exact rows; ",
  "explicit %.3fs, recursive %.3fs\n"), nrow(poolE), finite$rows,
  explicitGroupTime, finiteGroupTime))
