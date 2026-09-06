## Exactness gate for the Normal/lognormal sufficient-statistic M-step.
suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})

set.seed(82631)
d <- 4L; dEta <- 2L; n <- 600L
R <- matrix(c(1, .25, .42, -.18,
              .25, 1, -.20, .36,
              .42, -.20, 1, .12,
              -.18, .36, .12, 1), d, d, byrow = TRUE)
stopifnot(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) > 0)
structure <- rvinecopulib::cvine_structure(c(3, 1, 4, 2))
vine <- copulaGaussianRvineFromCor(R, structure)
margins <- list(copulaMarginNormal(.30), copulaMarginNormal(.42),
  copulaMarginCovariateLognormal(log(70), .21),
  copulaMarginCovariateNormal(92, 16))
joint <- copulaMarginsQuantile(rvinecopulib::rvinecop(n, vine), margins)
w <- rexp(n); w <- w / sum(w)

start <- list(copulaMarginNormal(.38), copulaMarginNormal(.34),
  copulaMarginCovariateLognormal(log(66), .28),
  copulaMarginCovariateNormal(98, 20))
fit <- copulaMaximiseGaussianFremMvnorm(joint, w, start, vine, d, dEta,
  maxit = 40L, withMu = TRUE)
stopifnot(identical(fit$backend, "gaussian-copula-frem-mvnormal-em"))

Y <- joint; Y[, 3] <- log(Y[, 3])
mu <- colSums(w * Y)
S <- crossprod(sweep(Y, 2L, mu, "-"), w * sweep(Y, 2L, mu, "-"))
fittedMu <- c(fit$delta[1:2], fit$margins[[3]]$parameters["meanlog"],
              fit$margins[[4]]$parameters["mean"])
fittedScale <- c(fit$margins[[1]]$parameters["sd"],
  fit$margins[[2]]$parameters["sd"], fit$margins[[3]]$parameters["sdlog"],
  fit$margins[[4]]$parameters["sd"])
fittedS <- diag(fittedScale) %*% copulaGaussianRvineCor(fit$vine) %*%
  diag(fittedScale)
stopifnot(max(abs(fittedMu - mu)) < 2e-8,
          max(abs(fittedS - S)) < 2e-8)

## Missing-covariate EM must improve the exact observed-data common Q.
missing <- joint
missing[seq(3, n, 5), 3] <- NA_real_
missing[seq(2, n, 7), 4] <- NA_real_
missingFit <- copulaMaximiseGaussianFremMvnorm(
  missing, w, start, vine, d, dEta, maxit = 80L, withMu = TRUE)
q0 <- sum(w * copulaGaussianFremLogPrior(missing, vine, start, dEta))
q1 <- sum(w * copulaGaussianFremLogPrior(
  sweep(missing, 2L, missingFit$delta, "-"), missingFit$vine,
  missingFit$margins, dEta))
stopifnot(is.finite(q1), q1 + 1e-10 >= q0)

## Moment compression must return the identical M-step, not merely a nearby Q.
literalFit <- copulaMaximiseGaussianFremMvnorm(
  missing, w, start, vine, d, dEta, maxit = 80L, withMu = TRUE,
  compressMoments = FALSE)
stopifnot(max(abs(literalFit$delta - missingFit$delta)) < 2e-8,
  max(abs(copulaGaussianRvineCor(literalFit$vine) -
    copulaGaussianRvineCor(missingFit$vine))) < 2e-8,
  abs(literalFit$value - missingFit$value) < 2e-8)

## A varying continuous design is part of the same quadratic statistic.  It
## must not force one compressed group per distinct design row.
X <- cbind(intercept = 1, exposure = seq(-1, 1, length.out = n))
locMap <- matrix(0, 2L, d)
locMap[1L, 1L] <- 1
locMap[2L, 2L] <- 1
betaTruth <- c(.18, -.11); betaStart <- c(0, 0)
designed <- missing + copulaLocation(X, betaTruth, locMap)
designedLiteral <- copulaMaximiseGaussianFremMvnorm(
  designed, w, start, vine, d, dEta, maxit = 80L, withMu = TRUE,
  X = X, locMap = locMap, beta0 = betaStart, compressMoments = FALSE)
designedFinite <- copulaMaximiseGaussianFremMvnorm(
  designed, w, start, vine, d, dEta, maxit = 80L, withMu = TRUE,
  X = X, locMap = locMap, beta0 = betaStart, compressMoments = TRUE)
compressedDesign <- copulaGaussianFremMomentCompress(
  { z <- designed; z[, 3L] <- log(z[, 3L]); z }, w, X)
stopifnot(compressedDesign$compressedRows < n / 10,
  max(abs(designedLiteral$beta - designedFinite$beta)) < 2e-8,
  max(abs(copulaGaussianRvineCor(designedLiteral$vine) -
    copulaGaussianRvineCor(designedFinite$vine))) < 2e-8,
  abs(designedLiteral$value - designedFinite$value) < 2e-8)

## The online statistic follows the same SA recursion as the explicit pool.
stats <- NULL; poolE <- matrix(numeric(), 0L, d); poolW <- numeric()
for (k in seq_len(8L)) {
  rows <- ((k - 1L) * 50L + 1L):(k * 50L)
  batch <- missing[rows, , drop = FALSE]
  gamma <- if (k <= 3L) 1 else 1 / (k - 2L)
  stats <- copulaGaussianFremStatsUpdate(stats, batch, NULL, gamma, 1L,
    start, dEta)
  poolW <- c((1 - gamma) * poolW, rep(gamma, nrow(batch)))
  poolE <- rbind(poolE, batch)
  keep <- poolW > 0; poolW <- poolW[keep]; poolE <- poolE[keep, , drop = FALSE]
}
finite <- copulaGaussianFremStatsMaterialize(stats, start)
explicitFit <- copulaMaximiseGaussianFremMvnorm(poolE, poolW, start, vine,
  d, dEta, maxit = 80L, withMu = TRUE, compressMoments = FALSE)
finiteFit <- copulaMaximiseGaussianFremMvnorm(finite$E, finite$w, start, vine,
  d, dEta, maxit = 80L, withMu = TRUE, compressMoments = FALSE)
stopifnot(max(abs(explicitFit$delta - finiteFit$delta)) < 2e-8,
  max(abs(copulaGaussianRvineCor(explicitFit$vine) -
    copulaGaussianRvineCor(finiteFit$vine))) < 2e-8,
  abs(explicitFit$value - finiteFit$value) < 2e-8)

## The recursive form also preserves a varying design without storing its raw
## rows. Compare it with the explicit weighted history.
statsX <- NULL; poolEX <- matrix(numeric(), 0L, d)
poolXX <- matrix(numeric(), 0L, ncol(X)); poolWX <- numeric()
for (k in seq_len(8L)) {
  rows <- ((k - 1L) * 50L + 1L):(k * 50L)
  gamma <- if (k <= 3L) 1 else 1 / (k - 2L)
  statsX <- copulaGaussianFremStatsUpdate(statsX,
    designed[rows, , drop = FALSE], X[rows, , drop = FALSE], gamma, 1L,
    start, dEta)
  poolWX <- c((1 - gamma) * poolWX, rep(gamma, length(rows)))
  poolEX <- rbind(poolEX, designed[rows, , drop = FALSE])
  poolXX <- rbind(poolXX, X[rows, , drop = FALSE])
  keep <- poolWX > 0; poolWX <- poolWX[keep]
  poolEX <- poolEX[keep, , drop = FALSE]; poolXX <- poolXX[keep, , drop = FALSE]
}
finiteX <- copulaGaussianFremStatsMaterialize(statsX, start)
explicitX <- copulaMaximiseGaussianFremMvnorm(poolEX, poolWX, start, vine,
  d, dEta, maxit = 80L, withMu = TRUE, X = poolXX, locMap = locMap,
  beta0 = betaStart, compressMoments = FALSE)
recursiveX <- copulaMaximiseGaussianFremMvnorm(finiteX$E, finiteX$w, start,
  vine, d, dEta, maxit = 80L, withMu = TRUE, X = finiteX$X,
  locMap = locMap, beta0 = betaStart, compressMoments = FALSE)
stopifnot(finiteX$rows < n / 10,
  max(abs(explicitX$beta - recursiveX$beta)) < 2e-8,
  max(abs(copulaGaussianRvineCor(explicitX$vine) -
    copulaGaussianRvineCor(recursiveX$vine))) < 2e-8,
  abs(explicitX$value - recursiveX$value) < 2e-8)

cat("Gaussian-FREM MVN sufficient-statistic M-step checks passed\n")
