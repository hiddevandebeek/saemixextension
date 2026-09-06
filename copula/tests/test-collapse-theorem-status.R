suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})

set.seed(280827L)
n <- 220L; d <- 3L; dEta <- 2L
R <- matrix(c(1, .3, .4, .3, 1, -.2, .4, -.2, 1), d, d)
structure <- rvinecopulib::cvine_structure(c(3, 1, 2))
vine <- copulaGaussianRvineFromCor(R, structure)
startVine <- copulaGaussianRvineFromCor(diag(d), structure)
w <- rep(1 / n, n)

normalMargins <- list(copulaMarginNormal(.4), copulaMarginNormal(.3),
  copulaMarginCovariateNormal(75, 12))
E <- copulaMarginsQuantile(rvinecopulib::rvinecop(n, vine), normalMargins)
full <- copulaMaximiseGaussianFrem(E, w, normalMargins, startVine, d, dEta,
  maxit = 100L)
stopifnot(identical(full$backend, "gaussian-copula-frem-mvnormal-em"),
  full$theorem8$theoremAligned, full$theorem8$exactMstep,
  identical(full$theorem8$route, "M5/Theorem5"))

missing <- E
missing[seq(3L, n, by = 7L), 3L] <- NA_real_
missingFit <- copulaMaximiseGaussianFrem(missing, w, normalMargins,
  startVine, d, dEta, maxit = 100L)
stopifnot(identical(missingFit$backend,
    "gaussian-copula-frem-mvnormal-em"),
  !missingFit$theorem8$theoremAligned,
  missingFit$theorem8$fixedMap,
  identical(missingFit$theorem8$route,
    "conditional-until-MAX1-MAX2-certified"))

gammaMargins <- list(copulaMarginNormal(.4), copulaMarginNormal(.3),
  copulaMarginCovariateGamma(2, 35))
Eg <- copulaMarginsQuantile(rvinecopulib::rvinecop(n, vine), gammaMargins)
gammaFit <- copulaMaximiseGaussianFrem(Eg, w, gammaMargins, startVine,
  d, dEta, maxit = 150L)
stopifnot(identical(gammaFit$backend,
    "gaussian-copula-frem-one-arbitrary-profile"),
  gammaFit$theorem8$exactObjectiveReduction,
  !gammaFit$theorem8$exactMstep,
  !gammaFit$theorem8$theoremAligned)

## Exact conditional augmentation restores a complete Gaussian statistic for
## missing continuous covariates.
set.seed(280828L)
nAug <- 40000L; rho <- .55
Raug <- matrix(c(1, rho, rho, 1), 2, 2)
vineAug <- copulaGaussianRvineFromCor(Raug,
  rvinecopulib::cvine_structure(c(2, 1)))
marginAug <- list(copulaMarginNormal(.4),
  copulaMarginCovariateNormal(75, 12))
etaAug <- matrix(rep(.7, nAug), ncol = 1L)
covAug <- matrix(NA_real_, nAug, 1L)
imputed <- copulaGaussianFremImputeMissingConditioning(
  etaAug, covAug, vineAug, marginAug, 1L)
zImputed <- stats::qnorm(marginAug[[2L]]$cdf(
  imputed[, 1L], marginAug[[2L]]$parameters))
expectedMean <- rho * (.7 / .4); expectedVariance <- 1 - rho^2
stopifnot(abs(mean(zImputed) - expectedMean) < .012,
  abs(stats::var(zImputed) - expectedVariance) < .012)

## The pool updater takes that augmentation automatically and therefore sends
## a complete sufficient statistic to the exact Gaussian M-step.
set.seed(280829L)
nPool <- 160L
joint <- copulaMarginsQuantile(rvinecopulib::rvinecop(nPool, vineAug), marginAug)
conditioning <- matrix(joint[, 2L], ncol = 1L)
conditioning[seq(2L, nPool, by = 3L), 1L] <- NA_real_
copulaSet(vineAug, margins = marginAug, conditioning = conditioning,
  warmStartOnActivate = FALSE, guard = FALSE,
  augmentMissingGaussian = TRUE)
on.exit(copulaClear(), add = TRUE)
copulaPoolUpdate(matrix(joint[, 1L], ncol = 1L), gamma = 1,
  nchains = 1L, subject = seq_len(nPool))
finitePool <- copulaGaussianFremStatsMaterialize(
  .cop$gaussianFremStats, .cop$margins)
stopifnot(!is.null(finitePool), !anyNA(finitePool$E))
poolFit <- copulaMaximiseGaussianFrem(finitePool$E,
  finitePool$w / sum(finitePool$w), marginAug, vineAug, 2L, 1L,
  maxit = 100L)
stopifnot(poolFit$theorem8$theoremAligned,
  poolFit$theorem8$exactMstep)
interceptDesign <- matrix(1, nrow(finitePool$E), 1L)
locationMap <- matrix(c(1, 0), 1L, 2L)
designedPoolFit <- copulaMaximiseGaussianFrem(finitePool$E,
  finitePool$w / sum(finitePool$w), marginAug, vineAug, 2L, 1L,
  maxit = 100L, X = interceptDesign, locMap = locationMap,
  beta0 = 0, betaFree = 1L)
stopifnot(designedPoolFit$theorem8$theoremAligned,
  designedPoolFit$theorem8$exactMstep)

cat("collapse-path Theorem status checks passed\n")
