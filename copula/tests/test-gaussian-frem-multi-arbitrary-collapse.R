## Exact Gaussian-block profile with two simultaneously curved margins.
suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})

set.seed(260827L)
n <- 450L; d <- 4L; dEta <- 2L
Rtruth <- matrix(c(1, .32, .50, -.18,
  .32, 1, -.14, .46,
  .50, -.14, 1, .24,
  -.18, .46, .24, 1), d, d, byrow = TRUE)
Rtruth <- stats::cov2cor(copulaGaussianFremEnsurePd(Rtruth, 1e-6))
structure <- rvinecopulib::cvine_structure(c(3, 1, 4, 2))
truthVine <- copulaGaussianRvineFromCor(Rtruth, structure)
truthMargins <- list(copulaMarginStudent(.48, 3.4),
  copulaMarginLaplace(.34),
  copulaMarginCovariateLognormal(log(72), .30),
  copulaMarginCovariateNormal(91, 16))
E <- copulaMarginsQuantile(rvinecopulib::rvinecop(n, truthVine), truthMargins)
E[, 1L] <- E[, 1L] + .10
E[, 2L] <- E[, 2L] - .07
w <- rep(1 / n, n)
startMargins <- list(copulaMarginStudent(.40, 7),
  copulaMarginLaplace(.42),
  copulaMarginCovariateLognormal(log(68), .38),
  copulaMarginCovariateNormal(88, 20))
startVine <- copulaGaussianRvineFromCor(diag(d), structure)

profileTime <- system.time(profile <-
  copulaMaximiseGaussianFremMultiArbitrary(E, w, startMargins, startVine,
    d, dEta, maxit = 350L))[["elapsed"]]
stopifnot(!is.null(profile), identical(profile$backend,
  "gaussian-copula-frem-multi-arbitrary-profile"))
profileLiteral <- sum(w * copulaGaussianFremLogPrior(
  sweep(E, 2L, profile$delta, "-"), profile$vine, profile$margins,
  dEta, "joint"))
stopifnot(abs(profile$value - profileLiteral) < 1e-8)

## The ordinary eta-intercept design is an equivalent parameterization.
X <- matrix(1, n, dEta)
locMap <- matrix(0, dEta, d); diag(locMap[, seq_len(dEta)]) <- 1
designed <- copulaMaximiseGaussianFremMultiArbitrary(
  E, w, startMargins, startVine, d, dEta, maxit = 350L,
  X = X, locMap = locMap, beta0 = numeric(dEta))
stopifnot(!is.null(designed),
  max(abs(designed$beta - profile$delta[seq_len(dEta)])) < 1e-8,
  abs(designed$value - profile$value) < 1e-8)

## The public dispatcher must choose the same exact profile automatically.
automatic <- copulaMaximiseGaussianFrem(E, w, startMargins, startVine,
  d, dEta, maxit = 350L, diagnostics = TRUE)
stopifnot(identical(automatic$backend,
  "gaussian-copula-frem-multi-arbitrary-profile"),
  abs(automatic$value - profile$value) < 1e-8)

## Compare with the former literal all-parameter fallback on the same Q.
oldDisable <- getOption("saemix.copula.disablePartialCollapse")
oldEcm <- getOption("saemix.copula.disableMarginalEcm")
options(saemix.copula.disablePartialCollapse = TRUE,
  saemix.copula.disableMarginalEcm = TRUE)
literalTime <- system.time(literal <- copulaMaximiseGaussianFrem(
  E, w, startMargins, startVine, d, dEta, maxit = 800L,
  diagnostics = TRUE))[["elapsed"]]
options(saemix.copula.disablePartialCollapse = oldDisable,
  saemix.copula.disableMarginalEcm = oldEcm)
stopifnot(profile$value + 3e-6 >= literal$value)
cat(sprintf(paste0("Two-curved profile %.3fs Q %.10f; literal %.3fs Q %.10f; ",
  "deltaQ %.3g\n"), profileTime, profile$value, literalTime,
  literal$value, profile$value - literal$value))

## Three curved coordinates exercise a nontrivial 3x3 internal correlation
## block; the sole Gaussian coordinate must still be reconstructed exactly.
set.seed(260828L)
truth3 <- list(copulaMarginStudent(.48, 3.4),
  copulaMarginLaplace(.34), copulaMarginCovariateGamma(1.6, 44),
  copulaMarginCovariateNormal(91, 16))
E3 <- copulaMarginsQuantile(
  rvinecopulib::rvinecop(n, truthVine), truth3)
E3[, 1L] <- E3[, 1L] + .10; E3[, 2L] <- E3[, 2L] - .07
start3 <- list(copulaMarginStudent(.40, 7),
  copulaMarginLaplace(.42), copulaMarginCovariateGamma(2.1, 35),
  copulaMarginCovariateNormal(88, 20))
profile3 <- copulaMaximiseGaussianFremMultiArbitrary(
  E3, w, start3, startVine, d, dEta, maxit = 400L)
stopifnot(!is.null(profile3), identical(profile3$backend,
  "gaussian-copula-frem-multi-arbitrary-profile"))
literal3 <- sum(w * copulaGaussianFremLogPrior(
  sweep(E3, 2L, profile3$delta, "-"), profile3$vine,
  profile3$margins, dEta, "joint"))
stopifnot(abs(profile3$value - literal3) < 1e-8)
cat(sprintf("Three-curved profile Q %.10f; literal identity %.3g\n",
  profile3$value, profile3$value - literal3))
