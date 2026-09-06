## Exact all-marginal ECM/profile path for Gaussian-copula FREM.
## The test combines every built-in non-affine margin and a generic beta margin,
## then compares the accelerated result with the literal common-Q engine.
suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})

set.seed(260828L)
n <- 320L; d <- 5L; dEta <- 2L
Rtruth <- matrix(c(
  1, .30, .42, -.16, .22,
  .30, 1, -.20, .38, -.12,
  .42, -.20, 1, .24, .31,
  -.16, .38, .24, 1, .18,
  .22, -.12, .31, .18, 1), d, d, byrow = TRUE)
Rtruth <- stats::cov2cor(copulaGaussianFremEnsurePd(Rtruth, 1e-6))
structure <- rvinecopulib::cvine_structure(c(3, 1, 5, 2, 4))
truthVine <- copulaGaussianRvineFromCor(Rtruth, structure)
betaMargin <- function(a, b) copulaMarginDistribution(
  "beta", c(shape1 = a, shape2 = b), c(.05, .05), c(100, 100),
  centered = FALSE, roles = c("shape", "shape"))
truthMargins <- list(
  copulaMarginStudent(.48, 3.5),
  copulaMarginLaplace(.35),
  copulaMarginCovariateGamma(1.55, 45),
  copulaMarginCovariateWeibull(1.15, 82),
  betaMargin(2.2, 5.5))
startMargins <- list(
  copulaMarginStudent(.39, 7),
  copulaMarginLaplace(.43),
  copulaMarginCovariateGamma(2.1, 35),
  copulaMarginCovariateWeibull(.90, 74),
  betaMargin(1.6, 4.2))
E <- copulaMarginsQuantile(
  rvinecopulib::rvinecop(n, truthVine), truthMargins)
E[, 1L] <- E[, 1L] + .11
E[, 2L] <- E[, 2L] - .08
w <- rep(1 / n, n)
startVine <- copulaGaussianRvineFromCor(diag(d), structure)

## The analytic hyperspherical correlation gradient is an exact derivative of
## the same Gaussian-copula correlation objective.
zProbe <- matrix(stats::rnorm(600L * d), 600L, d)
sProbe <- crossprod(zProbe) / nrow(zProbe)
angles <- copulaGaussianCorrelationAngles(Rtruth)
analytic <- copulaGaussianCorrelationAngleObjective(angles, sProbe)$gradient
numeric <- vapply(seq_along(angles), function(j) {
  h <- 1e-6; plus <- minus <- angles
  plus[j] <- plus[j] + h; minus[j] <- minus[j] - h
  (copulaGaussianCorrelationAngleObjective(plus, sProbe)$value -
    copulaGaussianCorrelationAngleObjective(minus, sProbe)$value) / (2 * h)
}, numeric(1))
stopifnot(max(abs(analytic - numeric)) < 2e-6)

oldSweeps <- getOption("saemix.copula.ecmMaxSweeps")
oldPartial <- getOption("saemix.copula.disablePartialCollapse")
oldEcm <- getOption("saemix.copula.disableMarginalEcm")
on.exit(options(saemix.copula.ecmMaxSweeps = oldSweeps,
  saemix.copula.disablePartialCollapse = oldPartial,
  saemix.copula.disableMarginalEcm = oldEcm), add = TRUE)
options(saemix.copula.ecmMaxSweeps = 20L,
  saemix.copula.disablePartialCollapse = FALSE,
  saemix.copula.disableMarginalEcm = FALSE)

ecmTime <- system.time(ecm <- copulaMaximiseGaussianFrem(
  E, w, startMargins, startVine, d, dEta, maxit = 120L,
  diagnostics = TRUE))[["elapsed"]]
stopifnot(identical(ecm$backend, "gaussian-copula-frem-marginal-ecm"))
ecmLiteral <- sum(w * copulaGaussianFremLogPrior(
  sweep(E, 2L, ecm$delta, "-"), ecm$vine, ecm$margins,
  dEta, "joint"))
stopifnot(abs(ecm$value - ecmLiteral) < 1e-8)
startQ <- sum(w * copulaGaussianFremLogPrior(
  E, startVine, startMargins, dEta, "joint"))
oneSweepTime <- system.time(oneSweep <- copulaMaximiseGaussianFrem(
  E, w, startMargins, startVine, d, dEta, maxit = 40L,
  diagnostics = TRUE, ecmSweeps = 1L))[["elapsed"]]
stopifnot(oneSweep$value + 1e-10 >= startQ,
  identical(oneSweep$ecmSweeps, 1L))

## Ordinary FREM eta-intercept design is the same parameterization and keeps
## the fitted location in beta rather than in the returned residual delta.
X <- matrix(1, n, dEta)
locMap <- matrix(0, dEta, d); diag(locMap[, seq_len(dEta)]) <- 1
designed <- copulaMaximiseGaussianFrem(E, w, startMargins, startVine,
  d, dEta, maxit = 120L, X = X, locMap = locMap,
  beta0 = numeric(dEta), diagnostics = TRUE)
stopifnot(identical(designed$backend,
    "gaussian-copula-frem-marginal-ecm"),
  max(abs(designed$beta - ecm$delta[seq_len(dEta)])) < 2e-6,
  abs(designed$value - ecm$value) < 2e-7)
repeatFit <- copulaMaximiseGaussianFrem(E, w, designed$margins, designed$vine,
  d, dEta, maxit = 120L, X = X, locMap = locMap,
  beta0 = designed$beta, diagnostics = TRUE)
cat(sprintf("Repeat ECM deltaQ %.3g; max beta change %.3g\n",
  repeatFit$value - designed$value,
  max(abs(repeatFit$beta - designed$beta))))
stopifnot(abs(repeatFit$value - designed$value) < 2e-7)

## Retained conditional particles repeat observed covariates. Their exact-value
## grouping must reduce Gamma/Weibull/beta evaluations to one per subject while
## retaining the literal common Q.
index <- rep(seq_len(n), each = 5L)
Erep <- E[index, ]; Erep[, seq_len(dEta)] <-
  Erep[, seq_len(dEta)] + matrix(rnorm(length(index) * dEta, sd = .01),
    ncol = dEta)
wrep <- rep(1 / length(index), length(index))
grouped <- copulaMaximiseGaussianFrem(Erep, wrep, startMargins, startVine,
  d, dEta, maxit = 80L, diagnostics = TRUE)
groupedLiteral <- sum(wrep * copulaGaussianFremLogPrior(
  sweep(Erep, 2L, grouped$delta, "-"), grouped$vine, grouped$margins,
  dEta, "joint"))
stopifnot(all(grouped$marginalEvaluationRows[(dEta + 1L):d] == n),
  abs(grouped$value - groupedLiteral) < 1e-8)

## Disable every collapse/profile path to recover the former literal joint
## optimizer on exactly the same model and objective.
options(saemix.copula.disablePartialCollapse = TRUE,
  saemix.copula.disableMarginalEcm = TRUE)
literalTime <- system.time(literal <- copulaMaximiseGaussianFrem(
  E, w, startMargins, startVine, d, dEta, maxit = 500L,
  diagnostics = TRUE))[["elapsed"]]
stopifnot(ecm$value + 2e-5 >= literal$value)

cat(sprintf(paste0(
  "All-marginal ECM %.3fs Q %.10f (%d sweeps); ",
  "one SAEM sweep %.3fs; literal %.3fs Q %.10f; ",
  "speedup %.2fx; deltaQ %.3g\n"),
  ecmTime, ecm$value, ecm$ecmSweeps, oneSweepTime,
  literalTime, literal$value,
  literalTime / ecmTime, ecm$value - literal$value))
