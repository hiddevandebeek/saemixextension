## Exact envelope-score acceleration for supported one-curved-margin profiles.
## Every accelerated result is compared with the identical profiled objective
## optimized using optim()'s numerical gradient.
suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})

set.seed(260826L)
n <- 500L; d <- 4L; dEta <- 2L
R <- matrix(c(1, .35, .48, -.20,
              .35, 1, -.16, .30,
              .48, -.16, 1, .22,
              -.20, .30, .22, 1), d, d, byrow = TRUE)
R <- stats::cov2cor(copulaGaussianFremEnsurePd(R, 1e-6))
structure <- rvinecopulib::cvine_structure(c(3, 1, 2, 4))
vine <- copulaGaussianRvineFromCor(R, structure)
startVine <- copulaGaussianRvineFromCor(diag(d), structure)
w <- rep(1 / n, n)

cases <- list(
  gamma = list(
    truth = list(copulaMarginNormal(.42), copulaMarginNormal(.31),
      copulaMarginCovariateGamma(1.4, 48),
      copulaMarginCovariateNormal(92, 17)),
    start = list(copulaMarginNormal(.38), copulaMarginNormal(.35),
      copulaMarginCovariateGamma(2.2, 35),
      copulaMarginCovariateNormal(88, 20))),
  weibull = list(
    truth = list(copulaMarginNormal(.42), copulaMarginNormal(.31),
      copulaMarginCovariateWeibull(.90, 80),
      copulaMarginCovariateNormal(92, 17)),
    start = list(copulaMarginNormal(.38), copulaMarginNormal(.35),
      copulaMarginCovariateWeibull(1.00, 75),
      copulaMarginCovariateNormal(88, 20))),
  laplace = list(
    truth = list(copulaMarginLaplace(.50), copulaMarginNormal(.31),
      copulaMarginCovariateLognormal(log(72), .30),
      copulaMarginCovariateNormal(92, 17)),
    start = list(copulaMarginLaplace(.40), copulaMarginNormal(.35),
      copulaMarginCovariateLognormal(log(68), .38),
      copulaMarginCovariateNormal(88, 20))),
  student = list(
    truth = list(copulaMarginStudent(.50, 3.3), copulaMarginNormal(.31),
      copulaMarginCovariateLognormal(log(72), .30),
      copulaMarginCovariateNormal(92, 17)),
    start = list(copulaMarginStudent(.40, 7), copulaMarginNormal(.35),
      copulaMarginCovariateLognormal(log(68), .38),
      copulaMarginCovariateNormal(88, 20))))

oldDisable <- getOption("saemix.copula.disableEnvelopeGradient")
on.exit(options(saemix.copula.disableEnvelopeGradient = oldDisable), add = TRUE)
result <- lapply(names(cases), function(family) {
  cat("Checking", family, "...\n")
  case <- cases[[family]]
  E <- copulaMarginsQuantile(rvinecopulib::rvinecop(n, vine), case$truth)
  if (family %in% c("laplace", "student")) E[, 1L] <- E[, 1L] + .11
  options(saemix.copula.disableEnvelopeGradient = FALSE)
  acceleratedTime <- system.time(accelerated <-
    copulaMaximiseGaussianFremOneArbitrary(E, w, case$start, startVine,
      d, dEta, maxit = 300L))[["elapsed"]]
  options(saemix.copula.disableEnvelopeGradient = TRUE)
  numericalTime <- system.time(numerical <-
    copulaMaximiseGaussianFremOneArbitrary(E, w, case$start, startVine,
      d, dEta, maxit = 300L))[["elapsed"]]
  stopifnot(!is.null(accelerated), !is.null(numerical))
  literal <- sum(w * copulaGaussianFremLogPrior(
    sweep(E, 2L, accelerated$delta, "-"), accelerated$vine,
    accelerated$margins, dEta, "joint"))
  parameterDifference <- max(abs(
    accelerated$margins[[if (family %in% c("gamma", "weibull")) 3L else 1L]]$parameters -
    numerical$margins[[if (family %in% c("gamma", "weibull")) 3L else 1L]]$parameters))
  stopifnot(abs(accelerated$value - literal) < 1e-8,
    abs(accelerated$value - numerical$value) < 3e-7,
    parameterDifference < 2e-3,
    max(abs(copulaGaussianRvineCor(accelerated$vine) -
      copulaGaussianRvineCor(numerical$vine))) < 2e-4)
  data.frame(family = family, accelerated_s = acceleratedTime,
    numerical_s = numericalTime, speedup = numericalTime / acceleratedTime,
    objective_difference = accelerated$value - numerical$value,
    parameter_difference = parameterDifference)
})
result <- do.call(rbind, result)
print(result, row.names = FALSE, digits = 6)
