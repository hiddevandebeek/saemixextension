## Automatic family selection uses a shared posterior measure and automatically
## dispatches algebraically Gaussianizable candidates to the MVN collapse.
suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})

set.seed(82741)
n <- 700L; dEta <- 1L
R <- matrix(c(1, .55, .55, 1), 2L, 2L)
vine <- copulaGaussianRvineFromCor(R,
  rvinecopulib::cvine_structure(c(2, 1)))
truth <- list(copulaMarginNormal(.30),
  copulaMarginCovariateLognormal(log(70), .55))
E <- copulaMarginsQuantile(rvinecopulib::rvinecop(n, vine), truth)
w <- rep(1 / n, n)
start <- list(copulaMarginNormal(.4),
  copulaMarginCovariateNormal(mean(E[, 2]), sd(E[, 2])))
selection <- copulaSelectGaussianFremMargins(E, w, vine, start,
  dEta = dEta, nSubjects = n, etaFamilies = "normal",
  covariateFamilies = c("normal", "lognormal"), cycles = 2L, maxit = 30L)
stopifnot(inherits(selection, "saemixGaussianFremMarginSelection"),
  identical(unname(selection$selectedFamilies), c("normal", "lognormal")),
  identical(selection$fit$backend, "gaussian-copula-frem-mvnormal-em"),
  all(selection$ranking$backend == "gaussian-copula-frem-mvnormal-em"),
  isTRUE(selection$commonMeasure))

candidates <- copulaGaussianFremMarginCandidates(
  cbind(rnorm(100), rgamma(100, 5, scale = 2)), rep(.01, 100), 1L)
stopifnot(identical(names(candidates[[1]]), c("normal", "student", "laplace")),
  identical(names(candidates[[2]]), c("normal", "lognormal", "gamma", "weibull")))

## The positive-margin fitter must move away from its Weibull starting shape
## even when an L-BFGS trial encounters a non-finite tail value.
set.seed(82742)
weibullData <- rweibull(180L, shape = 2.6, scale = 80)
weibullFit <- copulaFitCovariateMargin(weibullData, "weibull")
stopifnot(abs(weibullFit$parameters["shape"] - 2) > .05,
  abs(weibullFit$parameters["shape"] - 2.6) < .6)
cat("automatic Gaussian-copula FREM margin-selection checks passed\n")
