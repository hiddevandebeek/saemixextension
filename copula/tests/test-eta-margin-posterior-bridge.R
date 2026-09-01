source("copula/tests/helper-load.R")

log_add_exp <- function(a, b) {
  maximum <- pmax(a, b)
  maximum + log(exp(a - maximum) + exp(b - maximum))
}
log_normal_laplace <- function(y, etaSd, residualSd) {
  scale <- etaSd / sqrt(2)
  common <- -log(2 * scale) + residualSd^2 / (2 * scale^2)
  common + log_add_exp(
    -y / scale + pnorm(y / residualSd - residualSd / scale, log.p = TRUE),
    y / scale + pnorm(-y / residualSd - residualSd / scale, log.p = TRUE))
}

set.seed(930101L)
n <- 100L; samples <- 1500L; etaSd <- .45; residualSd <- .25
scale <- etaSd / sqrt(2)
eta <- sign(runif(n) - .5) * rexp(n, 1 / scale)
y <- eta + rnorm(n, 0, residualSd)
incumbentSd <- sqrt(max(mean(y^2) - residualSd^2, .02^2))
posteriorVariance <- 1 / (1 / incumbentSd^2 + 1 / residualSd^2)
posteriorMean <- posteriorVariance * y / residualSd^2
draw_pool <- function(seed) {
  set.seed(seed)
  array(matrix(rnorm(n * samples,
    mean = rep(posteriorMean, samples), sd = sqrt(posteriorVariance)),
    nrow = n, ncol = samples), c(n, 1L, samples))
}
training <- draw_pool(930102L)
validation <- draw_pool(930103L)
conditioning <- matrix(rnorm(n), ncol = 1L)
vine <- copulaGaussianRvineFromCor(diag(2L),
  rvinecopulib::cvine_structure(1:2))
state <- list(d = 2L, dEta = 1L, dConditioning = 1L,
  conditioning = conditioning, conditioningName = "COV",
  margins = list(copulaMarginNormal(incumbentSd),
    copulaMarginCovariateNormal(0, 1)), vine = vine)

fits <- lapply(c("normal", "laplace"), function(family)
  copulaFitPosteriorEtaCandidate(family, state, training, validation, 60L))
for (fit in fits) {
  candidateSd <- copulaMarginScales(fit$margins)[1L]
  candidateLogLik <- if (fit$families == "normal")
    sum(dnorm(y, 0, sqrt(candidateSd^2 + residualSd^2), log = TRUE)) else
    sum(log_normal_laplace(y, candidateSd, residualSd))
  incumbentLogLik <- sum(dnorm(y, 0,
    sqrt(incumbentSd^2 + residualSd^2), log = TRUE))
  directDelta <- candidateLogLik - incumbentLogLik
  stopifnot(abs(fit$validation$deltaLogLik - directDelta) <
    3 * fit$validation$mcse + .03)
}

cat("posterior eta-margin bridge identity checks passed\n")
