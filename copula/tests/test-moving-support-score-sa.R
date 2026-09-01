## Fixed-reference score-SA for continuous margins with parameter-dependent support.
source("copula/tests/helper-load.R")

## Constructor and probability-integral-transform checks.
uniform <- copulaMarginUniformCentered(.25)
gamma <- copulaMarginCenteredGamma(shape = 5, sd = .30)
covUniform <- copulaMarginCovariateUniform(40, 100)
genericUniform <- copulaMarginDistribution("unif",
  parameters = c(min = -2, max = 3), lower = c(-10, .1),
  upper = c(-.1, 10), roles = c("location", "scale"),
  scale = function(par) (par["max"] - par["min"]) / sqrt(12),
  scale_is_sd = TRUE, support_fixed = FALSE)
invisible(lapply(list(uniform, gamma, covUniform, genericUniform),
  copulaMarginValidate))
u <- seq(.01, .99, length.out = 101L)
stopifnot(max(abs(uniform$cdf(uniform$quantile(u, uniform$parameters),
    uniform$parameters) - u)) < 1e-12,
  max(abs(gamma$cdf(gamma$quantile(u, gamma$parameters),
    gamma$parameters) - u)) < 2e-10,
  all(vapply(list(uniform, gamma, covUniform, genericUniform),
    copulaMarginHasMovingSupport, logical(1))))

## Cached E-step prior algebra is exactly the same density as the public joint
## evaluator, including the moving-support Gamma margin.
cacheR <- matrix(c(1, .37, .37, 1), 2L)
cacheVine <- copulaGaussianRvineFromCor(cacheR,
  rvinecopulib::cvine_structure(c(2, 1)))
cacheMargins <- list(copulaMarginNormal(.27), gamma)
cacheEta <- copulaMarginsQuantile(
  rvinecopulib::rvinecop(150L, cacheVine), cacheMargins)
cacheKernel <- copulaGaussianFremContinuousPriorKernel(
  cacheVine, cacheMargins)
cacheLiteral <- -copulaGaussianFremLogPrior(cacheEta, cacheVine,
  cacheMargins, 2L, "joint")
stopifnot(max(abs(cacheKernel$negative(cacheEta) - cacheLiteral)) < 2e-12)

## Independent Fisher oracle for Y|U with X=Q_theta(U). This verifies the
## pathwise response term that a native-X complete score would miss.
y <- .18; sigma <- .22; half <- .48
grid <- (seq_len(200000L) - .5) / 200000L
logWeight <- stats::dnorm(y - half * (2 * grid - 1), sd = sigma, log = TRUE)
weight <- exp(logWeight - max(logWeight)); weight <- weight / sum(weight)
completeScore <- (y - half * (2 * grid - 1)) / sigma^2 * (2 * grid - 1)
posteriorScore <- sum(weight * completeScore)
observed <- function(h) mean(stats::dnorm(y - h * (2 * grid - 1), sd = sigma))
step <- 1e-5
observedScore <- (log(observed(half + step)) - log(observed(half - step))) /
  (2 * step)
stopifnot(abs(posteriorScore - observedScore) < 2e-7)

## Shared hybrid score: only the fixed-reference path coordinates are
## numerical; ordinary margin, dependence and residual coordinates remain
## analytic. Compare the complete internal score with the former global
## centered-difference oracle.
set.seed(300840L)
nOracle <- 70L
oracleR <- matrix(c(1, .32, .32, 1), 2L)
oracleVine <- copulaGaussianRvineFromCor(oracleR,
  rvinecopulib::cvine_structure(c(2, 1)))
oracleMargins <- list(copulaMarginNormal(.24),
  copulaMarginCenteredGamma(3.2, .38))
oracleZ <- matrix(rnorm(nOracle * 2L), ncol = 2L) %*% chol(oracleR)
oracleEta <- cbind(.24 * oracleZ[, 1L],
  oracleMargins[[2L]]$quantile(pnorm(oracleZ[, 2L]),
    oracleMargins[[2L]]$parameters))
oracleReference <- matrix(NA_real_, nOracle, 2L)
oracleReference[, 2L] <- oracleMargins[[2L]]$cdf(oracleEta[, 2L],
  oracleMargins[[2L]]$parameters)
oracleF <- exp(oracleEta[, 1L]) + .7 * exp(oracleEta[, 2L])
oracleY <- oracleF * (1 + .11 * rnorm(nOracle))
oracleResponse <- list(y = oracleY, f = oracleF, etype = rep(1L, nOracle),
  pres = c(0, .12), free = 2L,
  evaluate = function(phiRandom, candidateResidual) {
    f <- exp(phiRandom[, 1L]) + .7 * exp(phiRandom[, 2L])
    sd <- error(f, candidateResidual, rep(1L, length(f)))
    -.5 * ((oracleY - f) / sd)^2 - log(sd) - .5 * log(2 * pi)
  })
oracleResponse$gradient <- function(phiRandom, candidateResidual,
                                    coordinates, step) {
  answer <- matrix(0, nrow(phiRandom), ncol(phiRandom))
  for (coordinate in coordinates) {
    plus <- minus <- phiRandom
    plus[, coordinate] <- plus[, coordinate] + step
    minus[, coordinate] <- minus[, coordinate] - step
    answer[, coordinate] <- (oracleResponse$evaluate(plus,
      candidateResidual) - oracleResponse$evaluate(minus,
        candidateResidual)) / (2 * step)
  }
  answer
}
hybrid <- copulaGaussianFremPopulationScoreStep(oracleEta,
  rep(1 / nOracle, nOracle), oracleMargins, oracleVine, 2L, 2L, .35,
  scoreScale = .002, finiteDifference = 2e-5,
  response = oracleResponse, referenceUniform = oracleReference,
  analyticScore = TRUE)
global <- copulaGaussianFremPopulationScoreStep(oracleEta,
  rep(1 / nOracle, nOracle), oracleMargins, oracleVine, 2L, 2L, .35,
  scoreScale = .002, finiteDifference = 2e-5,
  response = oracleResponse, referenceUniform = oracleReference,
  analyticScore = FALSE)
hybridError <- max(abs(hybrid$score - global$score))
stopifnot(hybridError < 5e-4,
  identical(hybrid$scoreMethod, "hybrid-fixed-reference-path-score"),
  length(hybrid$numericalScoreComponents) < length(hybrid$score))

## End-to-end nonlinear fit with two different moving-support eta margins.
setwd(file.path(copulaTestRepo, "copula"))
source("R/simpleModels.R")
set.seed(300841L)
N <- 18L; R <- matrix(c(1, .28, .28, 1), 2L)
structure <- rvinecopulib::cvine_structure(c(2, 1))
margins <- list(copulaMarginUniformCentered(.25),
  copulaMarginCenteredGamma(5, .30))
vine <- copulaVineForMargins(
  copulaGaussianRvineFromCor(R, structure), margins)
latent <- matrix(rnorm(N * 2L), ncol = 2L) %*% chol(R)
eta <- cbind(margins[[1L]]$quantile(pnorm(latent[, 1L]),
    margins[[1L]]$parameters),
  margins[[2L]]$quantile(pnorm(latent[, 2L]),
    margins[[2L]]$parameters))
definition <- MODELS$iv1
psi <- sweep(exp(eta), 2L, definition$true, "*")
nt <- length(definition$times)
data <- data.frame(id = rep(seq_len(N), each = nt),
  dose = definition$dose, time = rep(definition$times, N))
prediction <- definition$f(psi, data$id, cbind(data$dose, data$time))
data$y <- pmax(prediction * (1 + .08 * rnorm(nrow(data))), 1e-7)

population <- gaussianCopulaFrem(etaMargins = margins,
  correlation = R, structure = structure,
  scoreScale = .004, scoreBurn = 5L)
control <- list(seed = 300842L, save = FALSE, save.graphs = FALSE,
  print = FALSE, displayProgress = FALSE, warnings = FALSE,
  nbiter.saemix = c(16L, 10L), nbiter.mcmc = c(1L, 1L, 1L, 0L),
  ll.is = FALSE, fim = FALSE, map = FALSE)
fit <- saemix(saemixModelFor(definition), saemixDataFor(data), control,
  population = population)
state <- copulaGet(fit)
stopifnot(inherits(fit, "SaemixObject"),
  identical(state$lastJoint$backend, "gaussian-copula-frem-score-sa"),
  grepl("fixed percentile", state$lastJoint$scoreTheory$movingSupportAugmentation),
  identical(state$lastJoint$scoreMethod,
    "hybrid-fixed-reference-path-score"),
  state$lastJoint$projectionCount == 0L,
  state$lastJoint$postFreezeBacktrackCount == 0L,
  state$lastJoint$postFreezeNoMoveCount == 0L,
  isTRUE(state$lastJoint$scoreTheory$runtimeConditionsObserved),
  all(is.finite(fit@results@fixed.effects)))

## Moving-support conditioning margin, including missing covariate values.
set.seed(300843L)
N2 <- 14L
R2 <- matrix(c(1, .15, .45, .15, 1, -.20, .45, -.20, 1), 3L)
structure2 <- rvinecopulib::cvine_structure(c(3, 1, 2))
latent2 <- matrix(rnorm(N2 * 3L), ncol = 3L) %*% chol(R2)
eta2 <- sweep(latent2[, 1:2, drop = FALSE], 2L, c(.25, .30), "*")
covMargin <- copulaMarginCovariateUniform(35, 105)
covariate <- covMargin$quantile(pnorm(latent2[, 3L]), covMargin$parameters)
psi2 <- sweep(exp(eta2), 2L, definition$true, "*")
data2 <- data.frame(id = rep(seq_len(N2), each = nt),
  dose = definition$dose, time = rep(definition$times, N2))
prediction2 <- definition$f(psi2, data2$id, cbind(data2$dose, data2$time))
data2$y <- pmax(prediction2 * (1 + .08 * rnorm(nrow(data2))), 1e-7)
conditioning2 <- matrix(covariate, ncol = 1L,
  dimnames = list(NULL, "BOUNDED_COV"))
conditioning2[c(4L, 11L), 1L] <- NA_real_
population2 <- gaussianCopulaFrem(etaSd = c(.25, .30),
  covariates = conditioning2, covariateMargins = list(covMargin),
  correlation = R2, structure = structure2,
  scoreScale = .003, scoreBurn = 4L)
control2 <- control; control2$seed <- 300844L
control2$nbiter.saemix <- c(12L, 8L)
fit2 <- saemix(saemixModelFor(definition), saemixDataFor(data2), control2,
  population = population2)
state2 <- copulaGet(fit2)
stopifnot(copulaMarginHasMovingSupport(state2$margins[[3L]]),
  anyNA(state2$conditioning),
  grepl("fixed percentile", state2$lastJoint$scoreTheory$movingSupportAugmentation),
  state2$lastJoint$postFreezeBacktrackCount == 0L,
  state2$lastJoint$postFreezeNoMoveCount == 0L,
  all(is.finite(fit2@results@fixed.effects)))

cat(sprintf(paste0("moving-support score-SA checks passed; Fisher error %.3g, ",
  "hybrid score error %.3g\n"),
  abs(posteriorScore - observedScore), hybridError))
