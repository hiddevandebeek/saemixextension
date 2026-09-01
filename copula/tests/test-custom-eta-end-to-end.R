source("copula/tests/helper-load.R")
suppressPackageStartupMessages(library(rvinecopulib))

asymmetric_margin <- function(sd) {
  p <- .25
  rawScale <- function(par)
    par["sd"] * p * (1 - p) / sqrt(1 - 2 * p + 2 * p^2)
  center <- function(par) {
    b <- rawScale(par); b * (1 - 2 * p) / (p * (1 - p))
  }
  logDensity <- function(x, par) {
    b <- rawScale(par); u <- x / b
    log(p * (1 - p) / b) - u * (p - as.numeric(u < 0))
  }
  cdf <- function(x, par) {
    b <- rawScale(par)
    ifelse(x < 0, p * exp((1 - p) * x / b),
      1 - (1 - p) * exp(-p * x / b))
  }
  quantile <- function(u, par) {
    b <- rawScale(par)
    ifelse(u < p, b / (1 - p) * log(u / p),
      -b / p * log((1 - u) / (1 - p)))
  }
  copulaMarginCenteredCustom("user-asymmetric-laplace", c(sd = sd),
    1e-6, 10, logDensity, cdf, quantile, center,
    function(par) par["sd"], support_fixed = TRUE,
    set_scale = function(par, value) { par["sd"] <- value; par },
    roles = "scale")
}

set.seed(593101L)
N <- 60L; times <- c(.5, 1, 2, 4, 8)
truthR <- matrix(c(1, .60, .60, 1), 2L)
structure <- cvine_structure(c(2, 1))
truthVine <- copulaGaussianRvineFromCor(truthR, structure)
truthMargins <- list(asymmetric_margin(.35),
  copulaMarginCovariateNormal(0, 1))
joint <- copulaMarginsQuantile(rvinecop(N, truthVine), truthMargins)
psi <- 3 * exp(joint[, 1L])
data <- data.frame(id = rep(seq_len(N), each = length(times)),
  time = rep(times, N))
prediction <- psi[data$id] * exp(-.22 * data$time)
data$y <- pmax(prediction * (1 + .10 * rnorm(nrow(data))), 1e-8)
conditioning <- matrix(joint[, 2L], ncol = 1L,
  dimnames = list(NULL, "COV"))

modelFunction <- function(psi, id, xidep)
  psi[id, 1L] * exp(-.22 * xidep[, 1L])
model <- saemixModel(model = modelFunction, modeltype = "structural",
  psi0 = matrix(2.9, 1L, dimnames = list(NULL, "A")),
  transform.par = 1, covariance.model = matrix(1, 1L, 1L),
  omega.init = matrix(.40^2, 1L), error.model = "proportional",
  error.init = .12, verbose = FALSE)
sxdata <- saemixData(name.data = data, header = TRUE, name.group = "id",
  name.predictors = "time", name.response = "y", verbose = FALSE)
R0 <- matrix(c(1, .40, .40, 1), 2L)
population <- copulaPopulation(
  copulaGaussianRvineFromCor(R0, structure),
  margins = list(asymmetric_margin(.42),
    copulaFitCovariateMargin(conditioning[, 1L], "normal")),
  scale = "transformed-additive",
  conditioning = list(values = conditioning, variableName = "COV"),
  populationAlgorithm = "score-sa", scoreScale = .10,
  scoreGainScale = .2, scoreGainPower = .8, scoreGainOffset = 30,
  scoreBurn = 50L, scoreFiniteDifference = 1e-4,
  scoreProjection = 24)
control <- list(seed = 593102L, save = FALSE, save.graphs = FALSE,
  print = FALSE, displayProgress = FALSE, warnings = FALSE,
  nbiter.saemix = c(500L, 1000L), nbiter.mcmc = c(2L, 2L, 2L, 0L),
  ll.is = FALSE, fim = FALSE, map = FALSE)
fit <- saemix(model, sxdata, control, population = population)
state <- copulaGet(fit)
rho <- copulaGaussianRvineCor(state$vine, 2L)[1L, 2L]
cat(sprintf(paste0("custom eta diagnostic: A=%.3f, sd=%.3f, rho=%.3f, ",
  "scoreAverage=%.4g\n"), fit@results@fixed.effects[1L],
  state$margins[[1L]]$parameters[["sd"]], rho,
  state$lastJoint$scoreAverageMax))
stopifnot(identical(state$margins[[1L]]$name, "user-asymmetric-laplace"),
  isTRUE(state$lastJoint$scoreTheory$runtimeConditionsObserved),
  state$lastJoint$postFreezeProjectionCount == 0L,
  state$lastJoint$postFreezeBacktrackCount == 0L,
  state$lastJoint$postFreezeNoMoveCount == 0L,
  is.finite(state$margins[[1L]]$parameters[["sd"]]),
  abs(state$margins[[1L]]$parameters[["sd"]] - .35) < .16,
  abs(rho - .60) < .25,
  abs(fit@results@fixed.effects[1L] - 3) < .45)
cat(sprintf(paste0("custom eta end-to-end fit passed: A=%.3f, sd=%.3f, ",
  "rho=%.3f\n"), fit@results@fixed.effects[1L],
  state$margins[[1L]]$parameters[["sd"]], rho))
