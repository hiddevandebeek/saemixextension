source("copula/tests/helper-load.R")

set.seed(311001L)
n <- 70L; truth <- c(V = 20, CL = 3.5)
R <- matrix(c(1, .3, .3, 1), 2L)
vine <- copulaGaussianRvineFromCor(R,
  rvinecopulib::cvine_structure(c(2, 1)))
marginsTruth <- list(copulaNaturalMarginLognormal(.22),
  copulaNaturalMarginGamma(2.5))
u <- rvinecopulib::rvinecop(n, vine)
typical <- matrix(rep(truth, each = n), nrow = n)
psi <- copulaNaturalMarginsQuantile(u, typical, marginsTruth)
times <- c(.1, .5, 1, 2, 4, 8, 16, 24)
data <- data.frame(id = rep(seq_len(n), each = length(times)),
  dose = 100, time = rep(times, n))
pk <- function(psi, id, xidep)
  xidep[, 1L] / psi[id, 1L] *
    exp(-psi[id, 2L] / psi[id, 1L] * xidep[, 2L])
prediction <- pk(psi, data$id, cbind(data$dose, data$time))
data$y <- pmax(prediction * (1 + .12 * rnorm(nrow(data))), 1e-8)

model <- saemixModel(model = pk, modeltype = "structural", description = "",
  psi0 = matrix(c(19, 3.2), nrow = 1L,
    dimnames = list(NULL, c("V", "CL"))),
  transform.par = c(1, 1), covariance.model = matrix(1, 2L, 2L),
  omega.init = diag(c(.25^2, .5^2)), error.model = "proportional",
  error.init = c(0, .12), verbose = FALSE)
dataset <- saemixData(name.data = data, header = TRUE, name.group = "id",
  name.predictors = c("dose", "time"), name.response = "y", verbose = FALSE)
population <- gaussianCopulaFrem(parameterMargins = list(
  copulaNaturalMarginLognormal(.3), copulaNaturalMarginGamma(4)),
  correlation = matrix(c(1, .1, .1, 1), 2L))
control <- list(seed = 311002L, save = FALSE, save.graphs = FALSE,
  print = FALSE, displayProgress = FALSE, warnings = FALSE,
  nbiter.saemix = c(100L, 60L), nbiter.mcmc = c(2L, 2L, 2L, 0L),
  ll.is = FALSE, fim = FALSE, map = FALSE)

standardPopulation <- gaussianCopulaFrem(etaSd = c(.25, .5),
  correlation = matrix(c(1, .1, .1, 1), 2L))
standardControl <- control; standardControl$nbiter.saemix <- c(60L, 40L)
standard <- saemix(model, dataset, standardControl,
  population = standardPopulation)
selection <- suppressWarnings(copulaSelectParameterMargins(standard,
  supports = c("positive", "positive"),
  candidates = list("lognormal", c("lognormal", "gamma")),
  posteriorDraws = 150L, max.iter = 150L, seed = 311010L,
  optimizerMaxit = 80L, minimumEssFraction = 1e-4,
  minimumPosteriorEssFraction = 1e-3, maximumMcse = 5))
stopifnot(inherits(selection, "saemixNaturalMarginSelection"),
  all(c("lognormal/lognormal", "lognormal/gamma") %in%
    selection$table$families),
  inherits(selection$population, "saemixPopulation"),
  all(is.finite(selection$table$validation_delta_loglik)))

fit <- saemix(model, dataset, control, population = population)
state <- copulaGet(fit)
stopifnot(identical(state$populationScale, "parameter"),
  identical(state$lastJoint$scoreMethod, "hybrid-fixed-reference-path-score"),
  isTRUE(state$lastJoint$scoreTheory$runtimeConditionsObserved),
  state$lastJoint$postFreezeProjectionCount == 0L,
  state$lastJoint$postFreezeBacktrackCount == 0L,
  state$lastJoint$postFreezeNoMoveCount == 0L,
  state$lastJoint$metricUpdateCount > 0L,
  grepl("information inverse", state$lastJoint$scoreTheory$scoreMetric),
  isTRUE(state$lastJoint$scoreTheory$metricFrozen),
  all(is.finite(fit@results@fixed.effects)),
  all(is.finite(state$sd)),
  identical(vapply(state$margins, `[[`, character(1), "name"),
    c("lognormal", "gamma")))

fit@options$nmc.is <- 500L
fit <- suppressWarnings(llisCopula.saemix(fit, defensive = .2,
  batch = 100L, seed = 311003L))
likelihood <- attr(fit, "saemix.copula.likelihood")
stopifnot(is.finite(fit@results@ll.is),
  is.finite(likelihood$se_loglik_total), likelihood$draws_used == 500L)

cat("natural parameter score-SA end-to-end check passed\n")
