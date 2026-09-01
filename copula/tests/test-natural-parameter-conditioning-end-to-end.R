suppressPackageStartupMessages({
  library(devtools)
  devtools::load_all("C:/package/saemix-copula", quiet = TRUE)
  library(rvinecopulib)
})

set.seed(902101)
n <- 36L; times <- c(.25, 1, 4, 12)
R <- matrix(c(1, .2, .1, .2, 1, .55, .1, .55, 1), 3L, 3L)
vine <- copulaGaussianRvineFromCor(R, rvinecopulib::cvine_structure(c(3, 1, 2)))
parameterMargins <- list(copulaNaturalMarginLognormal(.22),
  copulaNaturalMarginGamma(2.5))
covariateMargin <- copulaMarginCovariateGamma(2, 20)
z <- matrix(rnorm(n * 3L), n, 3L) %*% chol(R)
typical <- matrix(rep(c(10, 3), each = n), n)
psi <- copulaNaturalMarginsQuantile(pnorm(z[, 1:2]), typical, parameterMargins)
crp <- covariateMargin$quantile(pnorm(z[, 3]), covariateMargin$parameters)
pk <- function(psi, id, xidep) xidep[, 1] / psi[id, 1] *
  exp(-psi[id, 2] / psi[id, 1] * xidep[, 2])
data <- data.frame(id = rep(seq_len(n), each = length(times)), dose = 100,
  time = rep(times, n))
data$y <- pmax(pk(psi, data$id, cbind(data$dose, data$time)) *
  (1 + .12 * rnorm(nrow(data))), 1e-8)
model <- saemixModel(model = pk, modeltype = "structural",
  psi0 = matrix(c(9.5, 2.8), 1L, dimnames = list(NULL, c("V", "CL"))),
  transform.par = c(1, 1), covariance.model = matrix(1, 2L, 2L),
  omega.init = diag(c(.22^2, .4^2)), error.model = "proportional",
  error.init = c(0, .12), verbose = FALSE)
sxdata <- saemixData(name.data = data, header = TRUE, name.group = "id",
  name.predictors = c("dose", "time"), name.response = "y", verbose = FALSE)
population <- gaussianCopulaFrem(parameterMargins = parameterMargins,
  covariates = matrix(crp, ncol = 1, dimnames = list(NULL, "CRP")),
  covariateMargins = list(covariateMargin), correlation = R,
  structure = vine$structure, scoreBurn = 20L,
  gainScale = .12, gainPower = .8)
control <- list(seed = 902102, save = FALSE, save.graphs = FALSE, print = FALSE,
  displayProgress = FALSE, warnings = FALSE, nbiter.saemix = c(80, 40),
  nbiter.mcmc = c(2, 1, 1, 0), ll.is = FALSE, fim = FALSE, map = FALSE)
fit <- saemix(model, sxdata, control, population = population)
state <- copulaGet(fit)
stopifnot(identical(state$populationScale, "parameter"),
  state$dEta == 2L, state$dConditioning == 1L,
  state$lastJoint$postFreezeProjectionCount == 0L,
  state$lastJoint$postFreezeBacktrackCount == 0L,
  state$lastJoint$postFreezeNoMoveCount == 0L)
fit@options$nmc.is <- 300L
fit <- suppressWarnings(llisCopula.saemix(fit, defensive = 0, batch = 50L,
  seed = 902103))
stopifnot(is.finite(fit@results@ll.is),
  is.finite(attr(fit, "saemix.copula.likelihood")$se_loglik_total))

standardPopulation <- copulaPopulation(vine,
  margins = list(copulaMarginNormal(.22), copulaMarginNormal(.4),
    covariateMargin), scale = "transformed-additive",
  conditioning = list(values = matrix(crp, ncol = 1,
    dimnames = list(NULL, "CRP"))), scoreBurn = 20L,
  scoreGainScale = .12, scoreGainPower = .8)
standardControl <- control; standardControl$seed <- 902104
standard <- saemix(model, sxdata, standardControl, population = standardPopulation)
selection <- suppressWarnings(copulaSelectParameterMargins(standard,
  supports = c("positive", "positive"),
  candidates = list("lognormal", "gamma"), posteriorDraws = 80L,
  max.iter = 120L, optimizerMaxit = 40L, minimumEssFraction = 0,
  minimumPosteriorEssFraction = 0, maximumMcse = Inf, seed = 902105))
selectedArguments <- selection$population$arguments
stopifnot(identical(selectedArguments$populationScale, "parameter"),
  ncol(selectedArguments$conditioning$values) == 1L,
  length(selectedArguments$margins) == 3L,
  inherits(selectedArguments$margins[[2L]], "saemix_natural_parameter_margin"),
  inherits(selectedArguments$margins[[3L]], "saemix_copula_margin"))
cat("natural-parameter conditioning end-to-end test passed\n")
