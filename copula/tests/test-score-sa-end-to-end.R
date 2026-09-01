source("copula/tests/helper-load.R")
setwd(file.path(copulaTestRepo, "copula"))
source("R/simpleModels.R")

set.seed(280838L)
N <- 18L; modelDefinition <- MODELS$iv1
R <- matrix(c(1, .3, .5, .3, 1, -.25, .5, -.25, 1), 3L, 3L)
structure <- rvinecopulib::cvine_structure(c(3, 1, 2))
vine <- copulaGaussianRvineFromCor(R, structure)
margins <- list(copulaMarginNormal(.28), copulaMarginNormal(.34),
  copulaMarginCovariateGamma(1.5, 45))
joint <- copulaMarginsQuantile(rvinecopulib::rvinecop(N, vine), margins)
eta <- joint[, 1:2, drop = FALSE]; covariate <- joint[, 3L]
psi <- sweep(exp(eta), 2L, modelDefinition$true, "*")
nt <- length(modelDefinition$times)
dataset <- data.frame(id = rep(seq_len(N), each = nt),
  dose = modelDefinition$dose, time = rep(modelDefinition$times, N))
prediction <- modelDefinition$f(
  psi, dataset$id, cbind(dataset$dose, dataset$time))
dataset$y <- pmax(prediction * (1 + .08 * rnorm(nrow(dataset))), 1e-6)
conditioning <- matrix(covariate, ncol = 1L,
  dimnames = list(NULL, "COV"))
startMargins <- list(copulaMarginNormal(.25), copulaMarginNormal(.30),
  copulaFitCovariateMargin(covariate, candidates = "gamma"))

copulaSet(vine, margins = startMargins,
  conditioning = list(values = conditioning, variableName = "COV"),
  warmStartOnActivate = FALSE, familySet = NULL, refitEvery = 1L,
  guard = FALSE, populationAlgorithm = "score-sa", scoreScale = .02,
  scoreBurn = 5L)
control <- list(seed = 280839L, save = FALSE, save.graphs = FALSE,
  print = FALSE, displayProgress = FALSE,
  nbiter.saemix = c(15, 15), nbiter.mcmc = c(1, 1, 1, 0),
  warnings = FALSE, ll.is = FALSE, fim = FALSE, map = FALSE)
fit <- saemix(saemixModelFor(modelDefinition), saemixDataFor(dataset), control)
state <- copulaGet(fit)
stopifnot(inherits(fit, "SaemixObject"),
  identical(state$lastJoint$backend,
    "gaussian-copula-frem-score-sa"),
  all(is.finite(fit@results@fixed.effects)),
  all(is.finite(fit@results@respar)),
  all(is.finite(state$sd)),
  state$lastJoint$projectionCount == 0L,
  isTRUE(state$lastJoint$scoreTheory$metricFrozen),
  isTRUE(state$lastJoint$scoreTheory$proposalScaleFrozen),
  isTRUE(state$lastJoint$scoreTheory$blockScheduleFrozen),
  isTRUE(state$lastJoint$scoreTheory$runtimeConditionsObserved),
  identical(state$lastJoint$rwProposalFreezeIteration, 5L),
  !is.null(state$rwProposalFrozen),
  identical(state$rwBlockSizeFrozen, 2L),
  state$lastJoint$postFreezeBacktrackCount == 0L,
  state$lastJoint$postFreezeNoMoveCount == 0L,
  isTRUE(state$lastJoint$usedAverage))

selection <- suppressWarnings(copulaSelectEtaMargins(fit,
  candidates = c("normal", "laplace"), posteriorDraws = 50L,
  max.iter = 50L, optimizerMaxit = 20L,
  minimumEssFraction = 1e-4, maximumMcse = 5))
stopifnot(inherits(selection, "saemixEtaMarginSelection"),
  nrow(selection$table) == 4L,
  all(is.finite(selection$table$validation_delta_loglik)),
  inherits(selection$population, "saemixPopulation"),
  isTRUE(selection$diagnostics$selectionOutsideScoreSa),
  isTRUE(selection$diagnostics$fixedDuringFinalFit),
  is.finite(selection$diagnostics$trainingMcmc$minimumMcmcEssFraction),
  is.finite(selection$diagnostics$validationMcmc$minimumMcmcEssFraction))
refitControl <- control
refitControl$seed <- control$seed + 1L
refitControl$nbiter.saemix <- c(8L, 8L)
refit <- saemix(saemixModelFor(modelDefinition), saemixDataFor(dataset),
  refitControl, population = selection$population)
refitState <- copulaGet(refit)
stopifnot(identical(vapply(refitState$margins[seq_len(refitState$dEta)],
    `[[`, character(1), "name"), as.character(selection$families)),
  identical(refitState$populationAlgorithm, "score-sa"))

cat("score-SA full saemix end-to-end check passed\n")
