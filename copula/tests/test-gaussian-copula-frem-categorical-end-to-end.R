## Small nonlinear SAEM fit with an observed/missing binary FREM covariate.
source("copula/tests/helper-load.R")
setwd(file.path(copulaTestRepo, "copula"))
source("R/simpleModels.R")

set.seed(280834L)
N <- 14L; rho <- .55; probability <- .35
R <- matrix(c(1, .22, rho, .22, 1, -.18, rho, -.18, 1), 3L)
structure <- rvinecopulib::cvine_structure(c(3, 1, 2))
baseVine <- copulaGaussianRvineFromCor(R, structure)
margins <- list(copulaMarginNormal(.25), copulaMarginNormal(.30),
  copulaMarginBernoulli(probability))
vine <- copulaVineForMargins(baseVine, margins)
latent <- matrix(rnorm(N * 3L), ncol = 3L) %*% chol(R)
eta <- sweep(latent[, 1:2, drop = FALSE], 2L, c(.25, .30), "*")
sex <- as.numeric(latent[, 3L] > qnorm(1 - probability))

definition <- MODELS$iv1
psi <- sweep(exp(eta), 2L, definition$true, "*")
nt <- length(definition$times)
data <- data.frame(id = rep(seq_len(N), each = nt),
  dose = definition$dose, time = rep(definition$times, N))
prediction <- definition$f(psi, data$id, cbind(data$dose, data$time))
data$y <- pmax(prediction * (1 + .09 * rnorm(nrow(data))), 1e-7)
conditioning <- matrix(sex, ncol = 1L, dimnames = list(NULL, "SEX"))
conditioning[c(4L, 11L), 1L] <- NA_real_

control <- list(seed = 280835L, save = FALSE, save.graphs = FALSE,
  print = FALSE, displayProgress = FALSE, warnings = FALSE,
  nbiter.saemix = c(16L, 8L), nbiter.mcmc = c(1L, 1L, 1L, 0L),
  ll.is = FALSE, fim = FALSE, map = FALSE)
## The score route differentiates the exact observed category rectangle while
## eta remains the only missing state. Missing categories are integrated out.
copulaSet(vine, margins = margins,
  conditioning = list(values = conditioning, variableName = "SEX"),
  warmStartOnActivate = FALSE, familySet = NULL, refitEvery = 1L,
  guard = FALSE, populationAlgorithm = "score-sa", scoreScale = .01,
  scoreBurn = 5L)
scoreFit <- saemix(saemixModelFor(definition), saemixDataFor(data), control)
scoreState <- copulaGet(scoreFit)
stopifnot(identical(scoreState$lastJoint$backend,
    "gaussian-copula-frem-score-sa"),
  identical(scoreState$margins[[3L]]$type, "discrete"),
  identical(scoreState$lastJoint$scoreMethod,
    "global-centered-difference"),
  scoreState$lastJoint$projectionCount == 0L,
  isTRUE(scoreState$lastJoint$scoreTheory$proposalScaleFrozen),
  isTRUE(scoreState$lastJoint$scoreTheory$runtimeConditionsObserved),
  scoreState$lastJoint$postFreezeBacktrackCount == 0L,
  scoreState$lastJoint$postFreezeNoMoveCount == 0L,
  isTRUE(scoreState$lastJoint$usedAverage),
  all(is.finite(scoreFit@results@fixed.effects)))

scoreFit@options$nmc.is <- 150L
scored <- suppressWarnings(llisCopula.saemix(scoreFit, defensive = 0,
  batch = 50L, seed = 280836L))
diagnostic <- attr(scored, "saemix.copula.likelihood")
stopifnot(is.finite(scored@results@ll.is),
  is.finite(diagnostic$se_loglik_total), diagnostic$pit_clipped == 0L)
cat("categorical Gaussian-copula FREM end-to-end check passed\n")
