## End-to-end score-SA fit with two dependent categorical FREM covariates.
source("copula/tests/helper-load.R")
setwd(file.path(copulaTestRepo, "copula"))
source("R/simpleModels.R")

set.seed(280838L)
N <- 16L
R <- matrix(c(
  1, .18, .42, -.22,
  .18, 1, -.16, .28,
  .42, -.16, 1, .30,
  -.22, .28, .30, 1), 4L, byrow = TRUE)
stopifnot(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) > 0)
structure <- rvinecopulib::cvine_structure(c(3, 1, 2, 4))
margins <- list(copulaMarginNormal(.25), copulaMarginNormal(.30),
  copulaMarginBernoulli(.35),
  copulaMarginOrdinal(c(.25, .50, .25), labels = 0:2))
vine <- copulaVineForMargins(
  copulaGaussianRvineFromCor(R, structure), margins)
latent <- matrix(rnorm(N * 4L), ncol = 4L) %*% chol(R)
eta <- sweep(latent[, 1:2, drop = FALSE], 2L, c(.25, .30), "*")
sex <- as.numeric(latent[, 3L] > qnorm(.65))
stage <- as.numeric(cut(latent[, 4L],
  c(-Inf, qnorm(.25), qnorm(.75), Inf), labels = FALSE)) - 1L

definition <- MODELS$iv1
psi <- sweep(exp(eta), 2L, definition$true, "*")
nt <- length(definition$times)
data <- data.frame(id = rep(seq_len(N), each = nt),
  dose = definition$dose, time = rep(definition$times, N))
prediction <- definition$f(psi, data$id, cbind(data$dose, data$time))
data$y <- pmax(prediction * (1 + .09 * rnorm(nrow(data))), 1e-7)
conditioning <- cbind(SEX = sex, STAGE = stage)
conditioning[c(4L, 13L), 1L] <- NA_real_
conditioning[c(7L, 13L), 2L] <- NA_real_

copulaSet(vine, margins = margins,
  conditioning = list(values = conditioning,
    variableName = c("SEX", "STAGE")),
  warmStartOnActivate = FALSE, familySet = NULL, refitEvery = 1L,
  guard = FALSE, populationAlgorithm = "score-sa", scoreScale = .006,
  scoreBurn = 5L)
control <- list(seed = 280839L, save = FALSE, save.graphs = FALSE,
  print = FALSE, displayProgress = FALSE, warnings = FALSE,
  nbiter.saemix = c(14L, 8L), nbiter.mcmc = c(1L, 1L, 1L, 0L),
  ll.is = FALSE, fim = FALSE, map = FALSE)
fit <- saemix(saemixModelFor(definition), saemixDataFor(data), control)
state <- copulaGet(fit)

stopifnot(inherits(fit, "SaemixObject"),
  identical(state$lastJoint$backend, "gaussian-copula-frem-score-sa"),
  sum(vapply(state$margins, function(m)
    identical(m$type, "discrete"), logical(1))) == 2L,
  grepl("fixed-support",
    state$lastJoint$scoreTheory$categoricalAugmentation),
  grepl("rectangle-dependent random walks disabled",
    state$lastJoint$scoreTheory$individualKernel),
  isFALSE(state$lastJoint$scoreTheory$rectangleProbabilityUsedInScore),
  identical(state$lastJoint$scoreMethod, "global-centered-difference"),
  state$lastJoint$projectionCount == 0L,
  state$lastJoint$postFreezeBacktrackCount == 0L,
  state$lastJoint$postFreezeNoMoveCount == 0L,
  isTRUE(state$lastJoint$scoreTheory$proposalScaleFrozen),
  isTRUE(state$lastJoint$scoreTheory$runtimeConditionsObserved),
  isTRUE(state$lastJoint$usedAverage),
  all(is.finite(fit@results@fixed.effects)))

cat("multiple-categorical Gaussian-copula FREM end-to-end check passed\n")
