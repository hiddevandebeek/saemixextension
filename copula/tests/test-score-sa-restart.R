## Re-initialisation on expansion of the truncation set (Fort, Moulines,
## Schreck and Vihola 2016, Algorithm 2). A far start inside a deliberately
## small box must trigger at least one expansion with a restart, and the
## restarted recursion must still reach the same point as an unconstrained
## run from the same start.
source("copula/tests/helper-load.R")
suppressPackageStartupMessages(library(rvinecopulib))

set.seed(906201)
n <- 60L; times <- c(.25, 1, 4, 12)
R <- matrix(c(1, .2, .1, .2, 1, .55, .1, .55, 1), 3L, 3L)
vine <- copulaGaussianRvineFromCor(R, rvinecopulib::cvine_structure(c(3, 1, 2)))
z <- matrix(rnorm(n * 3L), n, 3L) %*% chol(R)
psi <- cbind(10 * exp(.22 * z[, 1]), 3 * exp(.30 * z[, 2]))
crp <- 40 + 15 * z[, 3]
pk <- function(psi, id, xidep) xidep[, 1] / psi[id, 1] *
  exp(-psi[id, 2] / psi[id, 1] * xidep[, 2])
data <- data.frame(id = rep(seq_len(n), each = length(times)), dose = 100,
  time = rep(times, n))
data$y <- pk(psi, data$id, cbind(data$dose, data$time)) *
  (1 + .12 * rnorm(nrow(data)))
sxdata <- saemixData(name.data = data, header = TRUE, name.group = "id",
  name.predictors = c("dose", "time"), name.response = "y", verbose = FALSE)
## Far start: typical values doubled, the covariate-margin mean doubled.
model <- saemixModel(model = pk, modeltype = "structural",
  psi0 = matrix(c(20, 6), 1L, dimnames = list(NULL, c("V", "CL"))),
  transform.par = c(1, 1), covariance.model = matrix(1, 2L, 2L),
  omega.init = diag(c(.3^2, .4^2)), error.model = "proportional",
  error.init = c(0, .12), verbose = FALSE)
population <- function(projection) copulaPopulation(vine,
  margins = list(copulaMarginNormal(.3), copulaMarginNormal(.4),
    copulaMarginCovariateNormal(80, 15)), scale = "transformed-additive",
  conditioning = list(values = matrix(crp, ncol = 1,
    dimnames = list(NULL, "CRP"))), scoreBurn = 20L,
  scoreGainPower = .8, scoreProjection = projection)
control <- list(seed = 906202, save = FALSE, save.graphs = FALSE,
  print = FALSE, displayProgress = FALSE, warnings = FALSE,
  nbiter.saemix = c(400, 100), nbiter.mcmc = c(2, 1, 1, 0), ll.is = FALSE,
  fim = FALSE, map = FALSE)

## Box of 30 units: the covariate mean must travel 40, so the box is hit once
## and doubled to 60.
small <- saemix(model, sxdata, control, population = population(30))
smallState <- copulaGet(small)$lastJoint
## Wide box (100 units): no expansion expected.
wide <- saemix(model, sxdata, control, population = population(100))
wideState <- copulaGet(wide)$lastJoint

stopifnot(smallState$restartCount >= 1L,
  smallState$scoreTheory$projectionExpansions >= 1L,
  wideState$restartCount == 0L)
smallFixed <- as.numeric(small@results@fixed.effects)[1:2]
wideFixed <- as.numeric(wide@results@fixed.effects)[1:2]
covariateMean <- function(fit) copulaGet(fit)$margins[[3L]]$parameters[["mean"]]
cat(sprintf("diagnostic: restarts %d, expansions %d, small V/CL %.2f/%.2f, wide %.2f/%.2f, CRP mean %.1f/%.1f (sample %.1f), small kiter %d
", smallState$restartCount, smallState$scoreTheory$projectionExpansions, smallFixed[1], smallFixed[2], wideFixed[1], wideFixed[2], covariateMean(small), covariateMean(wide), mean(crp), smallState$kiter))
stopifnot(all(abs(smallFixed / wideFixed - 1) < .15),
  abs(covariateMean(small) - covariateMean(wide)) < 5,
  abs(covariateMean(small) - mean(crp)) < 5)
cat(sprintf("restart test passed: %d restart(s), V %.2f/%.2f, CL %.2f/%.2f, CRP mean %.1f/%.1f\n",
  smallState$restartCount, smallFixed[1], wideFixed[1], smallFixed[2],
  wideFixed[2], covariateMean(small), covariateMean(wide)))
