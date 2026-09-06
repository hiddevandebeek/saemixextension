## The score recursion against plain saemix, on identical data.
##
## With Normal margins and a Gaussian copula the population is an ordinary
## multivariate normal on the transformed scale, which is exactly what plain
## saemix fits. So the same data can be fitted by two routes that share no M
## step, no parameterisation and no uncertainty machinery: saemix estimates
## omega directly in closed form, the copula route reaches the same
## distribution by a stochastic gradient on the marginal likelihood carrying
## marginal spreads and a correlation.
##
## They have to agree. This is the check that catches bad parameters, and it
## has already earned its place: it is what showed that the Fisher metric on
## this scale lands 28% low on the correlation, from any starting value and at
## any length, while the per-coordinate metric agrees to three decimals.
##
## The tolerances are roughly thirty times the differences actually observed,
## which leaves room for platform arithmetic while still being an order of
## magnitude tighter than the failure it was written to catch.

suppressPackageStartupMessages({
  library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)
  library(rvinecopulib)
})
source("C:/package/saemix-copula/copula/experiments/combined-natural-frem-study/functions.R")

nFail <- 0L
ok <- function(label, value, tolerance) {
  pass <- isTRUE(abs(value) <= tolerance)
  if (!pass) nFail <<- nFail + 1L
  cat(sprintf("%-46s %s  %+.5f (tolerance %.4f)\n", label,
    if (pass) "PASS" else "**FAIL**", value, tolerance))
}

set.seed(8301000L)
truth <- combined_truth(); n <- 250L
rho <- .45; sdV <- truth$sdlogV; sdCL <- truth$sdlogV * 1.4
correlation <- matrix(c(1, rho, rho, 1), 2L, 2L)
eta <- sweep(matrix(rnorm(n * 2L), n, 2L) %*% chol(correlation), 2L,
  c(sdV, sdCL), "*")
psi <- sweep(exp(eta), 2L, c(truth$V, truth$CL), "*")
dd <- data.frame(id = rep(seq_len(n), each = length(truth$times)),
  dose = truth$dose, time = rep(truth$times, n))
dd$y <- pmax(combined_pk(psi, dd$id, cbind(dd$dose, dd$time)) *
  (1 + truth$residual * rnorm(nrow(dd))), 1e-8)
sxdata <- combined_data(dd)

plain <- saemix(combined_model(), sxdata,
  list(seed = 8301011L, save = FALSE, save.graphs = FALSE, print = FALSE,
    displayProgress = FALSE, nbiter.saemix = c(400L, 400L), warnings = FALSE,
    ll.is = FALSE, fim = FALSE, map = FALSE))
omega <- plain@results@omega[1:2, 1:2, drop = FALSE]

## The copula route, started deliberately far from the answer so that the
## comparison is about where it arrives and not about where it began.
vine <- copulaGaussianRvineFromCor(matrix(c(1, .1, .1, 1), 2L, 2L),
  rvinecopulib::dvine_structure(1:2))
fit <- saemix(combined_model(), sxdata, combined_control(8301060L, 800L),
  population = copulaPopulation(vine,
    margins = list(copulaMarginNormal(.25), copulaMarginNormal(.30)),
    scale = "transformed-additive", populationAlgorithm = "score-sa",
    scoreScale = "auto", scoreBurn = 50L))
state <- copulaGet(fit)

ok("typical value of the first parameter",
  as.numeric(fit@results@fixed.effects)[1L] -
    as.numeric(plain@results@fixed.effects)[1L], .05)
ok("typical value of the second parameter",
  as.numeric(fit@results@fixed.effects)[2L] -
    as.numeric(plain@results@fixed.effects)[2L], .01)
ok("spread of the first random effect",
  unname(state$margins[[1L]]$parameters[1L]) - sqrt(omega[1L, 1L]), .005)
ok("spread of the second random effect",
  unname(state$margins[[2L]]$parameters[1L]) - sqrt(omega[2L, 2L]), .005)
ok("correlation between the random effects",
  copulaGaussianRvineCor(state$vine, 2L)[1L, 2L] -
    omega[1L, 2L] / sqrt(omega[1L, 1L] * omega[2L, 2L]), .01)
ok("residual error",
  as.numeric(fit@results@respar)[2L] - as.numeric(plain@results@respar)[2L],
  .001)

cat(if (nFail) sprintf("\n%d FAILURES\n", nFail) else
  "\nthe two routes agree\n")
if (nFail) quit(status = 1L)
