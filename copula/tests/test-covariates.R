## COVARIATE TESTS FOR THE JOINT COPULA M-STEP.
##
## A location-only joint step is not enough in a population model: every
## estimated covariate coefficient that shifts an eta-bearing parameter enters
## the copula prior and must be part of the same Q maximisation.  These tests
## first pin that algebra directly, then exercise it in a small SAEM nested-null
## fit where an all-Gaussian vine must agree with stock saemix.

suppressMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R")
source("R/simEta.R")
source("R/simpleModels.R")

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-60s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}

## C1 -- direct empirical-Q test with a non-Gaussian copula.  X follows
## saemix's parameter-specific column layout: intercept(V), covariate(V),
## intercept(CL).  locMap maps each coefficient to its eta column.
set.seed(801)
d <- 2L; n <- 5000L
R <- matrix(c(1, .55, .55, 1), 2, 2)
vine <- makeTwinVines(R, "gumbel")$alt
sdTrue <- c(.30, .36)
x <- rep(c(-1, 1), length.out = n)
X <- cbind(1, x, 1)
locMap <- matrix(c(1, 0,
                   1, 0,
                   0, 1), nrow = 3, byrow = TRUE)
betaTrue <- c(.20, .55, -.15)
eta <- rEtaVine(n, vine, sdTrue)
phi <- copulaLocation(X, betaTrue, locMap) + eta
w <- rep(1 / n, n)
beta0 <- c(.28, .05, -.08)
q0 <- sum(w * copulaLogPrior(phi - copulaLocation(X, beta0, locMap),
                             vine, sdTrue * 1.15))
jm <- copulaMaximiseJoint(phi, w, sdTrue * 1.15, vine, d, maxit = 100L,
                          X = X, locMap = locMap, beta0 = beta0,
                          betaFree = seq_along(beta0))
q1 <- sum(w * copulaLogPrior(phi - copulaLocation(X, jm$beta, locMap),
                             jm$vine, jm$sd))
ok("C1 joint Q-step recovers intercepts and covariate effect",
   max(abs(jm$beta - betaTrue)) < .035,
   sprintf("beta=%s maxerr=%.4f", paste(round(jm$beta, 3), collapse = ","),
           max(abs(jm$beta - betaTrue))))
ok("C1b joint covariate step increases empirical Q", q1 > q0 + .05,
   sprintf("Q %.5f -> %.5f", q0, q1))

## C2 -- coefficients excluded by saemix's fixed.estim contract must not move.
jmFixed <- copulaMaximiseJoint(phi, w, sdTrue, vine, d, maxit = 60L,
                               X = X, locMap = locMap, beta0 = beta0,
                               betaFree = c(1L, 3L))
ok("C2 fixed covariate coefficient is preserved exactly",
   identical(jmFixed$beta[2], beta0[2]),
   sprintf("%.6f", jmFixed$beta[2]))

## Small end-to-end Gaussian nested null with a covariate on CL.
simCov <- function(N, vine, betaCov = .40) {
  m <- MODELS$iv1
  z <- rep(c(-1, 1), length.out = N)
  eta <- rEtaVine(N, vine, rep(ETA_SD, 2))
  phi <- cbind(log(m$true[1]) + eta[, 1],
               log(m$true[2]) + betaCov * z + eta[, 2])
  psi <- exp(phi)
  nt <- length(m$times)
  dd <- data.frame(id = rep(seq_len(N), each = nt), dose = m$dose,
                   time = rep(m$times, N), z = rep(z, each = nt))
  f <- m$f(psi, dd$id, cbind(dd$dose, dd$time))
  dd$y <- pmax(f * (1 + .08 * rnorm(nrow(dd))), 1e-6)
  dd
}

covModel <- function() {
  m <- MODELS$iv1
  saemix::saemixModel(
    model = m$f, modeltype = "structural", description = "covariate null",
    psi0 = rbind(m$true, c(0, .20)), transform.par = c(1, 1),
    covariate.model = matrix(c(0, 1), nrow = 1), fixed.estim = c(1, 1),
    covariance.model = matrix(1, 2, 2),
    omega.init = diag(rep(ETA_SD^2, 2)), error.model = "proportional",
    verbose = FALSE)
}
covData <- function(dd) saemix::saemixData(
  name.data = dd, header = TRUE, name.group = "id",
  name.predictors = c("dose", "time"), name.response = "y",
  name.covariates = "z", verbose = FALSE)
ctl <- function(seed) list(
  seed = seed, save = FALSE, save.graphs = FALSE, print = FALSE,
  displayProgress = FALSE, nbiter.saemix = c(100, 60),
  nbiter.mcmc = c(2, 2, 2, 0), warnings = FALSE)

vineG <- etaVineGaussian(R)
fixedDiff <- covDiff <- sdDiff <- numeric(0)
for (rr in 1:2) {
  set.seed(900 + rr)
  dd <- simCov(50, vineG)
  copulaClear()
  fs <- saemix::saemix(covModel(), covData(dd), ctl(rr))
  copulaSet(vineG, rep(ETA_SD, 2), familySet = "gaussian",
            mode = "joint", fitFrom = 20L, refitEvery = 2L,
            populationAlgorithm = "common-q")
  fc <- saemix::saemix(covModel(), covData(dd), ctl(rr))
  st <- copulaGet(); copulaClear()
  bs <- as.numeric(fs@results@fixed.effects)
  bc <- as.numeric(fc@results@fixed.effects)
  fixedDiff <- c(fixedDiff, max(abs(bc[1:2] - bs[1:2])))
  covDiff <- c(covDiff, abs(bc[3] - bs[3]))
  sdDiff <- c(sdDiff, max(abs(st$sd - sqrt(diag(fs@results@omega)))))
}
ok("C3 Gaussian-vine fixed effects match stock with covariates",
   max(fixedDiff) < .08, sprintf("max abs=%.4f", max(fixedDiff)))
ok("C3b Gaussian-vine covariate effect matches stock",
   max(covDiff) < .08, sprintf("max abs=%.4f", max(covDiff)))
ok("C3c Gaussian-vine marginal SDs match stock with covariates",
   max(sdDiff) < .08, sprintf("max abs=%.4f", max(sdDiff)))

## C4 -- conditioning variables own natural-scale, non-centred margins.  The
## eta identification constraints must not leak into these observed variables.
set.seed(812)
naturalCov <- cbind(
  WT = rlnorm(1500, meanlog = log(72), sdlog = .18),
  eGFR = rgamma(1500, shape = 14, scale = 7))
naturalMargins <- copulaFitCovariateMargins(naturalCov)
ind4 <- rvinecopulib::vinecop_dist(list(
  list(rvinecopulib::bicop_dist("indep"),
       rvinecopulib::bicop_dist("indep"),
       rvinecopulib::bicop_dist("indep")),
  list(rvinecopulib::bicop_dist("indep"),
       rvinecopulib::bicop_dist("indep")),
  list(rvinecopulib::bicop_dist("indep"))),
  rvinecopulib::dvine_structure(1:4))
pop4 <- copulaPopulation(ind4,
  margins = list(copulaMarginNormal(.24), copulaMarginNormal(.30)),
  conditioning = naturalCov, conditioningMargins = "auto",
  warmStartOnActivate = FALSE, guard = FALSE)
copulaUsePopulation(pop4)
state4 <- copulaGet(); copulaClear()
ok("C4 automatic natural-scale conditioning margins are accepted",
   state4$dEta == 2L && state4$dConditioning == 2L &&
     all(vapply(state4$margins[1:2], `[[`, logical(1), "centered")) &&
     all(!vapply(state4$margins[3:4], `[[`, logical(1), "centered")),
   paste(vapply(state4$margins, `[[`, character(1), "name"),
               collapse = ","))
ok("C4b support-aware selection recovers lognormal and gamma margins",
   identical(vapply(naturalMargins, `[[`, character(1), "name"),
                    c(WT = "lognormal", eGFR = "gamma")))

cat(sprintf("\n%d failure(s)\n", nFail))
if (nFail > 0) quit(status = 1)
