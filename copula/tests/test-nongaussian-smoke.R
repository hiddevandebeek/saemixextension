## Small acceptance smoke: a correctly specified fixed Gumbel candidate should
## improve the joint upper tail over the Gaussian random-effects model without
## materially degrading ordinary population/scale estimates.

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

set.seed(1201)
m <- MODELS$iv1; d <- 2L; N <- 80L; sdv <- rep(ETA_SD, d)
R <- matrix(c(1, .70, .70, 1), 2, 2)
truthVine <- makeTwinVines(R, "gumbel")$alt
sim <- simModel(m, N, truthVine, propErr = .08)
dat <- saemixDataFor(sim$data); mod <- saemixModelFor(m)
ctl <- list(seed = 41, save = FALSE, save.graphs = FALSE, print = FALSE,
            displayProgress = FALSE, nbiter.saemix = c(100, 60),
            nbiter.mcmc = c(2, 2, 2, 0), warnings = FALSE)

copulaClear()
fitG <- saemix::saemix(mod, dat, ctl)
sdG <- sqrt(diag(fitG@results@omega))
RG <- cov2cor(fitG@results@omega)
tau0 <- 2 / pi * asin(RG[1, 2])
g0 <- rvinecopulib::bicop_dist("gumbel", parameters = 1 / (1 - tau0))
startVine <- etaVine(list(g0), d)
copulaSet(startVine, sdG, familySet = NULL, mode = "joint",
          fitFrom = 20L, refitEvery = 2L, guardAction = "warn")
fitC <- saemix::saemix(mod, dat, ctl)
st <- copulaGet(); copulaClear()

pairHigh <- function(vine, q = .10, n = 100000L) {
  u <- withSeed(77, rvinecopulib::rvinecop(n, vine))
  mean(u[, 1] > 1 - q & u[, 2] > 1 - q) / q
}
tTrue <- pairHigh(truthVine)
tG <- pairHigh(etaVineGaussian(RG))
tC <- pairHigh(st$vine)
eG <- abs(tG - tTrue); eC <- abs(tC - tTrue)
ok("G1 fixed Gumbel joint fit improves upper-tail recovery", eC < eG,
   sprintf("true=%.3f gaussian=%.3f copula=%.3f", tTrue, tG, tC))

truth <- c(m$true, sdv)
parG <- c(fitG@results@fixed.effects, sdG)
parC <- c(fitC@results@fixed.effects, st$sd)
errG <- mean(abs(parG - truth) / truth)
errC <- mean(abs(parC - truth) / truth)
ok("G2 non-Gaussian tail gain keeps parameter error comparable",
   errC <= errG + .03, sprintf("gaussian=%.4f copula=%.4f", errG, errC))
ok("G3 final non-Gaussian model remained discretely frozen",
   isTRUE(st$modelFrozen) && !isTRUE(st$collapseDetected))
ok("G4 final SAEM iteration performs a full direct-Q polish",
   isTRUE(st$lastJoint$final) && identical(st$lastJoint$conv, 0L),
   sprintf("iteration=%s convergence=%s", st$lastJoint$kiter, st$lastJoint$conv))

cat(sprintf("\n%d failure(s)\n", nFail))
if (nFail > 0) quit(status = 1)
