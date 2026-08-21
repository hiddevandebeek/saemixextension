## NESTED-NULL TESTS.
##
## An all-Gaussian vine IS MVN(0,Omega) (test-ordering.R T1), so the copula path
## must reproduce stock saemix.  If it does not, the copula path is wrong -- this
## is the strongest available check, because it compares the SAME estep/mstep
## code under two priors that are provably the same distribution.

suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)
                  library(mvtnorm)})
setwd("C:/package/saemix-copula/copula")
source("R/simData.R")

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-56s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}

d <- 4; sdv <- PK_SD
Om <- diag(sdv) %*% PK_R %*% diag(sdv)
vnG <- etaVineGaussian(PK_R)
copulaSet(vnG, sdv)

## N1 -- copulaUeta differs from the saemix quadratic form by a CONSTANT only.
## The kernels use differences, so a constant offset is harmless; a non-constant
## difference would mean a different prior.
set.seed(5); E <- rmvnorm(500, sigma = Om)
qf <- 0.5 * rowSums(E * (E %*% solve(Om)))
dv <- copulaUeta(E) - qf
ok("N1 copulaUeta - Gaussian quadratic form is constant",
   diff(range(dv)) < 1e-9, sprintf("spread=%.2e const=%.4f", diff(range(dv)), mean(dv)))

## N2 -- and that constant is exactly the MVN log normalising constant.
ok("N2 constant == -log N(0,Omega) normalisation",
   abs(mean(dv) - 0.5 * (d * log(2 * pi) + determinant(Om, TRUE)$modulus)) < 1e-9,
   sprintf("%.6f", mean(dv)))

## N3 -- copulaOmega() recovers Omega.
e3 <- max(abs(copulaOmega() - Om))
ok("N3 copulaOmega() == Omega for a Gaussian vine", e3 < 5e-3, sprintf("maxabs=%.2e", e3))

## N4 -- kernel-1 proposals are drawn from the right prior.
set.seed(6); P <- copulaRandEta(60000)
ok("N4 copulaRandEta covariance == Omega", max(abs(cov(P) - Om)) < 5e-3,
   sprintf("maxabs=%.2e", max(abs(cov(P) - Om))))

## N5 -- END TO END across 3 datasets: same data, same seed, copula path vs
## stock saemix.  Not bit-identical (kernel 1 consumes the RNG differently), so
## compare statistically.  Split by identifiability: the fixed effects and the
## well-determined etas must agree tightly; V2 (peripheral volume) is weakly
## identified in this design -- stock itself lands at sd 0.335 against a truth of
## 0.40 -- so there the requirement is that the copula path is no FURTHER from
## truth than stock, not that the two agree with each other.
ctl <- function(sd) list(seed = sd, save = FALSE, save.graphs = FALSE, print = FALSE,
                         displayProgress = FALSE, nbiter.saemix = c(150, 80),
                         nbiter.mcmc = c(2, 2, 2, 0), warnings = FALSE)
fixRel <- sdRel <- c(); errS <- errC <- c()
for (rr in 1:3) {
  set.seed(70 + rr); sD <- simPK(120, vnG)
  copulaClear()
  fS <- saemix::saemix(pkSaemixModel(), pkSaemixData(sD$data), ctl(rr))
  copulaSet(etaVineGaussian(PK_R), sdv, familySet = "gaussian")
  fC <- saemix::saemix(pkSaemixModel(), pkSaemixData(sD$data), ctl(rr))
  sdC <- copulaGet()$sd; copulaClear()
  fixS <- fS@results@fixed.effects; fixC <- fC@results@fixed.effects
  sdS <- sqrt(diag(fS@results@omega))
  fixRel <- c(fixRel, max(abs(fixC - fixS) / abs(fixS)))
  sdRel  <- c(sdRel,  max(abs(sdC[1:3] - sdS[1:3]) / sdS[1:3]))
  errS <- c(errS, abs(sdS[4] - PK_SD[4])); errC <- c(errC, abs(sdC[4] - PK_SD[4]))
}
ok("N5 fixed effects match stock (3 datasets)", max(fixRel) < 0.06,
   sprintf("max rel=%.3f", max(fixRel)))
ok("N5b sd of well-identified etas match stock", max(sdRel) < 0.10,
   sprintf("max rel=%.3f", max(sdRel)))
ok("N5c weakly-identified sd(V2) no worse than stock", mean(errC) <= mean(errS) + 0.02,
   sprintf("|err| stock=%.3f copula=%.3f (truth %.2f)", mean(errS), mean(errC), PK_SD[4]))
copulaSet(etaVineGaussian(PK_R), sdv, familySet = "gaussian")

## N6 -- the copula path recovers the correlation structure it was given.
e6 <- max(abs(cov2cor(copulaOmega()) - PK_R))
ok("N6 copula path recovers eta correlations", e6 < 0.25, sprintf("maxabs=%.3f", e6))

## N7 -- copulaOmega() must NOT disturb the caller's RNG stream.  It uses fixed
## common random numbers internally; a bare set.seed() there would reset R's
## global stream inside the SAEM loop and make the E-step replay identical draws
## every iteration (this was a real bug -- the fit ran away monotonically).
set.seed(123); a <- rnorm(5)
set.seed(123); invisible(copulaOmega()); b <- rnorm(5)
ok("N7 copulaOmega() preserves the RNG stream", isTRUE(all.equal(a, b)))
copulaClear()

cat(sprintf("\n%d failure(s)\n", nFail))
if (nFail > 0) quit(status = 1)
