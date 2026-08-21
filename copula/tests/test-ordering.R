## Ordering / Rosenblatt sanity tests -- exact (density-based), not Monte Carlo.
##
## Everything downstream (kernel-1 prior sampling, the U.eta prior density, the
## non-centered eta = T(z) map) depends on the eta column order agreeing with the
## vine's variable order.  A mismatch is SILENT: still a valid joint law, just
## not the one specified.

suppressMessages({library(rvinecopulib); library(VineCopula); library(mvtnorm)})
source("R/etaCopula.R"); source("R/simEta.R")

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-56s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}

set.seed(20260821)
d <- 4
R <- matrix(c(1, .55, .30, .10,
              .55, 1, .45, .20,
              .30, .45, 1, .40,
              .10, .20, .40, 1), d, d)
stopifnot(all(eigen(R)$values > 0))
u <- matrix(runif(400 * d), ncol = d); z <- qnorm(u)
dGaussCop <- mvtnorm::dmvnorm(z, sigma = R) / apply(dnorm(z), 1, prod)

## T1 -- the Gaussian vine IS MVN(0,R).  Validates partial-correlation
## extraction, canonical edge order, and column identity in one exact check.
vnG <- etaVineGaussian(R)
e1 <- max(abs(rvinecopulib::dvinecop(u, vnG) - dGaussCop) / dGaussCop)
ok("T1 Gaussian vine density == MVN copula density", e1 < 1e-9, sprintf("relerr=%.2e", e1))

## T2 -- eta-scale density matches dmvnorm with Omega = diag(sd) R diag(sd).
sdv <- c(.30, .35, .45, .25)
Om  <- diag(sdv) %*% R %*% diag(sdv)
eta <- matrix(rnorm(400 * d), ncol = d) %*% chol(Om)
e2 <- max(abs(dEtaVine(eta, vnG, sdv) - mvtnorm::dmvnorm(eta, sigma = Om, log = TRUE)))
ok("T2 dEtaVine == dmvnorm(Omega) for the Gaussian vine", e2 < 1e-8, sprintf("maxabs=%.2e", e2))

## T3 -- independent oracle: VineCopula, whose array convention is TRANSPOSED
## relative to rvinecopulib.  Built separately from the same partial
## correlations; densities must agree.  Catches a convention drift in either.
vcG <- VineCopula::D2RVine(1:d, family = rep(1, length(dvinePartialCor(R))),
                           par = dvinePartialCor(R))
e3 <- max(abs(VineCopula::RVinePDF(u, vcG) - dGaussCop) / dGaussCop)
ok("T3 VineCopula oracle agrees (transposed convention)", e3 < 1e-9, sprintf("relerr=%.2e", e3))

## T4 -- column identity under a deliberately asymmetric vine: dependence only
## on edge (1,2).  A reversed/transposed array would move it to (3,4).
flatA <- c(list(rvinecopulib::bicop_dist("gaussian", parameters = 0.9)),
           rep(list(rvinecopulib::bicop_dist("indep")), 5))
RA <- cor(rEtaVine(20000, etaVine(flatA, d), rep(1, d)))
offIdx <- cbind(c(1, 1, 2, 2, 3), c(3, 4, 3, 4, 4))
ok("T4 asymmetric vine keeps dependence on cols (1,2)",
   abs(RA[1, 2] - 0.9) < 0.02 && max(abs(RA[offIdx])) < 0.03,
   sprintf("r12=%.3f maxOther=%.3f", RA[1, 2], max(abs(RA[offIdx]))))

## T5 -- Rosenblatt round trip, on the ETA scale (the form the E-step needs).
tw <- makeTwinVines(R, family = "gumbel")
etaT <- rEtaVine(2000, tw$alt, sdv)
zz   <- zFromEta(etaT, tw$alt, sdv)
back <- etaFromZ(zz, tw$alt, sdv)
e5 <- max(abs(back - etaT))
ok("T5 etaFromZ o zFromEta == identity (non-Gaussian vine)", e5 < 1e-6, sprintf("maxabs=%.2e", e5))
ok("T5b Rosenblatt scores z are decorrelated N(0,1)",
   max(abs(cor(zz)[upper.tri(diag(d))])) < 0.06 && max(abs(apply(zz, 2, sd) - 1)) < 0.05,
   sprintf("maxcor=%.3f", max(abs(cor(zz)[upper.tri(diag(d))]))))

## T6 -- ORDER MATTERS FOR THE MAP, NOT FOR THE LAW.
## etaVineReversed() relabels a D-vine to the reversed variable order (edges
## reversed within each tree).  For EXCHANGEABLE pair copulas that is a genuine
## relabelling: identical density everywhere.  But the inverse-Rosenblatt map
## still sends the same z to a different eta -- true even for a plain MVN.
etaChk <- rEtaVine(2000, tw$gauss, sdv)
zfix   <- matrix(rnorm(2000 * d), ncol = d)
vnRevG <- etaVineReversed(etaVineFlat(tw$gauss), d)
dSame  <- max(abs(dEtaVine(etaChk, tw$gauss, sdv) - dEtaVine(etaChk, vnRevG, sdv)))
ok("T6 relabelled Gaussian vine: IDENTICAL density", dSame < 1e-10, sprintf("maxabs=%.2e", dSame))
mapDiff <- max(abs(etaFromZ(zfix, tw$gauss, sdv) - etaFromZ(zfix, vnRevG, sdv)))
ok("T6b relabelled vine: DIFFERENT map, even for MVN", mapDiff > 1e-3,
   sprintf("maxdiff=%.3f -- eta order must be pinned", mapDiff))

## T6c -- and reversal is a relabelling ONLY for exchangeable families.  The
## twin vine carries one rotation-90 edge (the negative-tau edge), and rotated
## copulas are not exchangeable, so the reversed vine is a DIFFERENT MODEL.
## Consequence for the design: with rotations in play the eta order is part of
## the model specification, not a free choice.  Assert the difference so a
## future refactor that quietly reorders etas fails here.
flatGum <- rep(list(rvinecopulib::bicop_dist("gumbel", parameters = 1.8)), d * (d - 1) / 2)
dGum <- max(abs(dEtaVine(etaChk, etaVine(flatGum, d), sdv) -
                dEtaVine(etaChk, etaVineReversed(flatGum, d), sdv)))
ok("T6c unrotated Gumbel vine reverses cleanly", dGum < 1e-10, sprintf("maxabs=%.2e", dGum))
rots <- vapply(etaVineFlat(tw$alt), function(b) b$rotation, numeric(1))
dRot <- max(abs(dEtaVine(etaChk, tw$alt, sdv) -
                dEtaVine(etaChk, etaVineReversed(etaVineFlat(tw$alt), d), sdv)))
ok("T6d ROTATED vine does NOT reverse (order is part of model)", dRot > 1e-3,
   sprintf("maxabs=%.2e, rotations=%s", dRot, paste(rots, collapse = ",")))

## T7 -- twin vines share edgewise Kendall tau exactly.
ok("T7 twin vines share edgewise Kendall tau", max(abs(tw$tau - tw$tauAlt)) < 1e-6,
   sprintf("maxdTau=%.2e", max(abs(tw$tau - tw$tauAlt))))

## T8 -- twins differ in SHAPE not STRENGTH: matched Kendall tau on the eta
## scale, but different tail dependence.
eg <- rEtaVine(40000, tw$gauss, sdv); ea <- rEtaVine(40000, tw$alt, sdv)
sp <- function(x) cor(x, method = "spearman")
dSp <- max(abs(sp(eg) - sp(ea))[upper.tri(diag(d))])
q <- 0.95
upG <- mean(eg[, 1] > quantile(eg[, 1], q) & eg[, 2] > quantile(eg[, 2], q)) / (1 - q)
upA <- mean(ea[, 1] > quantile(ea[, 1], q) & ea[, 2] > quantile(ea[, 2], q)) / (1 - q)
ok("T8 twins match on Spearman rho", dSp < 0.05, sprintf("maxdRho=%.3f", dSp))
ok("T8b twins differ in upper tail dependence", (upA - upG) > 0.10,
   sprintf("lambdaU_hat: gauss=%.3f gumbel=%.3f", upG, upA))

cat(sprintf("\n%d failure(s)\n", nFail))
if (nFail > 0) quit(status = 1)
