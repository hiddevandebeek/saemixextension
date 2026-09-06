## Every native parameter of a fixed parametric pair family belongs in Q.
## This pins the parameters that a Kendall-tau-only optimiser necessarily held
## fixed: Student df and both BB shape coordinates.

suppressMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R")

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-60s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}

fitPrior <- function(truth, start, n = 7000L, seed = 1L) {
  set.seed(seed)
  sdTrue <- c(.30, .42)
  E <- rEtaVine(n, etaVine(list(truth), 2), sdTrue)
  w <- rep(1 / n, n)
  v0 <- etaVine(list(start), 2)
  q0 <- sum(w * copulaLogPrior(E, v0, sdTrue * 1.15))
  jm <- copulaMaximiseJoint(E, w, sdTrue * 1.15, v0, 2, maxit = 200L)
  q1 <- sum(w * copulaLogPrior(E, jm$vine, jm$sd))
  list(jm = jm, q0 = q0, q1 = q1,
       par = as.numeric(etaVineFlat(jm$vine)[[1]]$parameters))
}

rt <- fitPrior(
  rvinecopulib::bicop_dist("t", parameters = c(.55, 5)),
  rvinecopulib::bicop_dist("t", parameters = c(.20, 20)), seed = 211)
ok("N1 Student rho and df are both optimized",
   length(rt$par) == 2L && abs(rt$par[2] - 20) > 2 && rt$q1 > rt$q0 + .01,
   sprintf("par=(%.3f,%.2f), dQ=%.4f", rt$par[1], rt$par[2], rt$q1 - rt$q0))
ok("N1b Student fit recovers both native parameters",
   abs(rt$par[1] - .55) < .08 && abs(rt$par[2] - 5) < 3,
   sprintf("par=(%.3f,%.2f)", rt$par[1], rt$par[2]))

rb <- fitPrior(
  rvinecopulib::bicop_dist("bb1", parameters = c(.8, 1.8)),
  rvinecopulib::bicop_dist("bb1", parameters = c(.2, 1.15)), seed = 212)
ok("N2 both BB1 parameters are optimized",
   length(rb$par) == 2L && all(abs(rb$par - c(.2, 1.15)) > .05) &&
     rb$q1 > rb$q0 + .01,
   sprintf("par=(%.3f,%.3f), dQ=%.4f", rb$par[1], rb$par[2], rb$q1 - rb$q0))
ok("N2b BB1 fit recovers both native parameters",
   max(abs(rb$par - c(.8, 1.8))) < .25,
   sprintf("par=(%.3f,%.3f)", rb$par[1], rb$par[2]))

## Parameter counting is the continuous dimension used by BIC later.
vf <- etaVine(list(
  rvinecopulib::bicop_dist("t", parameters = c(.4, 6)),
  rvinecopulib::bicop_dist("bb1", parameters = c(.5, 1.5)),
  rvinecopulib::bicop_dist("indep")), 3)
lay <- copulaParameterLayout(copulaPadFlat(vf, 3))
ok("N3 native parameter count includes t df and BB shape", length(lay$par) == 4L,
   sprintf("npar=%d", length(lay$par)))

cat(sprintf("\n%d failure(s)\n", nFail))
if (nFail > 0) quit(status = 1)
