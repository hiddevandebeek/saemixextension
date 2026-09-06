## The hybrid engine uses compiled fixed-family fits as Q-increasing proposals,
## then a simultaneous native-parameter polish. It must retain the direct joint
## Q target while reducing runtime in multi-parameter and higher-d settings.

suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R"); source("R/simEta.R"); source("R/simpleModels.R")

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-60s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}

bench <- function(vTruth, v0, sdTrue, n, seed) {
  set.seed(seed); E <- rEtaVine(n, vTruth, sdTrue); w <- rep(1 / n, n)
  q0 <- sum(w * copulaLogPrior(E, v0, sdTrue * 1.15))
  td <- system.time(direct <- copulaMaximiseJoint(E, w, sdTrue * 1.15, v0,
    length(sdTrue), maxit = 60L, engine = "direct"))["elapsed"]
  th <- system.time(hybrid <- copulaMaximiseJoint(E, w, sdTrue * 1.15, v0,
    length(sdTrue), maxit = 30L, engine = "hybrid", cycles = 2L,
    polishMaxit = 12L))["elapsed"]
  qd <- sum(w * copulaLogPrior(sweep(E, 2, direct$delta, "-"), direct$vine, direct$sd))
  qh <- sum(w * copulaLogPrior(sweep(E, 2, hybrid$delta, "-"), hybrid$vine, hybrid$sd))
  list(q0 = q0, direct = direct, hybrid = hybrid, qd = qd, qh = qh,
       td = as.numeric(td), th = as.numeric(th))
}

vt <- etaVine(list(rvinecopulib::bicop_dist("t", parameters = c(.55, 5))), 2)
v0t <- etaVine(list(rvinecopulib::bicop_dist("t", parameters = c(.2, 20))), 2)
bt <- bench(vt, v0t, c(.30, .42), 5000L, 401L)
ok("H1 hybrid Student step increases the exact common Q", bt$qh > bt$q0,
   sprintf("Q %.5f -> %.5f", bt$q0, bt$qh))
ok("H1b hybrid reaches the direct joint-Q solution", bt$qh >= bt$qd - 2e-4,
   sprintf("direct=%.6f hybrid=%.6f", bt$qd, bt$qh))
ok("H1c hybrid is faster for a multi-parameter pair", bt$th < bt$td,
   sprintf("direct=%.2fs hybrid=%.2fs", bt$td, bt$th))

d <- 4L; vg <- makeTwinVines(etaR(d), "gumbel")$alt
v0g <- etaVineGaussian(etaR(d))
bg <- bench(vg, v0g, rep(.30, d), 2500L, 402L)
ok("H2 hybrid d=4 step increases the exact common Q", bg$qh > bg$q0,
   sprintf("Q %.5f -> %.5f", bg$q0, bg$qh))
ok("H2b hybrid d=4 stays near the direct joint-Q solution", bg$qh >= bg$qd - 1e-3,
   sprintf("direct=%.6f hybrid=%.6f", bg$qd, bg$qh))
ok("H2c hybrid is faster in d=4", bg$th < bg$td,
   sprintf("direct=%.2fs hybrid=%.2fs", bg$td, bg$th))

cat(sprintf("\n%d failure(s)\n", nFail))
if (nFail > 0) quit(status = 1)
