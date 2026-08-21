## Unit test: does the joint (MLE) M-step maximise Q better than the IFM split?
## Deliberately started from a WRONG sd (1.25x) so both have work to do.
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R"); source("R/simEta.R"); source("R/simpleModels.R")
set.seed(3); d <- 4; sdT <- rep(0.30, d)
tw <- makeTwinVines(etaR(d), "gumbel"); tauT <- tw$tauAlt
E <- rEtaVine(4000, tw$alt, sdT); w <- rep(1 / nrow(E), nrow(E))
QT <- sum(w * copulaLogPrior(E, tw$alt, sdT))
cat(sprintf("%-6s %10s %10s %12s\n", "mode", "sd err", "tau err", "Q"))
for (md in c("pool", "sa", "joint")) {
  copulaSet(tw$alt, sdT * 1.25, familySet = "gumbel", mode = md, truncLvl = Inf)
  copulaPoolUpdate(E, 1, 1)
  copulaMstep(50, 0, 1, sdSS = sqrt(colMeans(E^2)), gamma = 1)
  st <- copulaGet()
  tt <- vapply(copulaPadFlat(st$vine, d), function(b)
    if (b$family == "indep") 0 else rvinecopulib::par_to_ktau(b), numeric(1))
  cat(sprintf("%-6s %10.5f %10.5f %12.4f\n", md,
      mean(abs(st$sd - sdT)), mean(abs(tt - tauT)),
      sum(w * copulaLogPrior(E, st$vine, st$sd))))
  copulaClear()
}
cat(sprintf("%-6s %10.5f %10.5f %12.4f   <- truth\n", "TRUTH", 0, 0, QT))
