## End-to-end smoke test of mode="joint" in a real SAEM fit.  Small on purpose:
## d=2, N=60, 2 reps -- enough to show it converges and tracks stock, not enough
## for a precise efficiency claim.
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R"); source("R/simEta.R"); source("R/simpleModels.R")
m <- MODELS$iv1; d <- 2; sdv <- rep(ETA_SD, d); truth <- c(m$true, sdv)
tw <- makeTwinVines(etaR(d), "gumbel")
ctl <- function(s) list(seed = s, save = FALSE, save.graphs = FALSE, print = FALSE,
                        displayProgress = FALSE, nbiter.saemix = c(120, 80),
                        nbiter.mcmc = c(2, 2, 2, 0), warnings = FALSE)
cat(sprintf("%-6s %-6s %-24s %-16s %8s\n", "arm", "truth", "fixed (10, 3)", "sd (.3,.3)", "relErr"))
for (arm in c("null", "alt")) {
  vT <- if (arm == "null") tw$gauss else tw$alt
  for (rr in 1:2) {
    set.seed(500 + rr); s <- simModel(m, 60, vT)
    dat <- saemixDataFor(s$data); mod <- saemixModelFor(m)
    copulaClear()
    fS <- try(saemix::saemix(mod, dat, ctl(rr)), silent = TRUE); if (inherits(fS, "try-error")) next
    sdS <- sqrt(diag(fS@results@omega))
    cat(sprintf("%-6s %-6s %-24s %-16s %8.4f\n", "stock", arm,
        paste(round(fS@results@fixed.effects, 3), collapse = ", "),
        paste(round(sdS, 3), collapse = ", "),
        mean(abs(c(fS@results@fixed.effects, sdS) - truth) / truth)))
    for (md in c("sa", "joint")) {
      copulaSet(etaVineGaussian(cov2cor(fS@results@omega)), sdS,
                familySet = c("gaussian", "gumbel", "clayton", "frank"),
                mode = md, refitEvery = 3L, fitFrom = 25L)
      fC <- try(saemix::saemix(mod, dat, ctl(rr)), silent = TRUE)
      if (inherits(fC, "try-error")) { cat("  ", md, "FAILED\n"); copulaClear(); next }
      st <- copulaGet(); copulaClear()
      cat(sprintf("%-6s %-6s %-24s %-16s %8.4f\n", md, arm,
          paste(round(fC@results@fixed.effects, 3), collapse = ", "),
          paste(round(st$sd, 3), collapse = ", "),
          mean(abs(c(fC@results@fixed.effects, st$sd) - truth) / truth)))
    }
  }
}
