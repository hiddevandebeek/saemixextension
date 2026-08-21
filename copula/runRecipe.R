## THE RECIPE TEST: does the robust default (TRUNCATED vine) keep the tail win?
##
## Truncation fixed the estimation failure under a Gaussian truth, but the
## Gumbel twin has non-Gaussian pair copulas on EVERY edge, including the higher
## trees that truncation discards.  So truncation might buy robustness by
## throwing away exactly the structure the copula is there to capture.
##
## Arms, all on the same data / seed, d=4, GUMBEL truth:
##   stock  multivariate-normal random effects (the incumbent)
##   full   all 3 trees, family search
##   t1     truncated after tree 1, family search
## Metric: joint tail probabilities against the truth (computed by direct
## simulation of the fitted vines, so no importance sampling is involved and the
## comparison is not affected by the low-ESS problem), plus parameter error.
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R"); source("R/simEta.R"); source("R/simpleModels.R")
NREP <- 6; N <- 150
m <- MODELS$iv2; d <- 4; sdv <- rep(ETA_SD, d); truth <- c(m$true, sdv)
tw <- makeTwinVines(etaR(d), "gumbel")
FAMSET <- c("gaussian", "clayton", "gumbel", "frank")
ctl <- function(sd) list(seed = sd, save = FALSE, save.graphs = FALSE, print = FALSE,
                         displayProgress = FALSE, nbiter.saemix = c(150, 100),
                         nbiter.mcmc = c(2, 2, 2, 0), warnings = FALSE)
corners <- function(vine, q = 0.10, nsim = 2e5) {
  u <- withSeed(31, rvinecopulib::rvinecop(nsim, vine))
  c(allHigh = mean(rowSums(u > 1 - q) == ncol(u)),
    allLow  = mean(rowSums(u < q) == ncol(u)),
    pairHigh= mean(u[, 1] > 1 - q & u[, 2] > 1 - q) / q)
}
cTrue <- corners(tw$alt)
cat(sprintf("TRUTH: allHigh=%.4f allLow=%.4f pairHigh=%.4f\n\n", cTrue[1], cTrue[2], cTrue[3]))
res <- list()
for (rr in seq_len(NREP)) {
  set.seed(7000 + rr)
  s <- simModel(m, N, tw$alt); dat <- saemixDataFor(s$data); mod <- saemixModelFor(m)
  copulaClear()
  fS <- try(saemix::saemix(mod, dat, ctl(rr)), silent = TRUE); if (inherits(fS, "try-error")) next
  sdS <- sqrt(diag(fS@results@omega)); vnS <- etaVineGaussian(cov2cor(fS@results@omega))
  row <- list(rep = rr, arm = "stock", est = c(fS@results@fixed.effects, sdS), cor = corners(vnS))
  res[[length(res) + 1]] <- row
  for (nm in c("full", "t1")) {
    copulaSet(etaVineGaussian(cov2cor(fS@results@omega)), sdS, familySet = FAMSET,
              mode = "sa", refitEvery = 3L, fitFrom = 30L,
              truncLvl = if (nm == "t1") 1L else Inf)
    fC <- try(saemix::saemix(mod, dat, ctl(rr)), silent = TRUE)
    if (inherits(fC, "try-error")) { copulaClear(); next }
    st <- copulaGet(); copulaClear()
    res[[length(res) + 1]] <- list(rep = rr, arm = nm,
      est = c(fC@results@fixed.effects, st$sd), cor = corners(st$vine))
  }
  cat(sprintf("rep %d done\n", rr)); flush.console()
}
saveRDS(res, "out/recipe.rds")
cat(sprintf("\n%-7s %9s %28s %28s %28s\n", "arm", "parErr", "allHigh (true .0114)",
            "allLow (true .0017)", "pairHigh (true .49)"))
for (a in c("stock", "full", "t1")) {
  k <- vapply(res, function(x) x$arm == a, logical(1)); if (!any(k)) next
  E <- do.call(rbind, lapply(res[k], function(x) abs(x$est - truth) / truth))
  C <- do.call(rbind, lapply(res[k], function(x) x$cor))
  cat(sprintf("%-7s %9.4f %13.4f (err %+.4f) %13.4f (err %+.4f) %13.4f (err %+.4f)\n",
      a, mean(E), mean(C[, 1]), mean(C[, 1]) - cTrue[1], mean(C[, 2]), mean(C[, 2]) - cTrue[2],
      mean(C[, 3]), mean(C[, 3]) - cTrue[3]))
}
cat("\nper-replicate |pairHigh error|:\n")
for (a in c("stock", "full", "t1")) {
  k <- vapply(res, function(x) x$arm == a, logical(1)); if (!any(k)) next
  e <- vapply(res[k], function(x) abs(x$cor[3] - cTrue[3]), numeric(1))
  cat(sprintf("  %-6s %s  mean=%.4f\n", a, paste(sprintf("%.3f", e), collapse = " "), mean(e)))
}
