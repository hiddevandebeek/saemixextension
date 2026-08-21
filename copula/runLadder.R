## WHERE DOES IT BREAK?
##
## For each model on the ladder, two truths:
##   null  Gaussian vine  -- the copula path should MATCH stock (same model)
##   alt   Gumbel vine    -- the copula path should RECOVER the dependence
## and for the null, stock vs copula on the same data and seed.
##
## d=2 is the minimal case: exactly ONE pair copula, so recovery is a
## one-parameter problem and nothing can hide.
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R"); source("R/simEta.R"); source("R/simpleModels.R")

a <- commandArgs(trailingOnly = TRUE)
NREP <- if (length(a) > 0) as.integer(a[1]) else 5
N    <- if (length(a) > 1) as.integer(a[2]) else 100
ctl <- function(sd) list(seed = sd, save = FALSE, save.graphs = FALSE, print = FALSE,
                         displayProgress = FALSE, nbiter.saemix = c(150, 100),
                         nbiter.mcmc = c(2, 2, 2, 0), warnings = FALSE)
tauOf <- function(v) vapply(unlist(v$pair_copulas, recursive = FALSE),
  function(b) if (b$family == "indep") 0 else rvinecopulib::par_to_ktau(b), numeric(1))
famOf <- function(v) paste(vapply(unlist(v$pair_copulas, recursive = FALSE),
  function(b) b$family, character(1)), collapse = ",")

res <- list()
for (nm in names(MODELS)) {
  m <- MODELS[[nm]]; d <- m$d
  R <- etaR(d); sdv <- rep(ETA_SD, d)
  tw <- makeTwinVines(R, "gumbel")
  truth <- c(m$true, sdv)
  for (arm in c("null", "alt")) {
    vTrue <- if (arm == "null") tw$gauss else tw$alt
    tauTrue <- tauOf(vTrue)
    for (rr in seq_len(NREP)) {
      set.seed(1000 * which(names(MODELS) == nm) + 100 * (arm == "alt") + rr)
      s <- simModel(m, N, vTrue)
      dat <- saemixDataFor(s$data); mod <- saemixModelFor(m)
      copulaClear()
      fS <- try(saemix::saemix(mod, dat, ctl(rr)), silent = TRUE)
      if (inherits(fS, "try-error")) next
      sdS <- sqrt(diag(fS@results@omega))
      ## copula arm: start from the stock correlation, Gaussian families only in
      ## the null (nothing to find) and a family search in the alt arm
      fams <- if (arm == "null") "gaussian" else c("gaussian", "clayton", "gumbel", "frank")
      copulaSet(etaVineGaussian(cov2cor(fS@results@omega)), sdS,
                familySet = fams, mode = "sa", refitEvery = 3L, fitFrom = 30L)
      fC <- try(saemix::saemix(mod, dat, ctl(rr)), silent = TRUE)
      if (inherits(fC, "try-error")) { copulaClear(); next }
      st <- copulaGet(); copulaClear()
      res[[length(res) + 1]] <- data.frame(
        model = nm, d = d, arm = arm, rep = rr,
        errStock = mean(abs(c(fS@results@fixed.effects, sdS) - truth) / truth),
        errCop   = mean(abs(c(fC@results@fixed.effects, st$sd) - truth) / truth),
        tauTrue  = paste(round(tauTrue, 3), collapse = ","),
        tauStock = paste(round(tauOf(etaVineGaussian(cov2cor(fS@results@omega))), 3), collapse = ","),
        tauCop   = paste(round(tauOf(st$vine), 3), collapse = ","),
        tauErrStock = mean(abs(tauOf(etaVineGaussian(cov2cor(fS@results@omega))) - tauTrue)),
        tauErrCop   = mean(abs(tauOf(st$vine) - tauTrue)),
        famCop = famOf(st$vine), stringsAsFactors = FALSE)
      cat(sprintf("%-6s d=%d %-4s rep%d  err s/c = %.4f/%.4f   tauErr s/c = %.3f/%.3f  %s\n",
          nm, d, arm, rr, res[[length(res)]]$errStock, res[[length(res)]]$errCop,
          res[[length(res)]]$tauErrStock, res[[length(res)]]$tauErrCop,
          res[[length(res)]]$famCop)); flush.console()
    }
  }
}
res <- do.call(rbind, res); saveRDS(res, "out/ladder.rds")
cat(sprintf("\n%-6s %-3s %-5s %9s %9s %7s %10s %10s\n",
            "model", "d", "arm", "errStock", "errCop", "ratio", "tauErrStk", "tauErrCop"))
for (nm in names(MODELS)) for (arm in c("null", "alt")) {
  k <- res$model == nm & res$arm == arm; if (!any(k)) next
  cat(sprintf("%-6s %-3d %-5s %9.4f %9.4f %7.2f %10.3f %10.3f\n", nm, res$d[k][1], arm,
      mean(res$errStock[k]), mean(res$errCop[k]),
      mean(res$errCop[k]) / mean(res$errStock[k]),
      mean(res$tauErrStock[k]), mean(res$tauErrCop[k])))
}
