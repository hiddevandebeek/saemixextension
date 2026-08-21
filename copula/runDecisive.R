## DIMENSION or IDENTIFIABILITY?
##
## Same d=4 structural model, same eta SDs, same correlation, same N -- only the
## DESIGN differs, so shrinkage changes and dimension does not:
##   iv2      10 obs/subj, 10% error -> V2 shrinkage 0.21
##   iv2rich  15 obs/subj,  5% error -> V2 shrinkage 0.03
## Truth is the Gaussian vine in both, so the copula path is fitting exactly the
## right model and any gap is pure estimator cost.
##
## If the gap is present in iv2 and absent in iv2rich, the cause is
## identifiability and the prescription is "put the vine only on etas the data
## informs".  If it is present in both, the cause is dimension and the
## prescription is a truncated vine.
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R"); source("R/simEta.R"); source("R/simpleModels.R")
NREP <- 5; N <- 100
ctl <- function(sd) list(seed = sd, save = FALSE, save.graphs = FALSE, print = FALSE,
                         displayProgress = FALSE, nbiter.saemix = c(150, 100),
                         nbiter.mcmc = c(2, 2, 2, 0), warnings = FALSE)
res <- list()
for (nm in c("iv2", "iv2rich")) {
  m <- MODELS[[nm]]; d <- m$d; sdv <- rep(ETA_SD, d)
  pe <- if (is.null(m$propErr)) 0.10 else m$propErr
  vn <- etaVineGaussian(etaR(d)); truth <- c(m$true, sdv)
  for (rr in seq_len(NREP)) {
    set.seed(3000 + rr)
    s <- simModel(m, N, vn, propErr = pe)
    dat <- saemixDataFor(s$data); mod <- saemixModelFor(m)
    copulaClear()
    fS <- try(saemix::saemix(mod, dat, ctl(rr)), silent = TRUE); if (inherits(fS, "try-error")) next
    sdS <- sqrt(diag(fS@results@omega))
    shr <- 1 - apply(fS@results@cond.mean.psi, 2, function(x) sd(log(x))) / sdS
    copulaSet(etaVineGaussian(cov2cor(fS@results@omega)), sdS, familySet = "gaussian",
              mode = "sa", refitEvery = 3L, fitFrom = 30L)
    fC <- try(saemix::saemix(mod, dat, ctl(rr)), silent = TRUE)
    if (inherits(fC, "try-error")) { copulaClear(); next }
    st <- copulaGet(); copulaClear()
    eS <- abs(c(fS@results@fixed.effects, sdS) - truth) / truth
    eC <- abs(c(fC@results@fixed.effects, st$sd) - truth) / truth
    res[[length(res) + 1]] <- data.frame(model = nm, rep = rr,
      errStock = mean(eS), errCop = mean(eC), maxShrink = max(shr),
      perStock = I(list(eS)), perCop = I(list(eC)), stringsAsFactors = FALSE)
    cat(sprintf("%-8s rep%d  err s/c = %.4f/%.4f  ratio=%.2f  maxShrink=%.3f\n",
        nm, rr, mean(eS), mean(eC), mean(eC) / mean(eS), max(shr))); flush.console()
  }
}
res <- do.call(rbind, res); saveRDS(res, "out/decisive.rds")
cat(sprintf("\n%-8s %9s %9s %7s %10s\n", "model", "errStock", "errCop", "ratio", "maxShrink"))
for (nm in c("iv2", "iv2rich")) {
  k <- res$model == nm; if (!any(k)) next
  cat(sprintf("%-8s %9.4f %9.4f %7.2f %10.3f\n", nm, mean(res$errStock[k]),
      mean(res$errCop[k]), mean(res$errCop[k]) / mean(res$errStock[k]), mean(res$maxShrink[k])))
  cat("  per-par stock :", paste(sprintf("%.3f", colMeans(do.call(rbind, res$perStock[k]))), collapse = " "), "\n")
  cat("  per-par copula:", paste(sprintf("%.3f", colMeans(do.call(rbind, res$perCop[k]))), collapse = " "), "\n")
}
