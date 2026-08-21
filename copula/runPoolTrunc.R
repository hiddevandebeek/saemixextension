## Is the nested-null gap caused by POOL TRUNCATION?
##
## Stock saemix carries Q-hat as a sufficient statistic, an exact exponentially
## weighted average over ALL past iterations.  The copula path must carry Q-hat
## as a particle pool, which has to be capped -- so at equal iteration count it
## averages over fewer draws and is noisier.  If that is the mechanism, raising
## poolMax should close the gap.  If it does not, the gap is the IFM M-step.
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)})
setwd("C:/package/saemix-copula/copula"); source("R/simData.R")
vnG <- etaVineGaussian(PK_R)
ctl <- function(sd) list(seed = sd, save = FALSE, save.graphs = FALSE, print = FALSE,
                         displayProgress = FALSE, nbiter.saemix = c(150, 80),
                         nbiter.mcmc = c(2, 2, 2, 0), warnings = FALSE)
truth <- c(PK_TRUE, PK_SD)
res <- list()
for (rr in 1:4) {
  set.seed(70 + rr); sD <- simPK(120, vnG)
  copulaClear()
  fS <- saemix::saemix(pkSaemixModel(), pkSaemixData(sD$data), ctl(rr))
  res[[length(res)+1]] <- data.frame(rep=rr, arm="stock",
    est=I(list(c(fS@results@fixed.effects, sqrt(diag(fS@results@omega))))))
  for (pm in c(40L, 250L)) {
    copulaSet(vnG, PK_SD, familySet = "gaussian", poolMax = pm)
    fC <- saemix::saemix(pkSaemixModel(), pkSaemixData(sD$data), ctl(rr))
    res[[length(res)+1]] <- data.frame(rep=rr, arm=paste0("cop", pm),
      est=I(list(c(fC@results@fixed.effects, copulaGet()$sd))))
    copulaClear()
  }
  cat("rep", rr, "done\n"); flush.console()
}
res <- do.call(rbind, res)
saveRDS(res, "out/poolTrunc.rds")
cat(sprintf("\n%-8s %s\n", "arm", "mean |relative error vs truth| per parameter"))
for (a in unique(res$arm)) {
  M <- do.call(rbind, res$est[res$arm == a])
  cat(sprintf("%-8s %s   overall=%.4f\n", a,
      paste(sprintf("%.3f", colMeans(abs(sweep(M, 2, truth, "-") / rep(truth, each = nrow(M))))), collapse=" "),
      mean(abs(sweep(M, 2, truth, "-") / rep(truth, each = nrow(M))))))
}
