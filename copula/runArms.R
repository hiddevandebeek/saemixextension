## Can the copula path match stock saemix in the NESTED NULL (Gaussian vine)?
## The null is the hardest case for it: the copula buys nothing and pays every
## cost, so any gap here is pure estimator inefficiency.  Arms:
##   pool-v1   Q-hat as a truncated weighted particle pool, sd from the pool
##   pool-SS   same, but the marginal sd taken from saemix's EXACT SA sufficient
##             statistic (Gaussian margins still have one even when the copula
##             does not) -- only the copula term uses the pool
##   sa        no pool at all: fit the vine to the CURRENT draws, then a
##             Robbins-Monro step on the pair-copula parameters in tau space
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)})
setwd("C:/package/saemix-copula/copula"); source("R/simData.R")
NREP <- 5
vnG <- etaVineGaussian(PK_R); truth <- c(PK_TRUE, PK_SD)
ctl <- function(sd) list(seed = sd, save = FALSE, save.graphs = FALSE, print = FALSE,
                         displayProgress = FALSE, nbiter.saemix = c(150, 120),
                         nbiter.mcmc = c(2, 2, 2, 0), warnings = FALSE)
arms <- list(
  stock     = NULL,
  `pool-v1` = list(mode = "pool", sdFromSS = FALSE),
  `pool-SS` = list(mode = "pool", sdFromSS = TRUE),
  `sa`      = list(mode = "sa",   sdFromSS = TRUE))
res <- list()
for (rr in seq_len(NREP)) {
  set.seed(70 + rr); sD <- simPK(120, vnG)
  dat <- pkSaemixData(sD$data)
  for (a in names(arms)) {
    copulaClear()
    if (!is.null(arms[[a]]))
      copulaSet(vnG, PK_SD, familySet = "gaussian", mode = arms[[a]]$mode,
                sdFromSS = arms[[a]]$sdFromSS)
    f <- try(saemix::saemix(pkSaemixModel(), dat, ctl(rr)), silent = TRUE)
    if (inherits(f, "try-error")) { cat(rr, a, "FAILED\n"); copulaClear(); next }
    est <- c(f@results@fixed.effects,
             if (is.null(arms[[a]])) sqrt(diag(f@results@omega)) else copulaGet()$sd)
    copulaClear()
    res[[length(res) + 1]] <- data.frame(rep = rr, arm = a, est = I(list(est)))
  }
  cat("rep", rr, "done\n"); flush.console()
}
res <- do.call(rbind, res); saveRDS(res, "out/arms.rds")
cat(sprintf("\n%-9s %-47s %s\n", "arm", "mean |rel err| : V1 CL Q V2 sdV1 sdCL sdQ sdV2", "overall"))
for (a in names(arms)) {
  M <- do.call(rbind, res$est[res$arm == a])
  if (!length(M)) next
  R <- abs(sweep(M, 2, truth, "-") / rep(truth, each = nrow(M)))
  cat(sprintf("%-9s %-47s %.4f\n", a,
      paste(sprintf("%.3f", colMeans(R)), collapse = " "), mean(R)))
}
