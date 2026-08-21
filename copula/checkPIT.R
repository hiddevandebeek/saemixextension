## Is the copula path's PIT actually misspecified?
## The copula block assumes eta has Gaussian margins, but the M-step sees
## POSTERIOR DRAWS -- a finite mixture over subjects, which need not be
## Gaussian-marginal even when the prior is.  Measure it directly.
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)})
setwd("C:/package/saemix-copula/copula"); source("R/simData.R")
vnG <- etaVineGaussian(PK_R)
set.seed(71); s <- simPK(120, vnG)
f <- saemix::saemix(pkSaemixModel(), pkSaemixData(s$data),
  list(seed = 1, save = FALSE, save.graphs = FALSE, print = FALSE,
       displayProgress = FALSE, nbiter.saemix = c(150, 100),
       nbiter.mcmc = c(2, 2, 2, 0), warnings = FALSE))
f <- saemix::conddist.saemix(f, nsamp = 20, plot = FALSE)
phi <- f@results@phi.samp; mp <- f@results@mean.phi
idx <- which(diag(f@results@omega) > 1e-10)
E <- do.call(rbind, lapply(seq_len(dim(phi)[3]), function(k) phi[, idx, k] - mp[, idx]))
sdE <- apply(E, 2, sd)
cat("POOLED conditional draws, n =", nrow(E), "\n")
cat(sprintf("%-4s %8s %8s %10s %12s\n", "eta", "skew", "exkurt", "SW p", "sd(hat)/sd(prior)"))
sdPrior <- sqrt(diag(f@results@omega))[idx]
for (j in seq_along(idx)) {
  x <- E[, j]; z <- (x - mean(x)) / sd(x)
  cat(sprintf("%-4d %8.3f %8.3f %10.2e %12.3f\n", j, mean(z^3), mean(z^4) - 3,
      shapiro.test(sample(x, min(4000, length(x))))$p.value, sdE[j] / sdPrior[j]))
}
## and per-SLICE (one draw per subject), which is what the M-step really averages
E1 <- phi[, idx, 1] - mp[, idx]
cat("\nSINGLE slice (one draw per subject), n =", nrow(E1), "\n")
for (j in seq_along(idx)) {
  x <- E1[, j]; z <- (x - mean(x)) / sd(x)
  cat(sprintf("%-4d %8.3f %8.3f %10.2e\n", j, mean(z^3), mean(z^4) - 3,
      shapiro.test(x)$p.value))
}
