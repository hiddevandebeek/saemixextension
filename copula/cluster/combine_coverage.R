## combine_coverage.R -- pool the coverage replicates into a calibration table.
##
##   Rscript combine_coverage.R [out_coverage]
##
## Reads completed replicate files only, never the .part files a running chunk
## is still writing. The table answers one question per parameter: does the
## standard error the fit reports match the spread the estimator actually has.
## A ratio near one and coverage near 95% is the claim; anything else is the
## finding.
`%||%` <- function(a, b) if (is.null(a)) b else a
dir <- commandArgs(trailingOnly = TRUE)[1L] %||% "out_coverage"

files <- list.files(dir, "^cov_seed[0-9]+\\.rds$", full.names = TRUE)
if (!length(files)) stop("no completed replicates in ", dir)
z <- lapply(files, readRDS)
ok <- vapply(z, function(r) isTRUE(r$ok), TRUE)
cat(sprintf("replicates: %d found, %d ok, %d failed\n", length(z), sum(ok),
  sum(!ok)))
if (any(!ok)) {
  msg <- vapply(z[!ok], function(r) r$message %||% "?", character(1))
  cat("failure messages:\n"); print(table(substr(msg, 1L, 70L)))
}
z <- z[ok]
if (length(z) < 3L) stop("too few successful replicates to summarise")

## Did the averaging actually run? This study exists because it did not, so it
## is checked rather than assumed, and a mixed set is not pooled: averaged and
## unaveraged replicates are two different estimators.
avg <- vapply(z, function(r) as.integer(r$averagedOver %||% 0L), integer(1))
cat(sprintf("information averaged over %d to %d iterations (zero in %d)\n",
  min(avg), max(avg), sum(avg == 0L)))
if (any(avg == 0L) && any(avg > 0L))
  warning("some replicates averaged and some did not -- do NOT pool these")

pars <- Reduce(intersect, lapply(z, function(r) names(r$estimate)))
E <- do.call(rbind, lapply(z, function(r) r$estimate[pars]))
S <- do.call(rbind, lapply(z, function(r) r$se[pars]))
n <- nrow(E)

## The truth, on the reported scale. Only the coordinates whose truth is known
## are scored for coverage; the rest still get a calibration ratio, which needs
## no truth at all.
TR <- z[[1L]]$truth
known <- c(V = TR$V, CL = TR$CL,
  V.lognormal.sdlog = TR$sdlogV, CL.lognormal.sdlog = TR$sdlogCL,
  V.median = TR$V, CL.median = TR$CL,
  V.sd.log = TR$sdlogV, CL.sd.log = TR$sdlogCL,
  correlation.1.2 = TR$rho, residual.2 = TR$residual)

empirical <- apply(E, 2L, stats::sd)
reported <- colMeans(S, na.rm = TRUE)
## Positional lookup, not by name: the reported coordinates include derived
## quantities with no known truth, and indexing a named vector by an absent
## name is an error rather than an NA.
truth <- unname(known[match(pars, names(known))])
inside <- vapply(seq_along(pars), function(j) {
  if (is.na(truth[j])) return(NA_real_)
  mean(abs(E[, j] - truth[j]) < 1.959964 * S[, j], na.rm = TRUE)
}, numeric(1))

out <- data.frame(parameter = pars, truth = truth,
  mean = colMeans(E), empirical_sd = empirical, mean_se = reported,
  se_over_sd = round(reported / empirical, 3),
  coverage = round(inside, 3),
  bias_mcse = empirical / sqrt(n), row.names = NULL)
print(out, row.names = FALSE, digits = 4)

cat(sprintf("\n%d replicates. Monte Carlo error on a 0.95 coverage estimate is %.3f,\n",
  n, sqrt(.95 * .05 / n)))
cat(sprintf("and the SD ratio carries about %.0f%% error, so read it to that precision.\n",
  100 / sqrt(2 * (n - 1))))
cat("A ratio below about 0.8 means the reported errors understate the spread.\n")

saveRDS(list(n = n, table = out, estimates = E, se = S, averaged = avg,
  truth = TR, design = list(nsubj = z[[1L]]$nsubj, times = z[[1L]]$times,
    iters = z[[1L]]$iters)), file.path(dir, "pooled_coverage.rds"))
cat(sprintf("\nwritten: %s\n", file.path(dir, "pooled_coverage.rds")))
