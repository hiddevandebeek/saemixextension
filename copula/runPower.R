## Power study: is a non-Gaussian eta dependence detectable from a standard
## Gaussian-prior SAEM fit?  Null calibrated by simulation, not by chi-square.
suppressMessages({library(devtools); load_all("..", quiet = TRUE)})
source("R/replicate.R")

args <- commandArgs(trailingOnly = TRUE)
NREP <- if (length(args) > 0) as.integer(args[1]) else 60
N    <- if (length(args) > 1) as.integer(args[2]) else 100
FAM  <- if (length(args) > 2) args[3] else "gumbel"
TAG  <- if (length(args) > 3) args[4] else sprintf("N%d_%s", N, FAM)

tw <- makeTwinVines(PK_R, FAM)
res <- list(); k <- 0L
for (arm in c("gauss", "alt")) for (r in seq_len(NREP)) {
  k <- k + 1L
  sd <- if (arm == "gauss") 10000L else 20000L
  out <- try(suppressWarnings(runReplicate(sd + r, arm, tw, N = N, nsamp = 3,
                                           nbiter = c(150, 80))), silent = TRUE)
  if (!inherits(out, "try-error")) res[[length(res) + 1]] <- out
  cat(sprintf("[%s] %s rep %d/%d\n", format(Sys.time(), "%H:%M:%S"), arm, r, NREP))
  flush.console()
}
res <- do.call(rbind, res)
saveRDS(res, file.path("out", paste0("power_", TAG, ".rds")))
cat("saved out/power_", TAG, ".rds  rows=", nrow(res), "\n", sep = "")
