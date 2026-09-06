## The paper's non-Gaussian eta demonstration, re-run.
##
## This is the second analysis: a paired study in which each replicate fits the
## conventional Gaussian-margin model, screens candidate parameter margins from
## that fit's posterior draws, and refits with whatever the screening chose. It
## is the study behind `data/natural-margin-selection/` and
## `figure2_natural_margin_200.png`, whose published counts over 200 replicates
## are 149 gamma, 17 Weibull, 34 lognormal and 27 retained standard, with the
## flexible arm's likelihood better in 163 and its gain exceeding twice the
## Monte Carlo error in 148.
##
## It is repeated because the estimator changed underneath those numbers: the
## information-averaging phase now actually runs, where before it never did, and
## the generalized gamma margin was rewritten. Point estimates are unaffected;
## the selection counts are not, because the retain-or-switch decision is a
## likelihood gain judged against a Monte Carlo error.
##
## The work is `ng_run_replicate` from the experiment's own functions.R,
## unchanged. This file only assigns replicates to chunks and points paths at a
## Linux box. Replicates seed themselves from their index (1060000 + 1000 * i),
## so chunks have disjoint seeds by construction, replicate i is the same
## simulated data as the published replicate i -- making this a paired
## comparison rather than an independent one -- and a rerun skips finished work.
##
## Environment: REPS, SEED0 (an index offset), TAG, OUTDIR, PKG or DEV,
## FUNCTIONS, NVPC.

suppressPackageStartupMessages({
  dev <- Sys.getenv("DEV", ""); pkg <- Sys.getenv("PKG", "")
  if (nzchar(dev)) devtools::load_all(dev, quiet = TRUE) else
    if (nzchar(pkg)) library(saemix, lib.loc = pkg) else library(saemix)
  library(rvinecopulib)
})
`%||%` <- function(x, y) if (is.null(x)) y else x

REPS  <- as.integer(Sys.getenv("REPS", "8"))
SEED0 <- as.integer(Sys.getenv("SEED0", "0"))
TAG   <- Sys.getenv("TAG", "c01")
OUT   <- Sys.getenv("OUTDIR", "out_eta_demonstration")
FUN   <- Sys.getenv("FUNCTIONS", "ng_functions.R")
NVPC  <- as.integer(Sys.getenv("NVPC", "100000"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(FUN)) stop("cannot find functions.R at ", FUN)
source(FUN)

cat(sprintf("[%s] eta demonstration: replicates %d..%d, nVpc = %d\n", TAG,
  SEED0 + 1L, SEED0 + REPS, NVPC))
for (k in seq_len(REPS)) {
  i <- SEED0 + k
  path <- file.path(OUT, sprintf("replicate_%03d.rds", i))
  if (file.exists(path)) { cat(sprintf("[%s] replicate %d already done\n",
    TAG, i)); next }
  started <- Sys.time()
  answer <- tryCatch({ ng_run_replicate(i, OUT, force = FALSE, nVpc = NVPC)
    "ok" }, error = function(e) paste("FAILED:", conditionMessage(e)))
  cat(sprintf("[%s] replicate %d %s (%.1f min)\n", TAG, i, answer,
    as.numeric(difftime(Sys.time(), started, units = "mins"))))
  flush.console()
}
cat(sprintf("[%s] chunk finished\n", TAG))
