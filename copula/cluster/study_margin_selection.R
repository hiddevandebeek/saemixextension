## The paper's margin-selection study, re-run.
##
## This is the second analysis: fit the conventional model, draw from its
## posterior, and choose the parameter margin from those draws -- the workflow
## a user actually follows. It is the study behind
## `data/natural-margin-selection/` and `figure2_natural_margin_200.png`, whose
## published counts are 149 gamma, 17 Weibull, 34 lognormal and 27 retained
## standard out of 200 replicates.
##
## It is being repeated because the estimator has changed since those numbers
## were produced: the information-averaging phase now actually runs, and the
## generalized gamma margin has been rewritten. Point estimates are unaffected,
## but anything downstream of a standard error is.
##
## The work itself is `combined_run_replicate` from the experiment's own
## functions.R, unchanged -- this file only assigns replicates to chunks and
## points the paths somewhere that exists on a Linux box. `combined_run_replicate`
## is already resumable and seeds itself from the replicate index, so chunks
## have disjoint seeds by construction and a rerun skips completed work.
##
## Environment:
##   REPS      replicates in this chunk               (default 10)
##   SEED0     chunk offset; covers replicate SEED0+1 .. SEED0+REPS
##             (an index offset, not a seed -- the seeds derive from the index)
##   TAG       chunk label, for messages
##   OUTDIR    where replicate_NNN.rds go
##   PKG       library path for the installed package
##   DEV       source tree for devtools, for local smoke tests
##   FUNCTIONS path to the experiment's functions.R
##   MODE      "final" (250 subjects, 1500 iterations) or "pilot"

suppressPackageStartupMessages({
  dev <- Sys.getenv("DEV", ""); pkg <- Sys.getenv("PKG", "")
  if (nzchar(dev)) devtools::load_all(dev, quiet = TRUE) else
    if (nzchar(pkg)) library(saemix, lib.loc = pkg) else library(saemix)
  library(rvinecopulib)
})
`%||%` <- function(x, y) if (is.null(x)) y else x

REPS  <- as.integer(Sys.getenv("REPS", "10"))
SEED0 <- as.integer(Sys.getenv("SEED0", "0"))
TAG   <- Sys.getenv("TAG", "c01")
OUT   <- Sys.getenv("OUTDIR", "out_margin_selection")
MODE  <- Sys.getenv("MODE", "final")
FUN   <- Sys.getenv("FUNCTIONS", "functions.R")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(FUN)) stop("cannot find functions.R at ", FUN)
source(FUN)

settings <- if (identical(MODE, "pilot"))
  list(nSubjects = 180L, nVpc = 20000L, iterations = 1000L,
    screenDraws = 400L, screenRetry = 1000L, likelihoodDraws = 1200L) else
  list(nSubjects = 250L, nVpc = 50000L, iterations = 1500L,
    screenDraws = 700L, screenRetry = 2000L, likelihoodDraws = 2500L)
## Iterations are overridable because 1500 is a floor rather than a margin: at
## 800 the agreement test fails outright, and 1500 is the smallest count at
## which this estimator has been seen to converge. Nothing else about the gain
## needs touching to raise it -- the preheat length is defined as total/5, so
## the ratio that matters stays at 0.2 whatever the total is, and a preheat too
## long relative to the run is what previously cost coverage (1000 of 1500 gave
## 53%, 200 of 1500 gave 96%).
if (nzchar(Sys.getenv("ITERS", "")))
  settings$iterations <- as.integer(Sys.getenv("ITERS"))

cat(sprintf("[%s] margin selection (%s): replicates %d..%d\n", TAG, MODE,
  SEED0 + 1L, SEED0 + REPS))
cat(sprintf("[%s] %d subjects, %d iterations, %d screening draws\n", TAG,
  settings$nSubjects, settings$iterations, settings$screenDraws))

for (k in seq_len(REPS)) {
  i <- SEED0 + k
  path <- file.path(OUT, sprintf("replicate_%03d.rds", i))
  if (file.exists(path)) { cat(sprintf("[%s] replicate %d already done\n",
    TAG, i)); next }
  started <- Sys.time()
  answer <- tryCatch({
    do.call(combined_run_replicate, c(list(replicate = i, resultRoot = OUT,
      force = FALSE), settings))
    "ok"
  }, error = function(e) paste("FAILED:", conditionMessage(e)))
  cat(sprintf("[%s] replicate %d %s (%.1f min)\n", TAG, i, answer,
    as.numeric(difftime(Sys.time(), started, units = "mins"))))
  flush.console()
}
cat(sprintf("[%s] chunk finished\n", TAG))
