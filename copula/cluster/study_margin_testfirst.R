## The proposed test-first workflow, evaluated end to end.
##
## The published 500-replicate study (study_margin_selection.R) ranks candidate
## clearance families straight from the incumbent's posterior draws and holds
## volume lognormal. Both reviews of 7 September 2026 pointed out that this is
## the ranking-and-refit route rather than the workflow the paper proposes, in
## which the score test of the latent margins runs first and a coordinate is
## offered to the ranking only where a direction rejects.
##
## This chunk driver runs that workflow on the same datasets: the generator,
## the seeds, the standard arm and every reported endpoint are those of
## study_margin_selection.R, so replicate i of the two studies is the same
## dataset and the same Gaussian comparator, and the two workflows can be
## compared directly.
##
## Environment:
##   REPS      replicates in this chunk               (default 10)
##   SEED0     chunk offset; covers replicate SEED0+1 .. SEED0+REPS
##   TAG       chunk label, for messages
##   OUTDIR    where replicate_NNN.rds go
##   PKG       library path for the installed package
##   DEV       source tree for devtools, for local smoke tests
##   FUNCTIONS path to the experiment's functions.R
##   TESTFUNCTIONS path to functions_testfirst.R
##   MODE      "final" (250 subjects) or "pilot"
##   ITERS     override the iteration count
##   BOOTSTRAP bootstrap replicates for the score test (default 59)

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
OUT   <- Sys.getenv("OUTDIR", "out_margin_testfirst")
MODE  <- Sys.getenv("MODE", "final")
FUN   <- Sys.getenv("FUNCTIONS", "functions.R")
TFUN  <- Sys.getenv("TESTFUNCTIONS", "functions_testfirst.R")
BOOT  <- as.integer(Sys.getenv("BOOTSTRAP", "59"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(FUN)) stop("cannot find functions.R at ", FUN)
if (!file.exists(TFUN)) stop("cannot find functions_testfirst.R at ", TFUN)
source(FUN); source(TFUN)

settings <- if (identical(MODE, "pilot"))
  list(nSubjects = 180L, nVpc = 20000L, iterations = 1000L,
    screenDraws = 400L, screenRetry = 1000L, likelihoodDraws = 1200L) else
  list(nSubjects = 250L, nVpc = 50000L, iterations = 1500L,
    screenDraws = 700L, screenRetry = 2000L, likelihoodDraws = 2500L)
if (nzchar(Sys.getenv("ITERS", "")))
  settings$iterations <- as.integer(Sys.getenv("ITERS"))
settings$bootstrap <- BOOT

cat(sprintf("[%s] test-first workflow (%s): replicates %d..%d\n", TAG, MODE,
  SEED0 + 1L, SEED0 + REPS))
cat(sprintf("[%s] %d subjects, %d iterations, %d screening draws, %d bootstrap\n",
  TAG, settings$nSubjects, settings$iterations, settings$screenDraws, BOOT))

for (k in seq_len(REPS)) {
  i <- SEED0 + k
  path <- file.path(OUT, sprintf("replicate_%03d.rds", i))
  if (file.exists(path)) { cat(sprintf("[%s] replicate %d already done\n",
    TAG, i)); next }
  started <- Sys.time()
  answer <- tryCatch({
    do.call(combined_run_replicate_test_first, c(list(replicate = i,
      resultRoot = OUT, force = FALSE), settings))
    "ok"
  }, error = function(e) paste("FAILED:", conditionMessage(e)))
  cat(sprintf("[%s] replicate %d %s (%.1f min)\n", TAG, i, answer,
    as.numeric(difftime(Sys.time(), started, units = "mins"))))
  flush.console()
}
cat(sprintf("[%s] chunk finished\n", TAG))
