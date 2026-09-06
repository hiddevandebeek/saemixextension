## Interval coverage for the score recursion, re-measured.
##
## The published 96% was obtained while the information-averaging phase never
## ran, so every standard error in it came from a single noisy Istar rather
## than an average over the tail. That is now fixed, and the figure has to be
## re-established before it can be quoted.
##
## One replicate per file, named by seed, so a chunk that dies loses only the
## replicate it was in the middle of and a rerun skips what is already there.
## Nothing here reads a file another chunk writes.
##
## Environment:
##   REPS      replicates in this chunk           (default 10)
##   SEED0     first seed; the chunk covers SEED0+1 .. SEED0+REPS
##   TAG       chunk label, used only in messages
##   NSUBJ     subjects                           (default 400)
##   TIMES     comma-separated sampling times     (default "2,24")
##   ITERS     total saemix iterations            (default 1500)
##   OUTDIR    where the per-replicate rds go     (default "out_coverage")
##   PKG       library path for the installed package, or "" for the default

## DEV points at a source tree and uses devtools, which is how this is smoke
## tested locally before deploying; PKG points at an installed library, which
## is how it runs on the cluster. Neither set falls back to the default library.
suppressPackageStartupMessages({
  dev <- Sys.getenv("DEV", "")
  pkg <- Sys.getenv("PKG", "")
  if (nzchar(dev)) devtools::load_all(dev, quiet = TRUE) else
    if (nzchar(pkg)) library(saemix, lib.loc = pkg) else library(saemix)
  library(rvinecopulib)
})
`%||%` <- function(x, y) if (is.null(x)) y else x

REPS  <- as.integer(Sys.getenv("REPS", "10"))
SEED0 <- as.integer(Sys.getenv("SEED0", "1000"))
TAG   <- Sys.getenv("TAG", "c01")
NSUBJ <- as.integer(Sys.getenv("NSUBJ", "400"))
ITERS <- as.integer(Sys.getenv("ITERS", "1500"))
OUT   <- Sys.getenv("OUTDIR", "out_coverage")
TIMES <- as.numeric(strsplit(Sys.getenv("TIMES", "2,24"), ",")[[1L]])
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## The truth. Fixed here rather than read from anywhere, so a replicate file is
## self-describing and two chunks cannot disagree about what they simulated.
TRUTH <- list(V = 10, CL = 3.2, sdlogV = .22, sdlogCL = .30, rho = .2,
  residual = .12, dose = 100)

oneCompartment <- function(psi, id, xidep) {
  V <- psi[id, 1L]; CL <- psi[id, 2L]
  dose <- xidep[, 1L]; t <- xidep[, 2L]
  dose / V * exp(-(CL / V) * t)
}

model <- function() saemixModel(model = oneCompartment,
  modeltype = "structural", description = "one-compartment IV bolus",
  psi0 = matrix(c(TRUTH$V, TRUTH$CL), 1L,
    dimnames = list(NULL, c("V", "CL"))),
  transform.par = c(1, 1), covariance.model = diag(1, 2L),
  omega.init = diag(c(TRUTH$sdlogV, TRUTH$sdlogCL)^2),
  error.model = "proportional", error.init = c(0, TRUTH$residual),
  verbose = FALSE)

control <- function(seed) list(seed = seed, save = FALSE, save.graphs = FALSE,
  print = FALSE, displayProgress = FALSE, warnings = FALSE,
  nbiter.saemix = c(ITERS - 500L, 500L), nbiter.mcmc = c(2, 2, 2, 0),
  ll.is = FALSE, fim = FALSE, map = FALSE)

simulate <- function(seed) {
  set.seed(seed)
  ## Correlated etas on the log scale, so the copula correlation is a real
  ## parameter with a known truth rather than a nuisance fixed at zero.
  z <- matrix(rnorm(2L * NSUBJ), NSUBJ, 2L)
  z[, 2L] <- TRUTH$rho * z[, 1L] + sqrt(1 - TRUTH$rho^2) * z[, 2L]
  psi <- cbind(TRUTH$V * exp(TRUTH$sdlogV * z[, 1L]),
               TRUTH$CL * exp(TRUTH$sdlogCL * z[, 2L]))
  dd <- data.frame(id = rep(seq_len(NSUBJ), each = length(TIMES)),
    dose = TRUTH$dose, time = rep(TIMES, NSUBJ))
  dd$y <- pmax(oneCompartment(psi, dd$id, cbind(dd$dose, dd$time)) *
    (1 + TRUTH$residual * rnorm(nrow(dd))), 1e-8)
  saemixData(name.data = dd, header = TRUE, name.group = "id",
    name.predictors = c("dose", "time"), name.response = "y", verbose = FALSE)
}

runOne <- function(seed) {
  data <- simulate(seed)
  vine <- copulaGaussianRvineFromCor(
    matrix(c(1, TRUTH$rho, TRUTH$rho, 1), 2L, 2L),
    rvinecopulib::dvine_structure(1:2))
  started <- Sys.time()
  fit <- saemix(model(), data, control(seed + 500000L),
    population = copulaPopulation(vine, margins = list(
        copulaNaturalMarginLognormal(TRUTH$sdlogV),
        copulaNaturalMarginLognormal(TRUTH$sdlogCL)),
      scale = "parameter", populationAlgorithm = "score-sa",
      scoreScale = "auto", scoreBurn = 50L))
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
  tab <- copulaStandardErrors(fit)
  state <- copulaGet(fit)
  list(seed = seed, ok = TRUE, elapsed = elapsed,
    estimate = stats::setNames(tab$estimate, tab$parameter),
    se = stats::setNames(tab$se, tab$parameter),
    ## Whether the averaging actually ran, and over how many iterations. This
    ## is the whole reason the study is being repeated, so it is recorded per
    ## replicate rather than assumed.
    averagedOver = state$scoreState$fisherAverageCount %||% 0L,
    smoothingSource = state$scoreSmoothingSource %||% NA_character_,
    truth = TRUTH, nsubj = NSUBJ, times = TIMES, iters = ITERS)
}

cat(sprintf("[%s] coverage: %d replicates, seeds %d..%d, n=%d, %d obs, %d iter\n",
  TAG, REPS, SEED0 + 1L, SEED0 + REPS, NSUBJ, length(TIMES), ITERS))
for (k in seq_len(REPS)) {
  seed <- SEED0 + k
  path <- file.path(OUT, sprintf("cov_seed%06d.rds", seed))
  if (file.exists(path)) { cat(sprintf("[%s] seed %d already done\n", TAG, seed))
    next }
  answer <- tryCatch(runOne(seed), error = function(e)
    list(seed = seed, ok = FALSE, message = conditionMessage(e)))
  ## Write to a temporary name and rename, so a reader can never see a
  ## half-written file and a killed chunk leaves no corrupt result behind.
  tmp <- paste0(path, ".part")
  saveRDS(answer, tmp); file.rename(tmp, path)
  cat(sprintf("[%s] seed %d %s%s\n", TAG, seed,
    if (isTRUE(answer$ok)) "ok" else paste("FAILED:", answer$message),
    if (isTRUE(answer$ok)) sprintf(" (%.0fs, averaged over %d)",
      answer$elapsed, answer$averagedOver) else ""))
  flush.console()
}
cat(sprintf("[%s] chunk finished\n", TAG))
