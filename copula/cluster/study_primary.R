## study_primary.R -- the four-scenario primary study (400 frozen datasets,
## N = 100) and the ten N = 300 Gamma datasets, refitted with the current
## estimator so that the numbers in the paper come from the estimator the
## paper defines (audit item A1, 2026-09-05).
##
## Both arms use the likelihood-score recursion of `copulaPopulation`
## (`populationAlgorithm = "score-sa"`, three-phase gain, Fisher
## preconditioning, psi route). The Gaussian arm declares Normal covariate
## margins; the flexible arm declares the generating families and is warm
## started from the Gaussian fit, as in the joint study. Both are re-evaluated
## with an independent importance sample whose per-subject contributions are
## kept so that the non-nested (Vuong) test can be computed afterwards.
##
## Environment (launch.sh conventions):
##   REPS, SEED0   work items SEED0+1 .. SEED0+REPS of the work list
##                 (1..400 = manifest cells, 401..410 = N = 300 cells)
##   TAG, OUTDIR   chunk label and result directory
##   PKG | DEV     installed library or source tree
##   DATA          directory with the 400 frozen datasets and manifest.csv
##   DATA300       directory with the N = 300 datasets (optional)
##   ITERS         total iterations (default 1500), NMC importance draws
##                 (default 5000)
suppressPackageStartupMessages({
  dev <- Sys.getenv("DEV", ""); pkg <- Sys.getenv("PKG", "")
  if (nzchar(dev)) devtools::load_all(dev, quiet = TRUE) else
    if (nzchar(pkg)) library(saemix, lib.loc = pkg) else library(saemix)
  library(rvinecopulib)
})
`%||%` <- function(x, y) if (is.null(x)) y else x
S <- function(name) get(name, envir = asNamespace("saemix"))

REPS  <- as.integer(Sys.getenv("REPS", "10"))
SEED0 <- as.integer(Sys.getenv("SEED0", "0"))
TAG   <- Sys.getenv("TAG", "c01")
OUT   <- Sys.getenv("OUTDIR", "out_primary")
DATA  <- Sys.getenv("DATA", "primary_data")
DATA300 <- Sys.getenv("DATA300", "primary_data_n300")
ITERS <- as.integer(Sys.getenv("ITERS", "1500"))
NMC   <- as.integer(Sys.getenv("NMC", "5000"))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

manifest <- read.csv(file.path(DATA, "manifest.csv"), stringsAsFactors = FALSE)
work <- data.frame(cell = manifest$cell, path = file.path(DATA,
  paste0(manifest$cell, ".rds")), stringsAsFactors = FALSE)
if (dir.exists(DATA300)) {
  files <- sort(list.files(DATA300, "[.]rds$"))
  work <- rbind(work, data.frame(cell = sub("[.]rds$", "", files),
    path = file.path(DATA300, files), stringsAsFactors = FALSE))
}

pk <- function(psi, id, xidep) {
  V <- psi[id, 1L]; CL <- psi[id, 2L]
  xidep[, 1L] / V * exp(-(CL / V) * xidep[, 2L])
}
make_model <- function(fixed = c(V = 9, CL = 2.7), residual = .12,
                       omega = matrix(c(.0484, .003, .003, .0784), 2L))
  saemixModel(model = pk, modeltype = "structural",
    description = "primary study", psi0 = matrix(fixed, 1L,
      dimnames = list(NULL, c("V", "CL"))),
    transform.par = c(1, 1), covariance.model = matrix(1, 2L, 2L),
    omega.init = omega, error.model = "proportional",
    error.init = residual, verbose = FALSE)
make_population <- function(vine, margins, conditioning)
  copulaPopulation(vine, margins = margins, scale = "transformed-additive",
    conditioning = list(values = conditioning,
      variableName = c("COV1", "COV2")),
    populationAlgorithm = "score-sa", scoreScale = "auto",
    scoreBurn = 100L, scoreGainPower = .8, scoreFiniteDifference = 1e-4,
    scoreProjection = 24)
make_control <- function(seed, iterations = ITERS) list(seed = seed,
  save = FALSE, save.graphs = FALSE, print = FALSE, displayProgress = FALSE,
  warnings = FALSE, nbiter.saemix = c(iterations - 500L, 500L),
  nbiter.mcmc = c(2L, 2L, 2L, 0L), nmc.is = NMC, nu.is = 4L, ll.is = FALSE,
  fim = FALSE, map = FALSE)

state_ok <- function(state) {
  last <- state$lastJoint
  isTRUE(last$scoreTheory$runtimeConditionsObserved) &&
    identical(last$postFreezeProjectionCount, 0L) &&
    identical(last$postFreezeBacktrackCount, 0L) &&
    identical(last$postFreezeNoMoveCount, 0L) &&
    is.finite(last$scoreAverageMax)
}
covariate_loglik <- function(state, conditioning) {
  evaluated <- S("copulaGaussianFremEvaluateMargins")(conditioning,
    state$margins[3:4])
  if (!all(evaluated$valid)) return(NA_real_)
  R <- S("copulaGaussianRvineCor")(state$vine, 4L)[3:4, 3:4, drop = FALSE]
  sum(S("copulaGaussianLogDensity")(evaluated$z, R) -
    rowSums(stats::dnorm(evaluated$z, log = TRUE)) +
    rowSums(evaluated$logMargin))
}
trace_vector <- function(entry)
  c(entry$beta, entry$margin, entry$copula, entry$residual)
terminal_drift <- function(state, window = 500L) {
  raw <- state$trace[!vapply(state$trace, function(x) isTRUE(x$final),
    logical(1))]
  if (length(raw) < window) return(NA_real_)
  values <- do.call(rbind, lapply(tail(raw, window), trace_vector))
  half <- window %/% 2L
  previous <- colMeans(values[seq_len(half), , drop = FALSE])
  current <- colMeans(values[half + seq_len(window - half), , drop = FALSE])
  max(abs(current - previous) / pmax(abs(current), 1))
}
trace_frame <- function(state, arm, cell) {
  entries <- state$trace
  if (!length(entries)) return(NULL)
  keep <- unique(c(seq(1L, length(entries), by = 10L), length(entries)))
  rows <- lapply(entries[keep], function(entry) {
    head <- data.frame(arm = arm, cell = cell, kiter = entry$kiter,
      gamma = entry$gamma, value = entry$value, scoreMax = entry$scoreMax,
      residual = entry$residual)
    beta <- entry$beta
    if (length(beta)) names(beta) <- paste0("beta", seq_along(beta))
    margin <- entry$margin
    if (length(margin)) names(margin) <- paste0("margin_",
      if (is.null(names(margin))) seq_along(margin) else names(margin))
    copula <- entry$copula
    if (length(copula)) names(copula) <- paste0("copula", seq_along(copula))
    extra <- c(beta, margin, copula)
    if (!length(extra)) return(head)
    cbind(head, as.data.frame(as.list(extra)))
  })
  columns <- unique(unlist(lapply(rows, names)))
  rows <- lapply(rows, function(r) {
    for (missing in setdiff(columns, names(r))) r[[missing]] <- NA_real_
    r[, columns, drop = FALSE]
  })
  do.call(rbind, rows)
}

fit_arm <- function(arm, sim, sxdata, seed, llseed, start = NULL) {
  covFamilies <- if (arm == "standard") c("normal", "normal") else
    sim$covariateFamilies
  covariateMargins <- lapply(1:2, function(j)
    copulaFitCovariateMargin(sim$conditioning[, j], covFamilies[j]))
  if (is.null(start)) {
    margins <- c(list(copulaMarginNormal(.25), copulaMarginNormal(.30)),
      covariateMargins)
    R0 <- matrix(.04, 4L, 4L); diag(R0) <- 1
    vine <- copulaGaussianRvineFromCor(R0, sim$vineStructure)
    model <- make_model()
  } else {
    margins <- c(start$state$margins[1:2], covariateMargins)
    vine <- start$state$vine
    fixed <- as.numeric(start$fit@results@fixed.effects)[1:2]
    names(fixed) <- c("V", "CL")
    omega <- as.matrix(start$fit@results@omega[1:2, 1:2, drop = FALSE])
    if (any(!is.finite(omega)) ||
        inherits(try(chol(omega), silent = TRUE), "try-error"))
      omega <- diag(c(.25^2, .30^2))
    residual <- tail(as.numeric(start$fit@results@respar), 1L)
    model <- make_model(fixed, residual, omega)
  }
  elapsed <- system.time(fit <- saemix(model, sxdata,
    make_control(as.integer(seed)),
    population = make_population(vine, margins, sim$conditioning)))[["elapsed"]]
  llElapsed <- system.time(fit <- S("llisCopula.saemix")(fit, defensive = 0,
    batch = 100L, seed = as.integer(llseed)))[["elapsed"]]
  state <- copulaGet(fit); diagnostic <- attr(fit, "saemix.copula.likelihood")
  last <- state$lastJoint
  scales <- S("copulaMarginScales")(state$margins)
  R <- S("copulaGaussianRvineCor")(state$vine, 4L)
  D <- diag(c(scales[1:2], 1, 1), 4L)
  Sigma <- D %*% R %*% D
  B <- Sigma[1:2, 3:4, drop = FALSE] %*% solve(Sigma[3:4, 3:4, drop = FALSE])
  Omega <- Sigma[1:2, 1:2, drop = FALSE] - B %*% Sigma[3:4, 1:2, drop = FALSE]
  fixed <- as.numeric(fit@results@fixed.effects)[1:2]
  ll <- as.numeric(fit@results@ll.is); k <- as.integer(fit@results@npar.est)
  summary <- data.frame(cell = sim$cell, scenario = sim$scenario, arm = arm,
    V = fixed[1L], CL = fixed[2L], residual_sd = tail(fit@results@respar, 1L),
    eta_V_sd = scales[1L], eta_CL_sd = scales[2L],
    beta_cov1_V = B[1L, 1L], beta_cov2_V = B[1L, 2L],
    beta_cov1_CL = B[2L, 1L], beta_cov2_CL = B[2L, 2L],
    omega_V = Omega[1L, 1L], omega_V_CL = Omega[1L, 2L],
    omega_CL = Omega[2L, 2L], covariate_rho = cov2cor(Sigma[3:4, 3:4])[1L, 2L],
    log_likelihood_raw = ll, deviance_raw = -2 * ll, k = k,
    AIC_common = -2 * ll + 2 * k,
    BIC_N_common = -2 * ll + log(nrow(sim$conditioning)) * k,
    likelihood_mcse = diagnostic$se_loglik_total, ess_min = diagnostic$ess_min,
    covariate_log_likelihood = covariate_loglik(state, sim$conditioning),
    runtime_seconds = elapsed, likelihood_seconds = llElapsed,
    iterations = ITERS, backend = last$backend %||% NA_character_,
    score_average_max = last$scoreAverageMax %||% NA_real_,
    post_freeze_events = (last$postFreezeProjectionCount %||% NA_integer_) +
      (last$postFreezeBacktrackCount %||% NA_integer_) +
      (last$postFreezeNoMoveCount %||% NA_integer_),
    runtime_conditions_observed =
      isTRUE(last$scoreTheory$runtimeConditionsObserved),
    state_ok = state_ok(state), terminal_drift = terminal_drift(state),
    cov1_family = state$margins[[3L]]$name,
    cov2_family = state$margins[[4L]]$name,
    cov1_parameters = paste(names(state$margins[[3L]]$parameters),
      signif(state$margins[[3L]]$parameters, 12), collapse = ";"),
    cov2_parameters = paste(names(state$margins[[4L]]$parameters),
      signif(state$margins[[4L]]$parameters, 12), collapse = ";"),
    stringsAsFactors = FALSE)
  list(fit = fit, state = state, summary = summary, R = R, B = B,
    Omega = Omega,
    marginParameters = lapply(state$margins, `[[`, "parameters"),
    perSubjectLogLik = as.numeric(diagnostic$per_subject_loglik),
    trace = trace_frame(state, arm, sim$cell))
}

run_cell <- function(item) {
  path <- file.path(OUT, paste0(item$cell, ".rds"))
  if (file.exists(path)) return("already done")
  sim <- readRDS(item$path)
  ## The N = 300 Gamma records carry no manifest row: derive their seeds from
  ## their position in the work list and declare the families they were
  ## generated from (Gamma COV1, Normal COV2, as for gamma_normal).
  if (is.null(sim$seeds)) {
    index <- match(item$cell, work$cell)
    sim$seeds <- data.frame(data_seed = 1000000L + index * 100L + 1L,
      standard_seed = 1000000L + index * 100L + 12L,
      flexible_seed = 1000000L + index * 100L + 22L,
      standard_ll_seed = 1000000L + index * 100L + 32L,
      flexible_ll_seed = 1000000L + index * 100L + 42L)
  }
  if (is.null(sim$scenario)) sim$scenario <- "gamma_normal"
  if (is.null(sim$covariateFamilies))
    sim$covariateFamilies <- c("gamma", "normal")
  seeds <- sim$seeds
  sxdata <- saemixData(name.data = sim$data, header = TRUE, name.group = "id",
    name.predictors = c("dose", "time"), name.response = "y", verbose = FALSE)
  standard <- fit_arm("standard", sim, sxdata, seeds$standard_seed,
    seeds$standard_ll_seed)
  flexible <- fit_arm("flexible", sim, sxdata, seeds$flexible_seed,
    seeds$flexible_ll_seed, start = standard)
  strip <- function(arm) arm[setdiff(names(arm), c("fit", "state"))]
  result <- list(status = "complete", schema = 2L, cell = sim$cell,
    scenario = sim$scenario, N = nrow(sim$conditioning),
    covariateFamilies = sim$covariateFamilies, truth = sim$truth,
    truthR = sim$truth$R, seeds = seeds, iterations = ITERS,
    importanceDraws = NMC, standard = strip(standard),
    flexible = strip(flexible), completed = Sys.time(),
    rVersion = R.version.string)
  temporary <- paste0(path, ".tmp-", Sys.getpid()); saveRDS(result, temporary)
  if (!file.rename(temporary, path)) stop("atomic rename failed")
  sprintf("ok (%.0f + %.0f s, dLL %.2f)", standard$summary$runtime_seconds,
    flexible$summary$runtime_seconds,
    flexible$summary$log_likelihood_raw - standard$summary$log_likelihood_raw)
}

cat(sprintf("[%s] primary study: work items %d..%d of %d, %d iterations, %d draws\n",
  TAG, SEED0 + 1L, SEED0 + REPS, nrow(work), ITERS, NMC))
for (k in seq_len(REPS)) {
  i <- SEED0 + k
  if (i > nrow(work)) break
  started <- Sys.time()
  answer <- tryCatch(run_cell(work[i, ]), error = function(e)
    paste("FAILED:", conditionMessage(e)))
  cat(sprintf("[%s] %s %s (%.1f min)\n", TAG, work$cell[i], answer,
    as.numeric(difftime(Sys.time(), started, units = "mins"))))
  flush.console()
}
cat(sprintf("[%s] chunk finished\n", TAG))
