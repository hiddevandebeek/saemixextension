## The candidate grid may be dispatched across mirai daemons. Candidate fitting
## is deterministic BFGS over pre-generated posterior pools with no RNG, so the
## parallel result must be bit-identical to the serial one.

suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet=TRUE)})
suppressPackageStartupMessages(library(rvinecopulib))

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-58s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}

if (!requireNamespace("mirai", quietly = TRUE)) {
  cat("mirai not installed; parallel dispatch test skipped\n"); quit(save = "no")
}

root <- "C:/package/saemix-copula/copula/experiments/combined-natural-frem-study"
source(file.path(root, "functions.R"))
simulation <- combined_simulate(1401001L, 100L)
sxdata <- combined_data(simulation$data)
standard <- saemix(combined_model(), sxdata, combined_control(1401011L, 500L),
  population = combined_population(combined_initial_vine(),
    combined_standard_margins(simulation$crp), simulation$crp,
    "transformed-additive"))
start <- combined_fit_start(standard); state <- copulaGet(standard)
incumbent <- saemix(combined_model(start$fixed, start$residual, start$omega),
  sxdata, combined_control(1401021L, 500L),
  population = combined_population(state$vine,
    c(state$margins[1:2], list(combined_flexible_covariate_margin(simulation$crp))),
    simulation$crp, "transformed-additive"))

arguments <- list(object = incumbent, supports = c("positive", "positive"),
  candidates = list(c("lognormal", "gamma", "weibull"),
                    c("lognormal", "gamma", "weibull")),
  posteriorDraws = 200L, max.iter = 200L, seed = 1401031L,
  optimizerMaxit = 200L, minimumEssFraction = 0,
  minimumPosteriorEssFraction = 0, maximumMcse = Inf,
  ## the product grid is pinned here so the dispatch check has a fixed,
  ## deterministic set of candidates; coordinate search is covered by
  ## test-margin-search-modes.R
  search = "product")

ok("no daemons means serial dispatch", !saemix:::copulaMiraiDaemons())
serial <- suppressWarnings(do.call(copulaSelectParameterMargins, arguments))
ok("grid is the full product", nrow(serial$table) == 9L)

mirai::daemons(4L, dispatcher = TRUE)
on.exit(try(mirai::daemons(0L), silent = TRUE), add = TRUE)
mirai::everywhere({
  suppressMessages({library(devtools)
    load_all("C:/package/saemix-copula", quiet = TRUE)})
})
ok("daemons detected", saemix:::copulaMiraiDaemons())
parallel <- suppressWarnings(do.call(copulaSelectParameterMargins, arguments))
mirai::daemons(0L)

ok("selected families identical",
  identical(serial$families, parallel$families),
  paste(parallel$families, collapse = "/"))
numeric <- c("training_delta_loglik", "validation_delta_loglik",
  "validation_bic_advantage", "validation_mcse",
  "minimum_ess_fraction", "median_ess_fraction")
difference <- max(abs(as.matrix(serial$table[, numeric]) -
  as.matrix(parallel$table[, numeric])))
ok("candidate table is bit-identical", difference == 0,
  sprintf("max diff %.0e", difference))
ok("full table identical", isTRUE(all.equal(serial$table, parallel$table,
  tolerance = 0)))

## the two quasi-Newton methods must agree on the selection
bfgs <- suppressWarnings(do.call(copulaSelectParameterMargins,
  c(arguments, list(optimiser = "BFGS"))))
lbfgs <- suppressWarnings(do.call(copulaSelectParameterMargins,
  c(arguments, list(optimiser = "L-BFGS-B"))))
ok("BFGS and L-BFGS-B select the same families",
  identical(bfgs$families, lbfgs$families),
  paste(lbfgs$families, collapse = "/"))
## The two searches stop at slightly different points inside the same basin, so
## the validation BIC advantage agrees to roughly 1e-4 rather than exactly. That
## is orders of magnitude below the Monte Carlo error gate the workflow already
## applies to this quantity (maximumMcse defaults to 0.5 BIC units), so it
## cannot change a ranking; the selected families are asserted separately.
gap <- max(abs(bfgs$table$validation_bic_advantage -
  lbfgs$table$validation_bic_advantage))
ok("attained BIC advantage agrees well inside MCSE", gap < 1e-2,
  sprintf("max diff %.2e", gap))
ok("L-BFGS-B is the default",
  identical(formals(copulaSelectParameterMargins)$optimiser, "L-BFGS-B"))
ok("unknown optimiser rejected", inherits(try(
  saemix:::copulaCandidateOptimise(c(0, 0), function(v) sum(v^2), 50L,
    "Nelder-Mead"), silent = TRUE), "try-error"))

## a single-candidate grid must never dispatch
one <- saemix:::copulaMapCandidates(1L, function(i) i * 2L)
ok("single candidate evaluated serially", identical(one, list(2L)))

if (nFail) stop(sprintf("candidate grid parallel checks failed: %d", nFail))
cat("candidate grid parallel dispatch checks passed\n")
