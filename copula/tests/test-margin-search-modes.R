## Family search strategy. The product grid enumerates K^d combinations, which
## is unusable beyond a couple of eta coordinates; coordinate ascent varies one
## coordinate at a time, K*d per sweep. Every comparison is the same exact
## likelihood-ratio identity, so the two must agree on the families they
## return -- only the set of visited combinations differs.

suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet=TRUE)})
suppressPackageStartupMessages(library(rvinecopulib))

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-56s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}

## --- the option helper validates support, the grid is its product ---
opts <- saemix:::copulaNaturalFamilyOptions(NULL, c("positive", "positive"))
ok("options are per coordinate", length(opts) == 2L &&
  all(vapply(opts, function(o) "gamma" %in% o, logical(1))))
grid <- saemix:::copulaNaturalFamilyGrid(NULL, c("positive", "positive"), 81L)
ok("grid is the product of the options",
  nrow(grid) == prod(lengths(opts)), sprintf("%d rows", nrow(grid)))
ok("grid still refuses an oversized product", inherits(try(
  saemix:::copulaNaturalFamilyGrid(NULL, rep("positive", 5L), 81L),
  silent = TRUE), "try-error"))
ok("incompatible family still rejected", inherits(try(
  saemix:::copulaNaturalFamilyOptions(list("gamma"), "unit"),
  silent = TRUE), "try-error"))

## --- coordinate ascent visits far fewer combinations as d grows ---
perSweep <- function(K, d) K + (d - 1L) * (K - 1L)
for (d in c(2L, 3L, 5L)) {
  ok(sprintf("d=%d: coordinate %d per sweep vs product %d", d,
    perSweep(3L, d), 3L^d), perSweep(3L, d) < 3L^d)
}

## --- the two searches agree on a real fit ---
root <- "C:/package/saemix-copula/copula/experiments/combined-natural-frem-study"
source(file.path(root, "functions.R"))
simulation <- combined_simulate(1401001L, 200L)
sxdata <- combined_data(simulation$data)
cov <- combined_flexible_covariate_margin(simulation$crp)
standard <- saemix(combined_model(), sxdata, combined_control(1401011L, 500L),
  population = combined_population(combined_initial_vine(),
    combined_standard_margins(simulation$crp), simulation$crp,
    "transformed-additive"))
st <- copulaGet(standard); s <- combined_fit_start(standard)
incumbent <- saemix(combined_model(s$fixed, s$residual, s$omega), sxdata,
  combined_control(1401021L, 500L),
  population = combined_population(st$vine,
    c(st$margins[1:2], list(cov)), simulation$crp, "transformed-additive"))
args <- list(object = incumbent, supports = c("positive", "positive"),
  candidates = list(c("lognormal", "gamma", "weibull"),
    c("lognormal", "gamma", "weibull")),
  posteriorDraws = 400L, max.iter = 300L, seed = 1401031L,
  optimizerMaxit = 200L, minimumEssFraction = 0,
  minimumPosteriorEssFraction = 0, maximumMcse = Inf)
product <- suppressWarnings(do.call(copulaSelectParameterMargins,
  c(args, list(search = "product"))))
coordinate <- suppressWarnings(do.call(copulaSelectParameterMargins,
  c(args, list(search = "coordinate"))))
ok("both searches select the same families",
  identical(product$families, coordinate$families),
  paste(coordinate$families, collapse = "/"))
gap <- max(product$table$validation_bic_advantage) -
  max(coordinate$table$validation_bic_advantage)
ok("coordinate ascent loses no BIC advantage", gap <= 1e-8,
  sprintf("gap %.2e", gap))
ok("coordinate ascent evaluates fewer combinations",
  nrow(coordinate$table) <= nrow(product$table),
  sprintf("%d vs %d", nrow(coordinate$table), nrow(product$table)))
ok("coordinate is the default",
  identical(eval(formals(copulaSelectParameterMargins)$search)[1], "coordinate"))
ok("Pareto k is reported in the table",
  "maximum_pareto_k" %in% names(coordinate$table))

if (nFail) stop(sprintf("margin search checks failed: %d", nFail))
cat("margin search mode checks passed\n")
