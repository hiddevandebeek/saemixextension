## Memoised hot-path quantities must return exactly what recomputing returns,
## and must invalidate when anything they depend on changes. A stale cache here
## would be silent, so each cache is checked for both hit and miss behaviour.

suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet=TRUE)})
suppressPackageStartupMessages(library(rvinecopulib))

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-56s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}
set.seed(8021L)

## ---- covariate probability integral transform ----
N <- 300L
crp <- rgamma(N, 2, scale = 20)
cond <- matrix(crp, ncol = 1L)
gammaMargin <- list(copulaFitCovariateMargin(crp, "gamma"))
weibullMargin <- list(copulaFitCovariateMargin(crp, "weibull"))

reference <- saemix:::copulaGaussianFremEvaluateMargins(cond, gammaMargin)
first <- saemix:::copulaCovariateTransform(cond, gammaMargin)
second <- saemix:::copulaCovariateTransform(cond, gammaMargin)
ok("covariate transform matches direct evaluation", identical(first, reference))
ok("repeat call returns the same object", identical(second, reference))

other <- saemix:::copulaCovariateTransform(cond, weibullMargin)
ok("different family invalidates the cache",
  identical(other, saemix:::copulaGaussianFremEvaluateMargins(cond,
    weibullMargin)) && !identical(other$z, reference$z))

bumped <- gammaMargin
bumped[[1L]]$parameters[1L] <- bumped[[1L]]$parameters[1L] * 1.05
ok("changed parameters invalidate the cache",
  identical(saemix:::copulaCovariateTransform(cond, bumped),
    saemix:::copulaGaussianFremEvaluateMargins(cond, bumped)))

cond2 <- matrix(rgamma(N, 2, scale = 20), ncol = 1L)
ok("changed covariate values invalidate the cache",
  identical(saemix:::copulaCovariateTransform(cond2, gammaMargin),
    saemix:::copulaGaussianFremEvaluateMargins(cond2, gammaMargin)))
ok("original key still reproduces the original answer",
  identical(saemix:::copulaCovariateTransform(cond, gammaMargin), reference))

## ---- vine correlation decode ----
d <- 4L
A <- matrix(rnorm(d * d), d); R <- cov2cor(crossprod(A) + diag(d))
vine <- copulaGaussianRvineFromCor(R, rvinecopulib::dvine_structure(seq_len(d)))
ok("decode matches the correlation it was built from",
  max(abs(copulaGaussianRvineCor(vine, d) - R)) < 1e-10)
ok("cached decode is stable across repeats",
  identical(copulaGaussianRvineCor(vine, d), copulaGaussianRvineCor(vine, d)))
B <- matrix(rnorm(d * d), d); R2 <- cov2cor(crossprod(B) + diag(d))
vine2 <- copulaGaussianRvineFromCor(R2, rvinecopulib::dvine_structure(seq_len(d)))
ok("a different vine is not served from the cache",
  max(abs(copulaGaussianRvineCor(vine2, d) - R2)) < 1e-10)
ok("returning to the first vine still decodes correctly",
  max(abs(copulaGaussianRvineCor(vine, d) - R)) < 1e-10)

## ---- response likelihood layout ----
root <- "C:/package/saemix-copula/copula/experiments/combined-natural-frem-study"
source(file.path(root, "functions.R"))
simulation <- combined_simulate(1401001L, 120L)
fit <- saemix(combined_model(), combined_data(simulation$data),
  combined_control(1401011L, 500L),
  population = combined_population(combined_initial_vine(),
    combined_standard_margins(simulation$crp), simulation$crp,
    "transformed-additive"))
phi <- fit["results"]["mean.phi"]
layout <- copulaResponseLayout(fit, 1L)
worst <- 0
for (r in 1:5) {
  p <- phi + matrix(rnorm(length(phi), 0, .2), nrow(phi))
  worst <- max(worst, max(abs(
    copulaResponseLogLikBatch(fit, p, 1L) -
    copulaResponseLogLikBatch(fit, p, 1L, layout))))
}
ok("reused layout gives identical log-likelihood", worst == 0,
  sprintf("max diff %.0e", worst))

if (nFail) stop(sprintf("hot-path cache checks failed: %d", nFail))
cat("hot-path cache checks passed\n")
