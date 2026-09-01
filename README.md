# Flexible-margin Gaussian-copula FREM for saemix

This fork extends `saemix` with a fixed Gaussian-copula full random-effects
model (FREM) whose random-effect and covariate margins can be declared
separately from their dependence.

The supported estimator uses controlled-MCMC likelihood-score stochastic
approximation for one fixed model specification. The all-Normal member is
exactly ordinary Gaussian FREM. Calling `saemix()` without `population`, or
with `population = NULL`, retains the stock Gaussian workflow.

## Minimal example

```r
population <- gaussianCopulaFrem(
  etaSd = c(0.25, 0.30),
  covariates = subjectData[, c("WT", "eGFR")]
)

fit <- saemix(model, data, control, population = population)
```

The default `covariateMargins = "auto"` selects a continuous family for each
covariate from Normal, lognormal, Gamma, and Weibull using marginal AIC on its
observed values. The selected families are then fixed during the likelihood-
score fit.

For a prespecified or exactly reproducible analysis, provide the margins
explicitly:

```r
population <- gaussianCopulaFrem(
  etaSd = c(0.25, 0.30),
  covariates = subjectData[, c("WT", "eGFR")],
  covariateMargins = list(
    copulaMarginCovariateLognormal(log(70), 0.20),
    copulaMarginCovariateGamma(7.5, 12)
  )
)

```

Continuous, moving-support, binary, ordinal and explicitly ordered categorical
covariate margins are supported under the conditions documented in the
package. The fitted Gaussian-copula FREM can be converted to an exact
conditional full-covariate representation with `copulaFremToFfem()`.

Automatic detection currently applies to continuous covariates only. Binary,
ordinal, and other categorical variables must be declared explicitly because
their scientific meaning cannot be inferred safely from numeric codes. Eta
family selection is also explicit: `etaSd` declares Normal eta margins, while
non-Normal eta margins are supplied through `etaMargins`. This keeps discrete
model selection separate from the fixed-model estimator used in the paper.

## Installation

```r
install.packages("remotes")
remotes::install_github(
  "hiddevandebeek/saemixextension",
  ref = "flexible-margin-frem"
)
```

## Validation

From the repository root:

```powershell
Rscript --vanilla copula/tests/test-minimal-gaussian-copula-frem-api.R
Rscript --vanilla copula/tests/test-score-sa-analytic-gradient.R
Rscript --vanilla copula/tests/test-score-sa-fisher-oracle.R
Rscript --vanilla copula/tests/test-score-sa-end-to-end.R
Rscript --vanilla copula/tests/test-score-sa-gaussian-nested-null.R
Rscript --vanilla copula/tests/test-frem-to-ffem.R
```

## Manuscript and scope

The medium-format manuscript and its figures are under
[`copula/manuscript/medium-form`](copula/manuscript/medium-form). The detailed
implementation boundary is documented in [`copula/README.md`](copula/README.md)
and [`copula/CODE-CLEANUP.md`](copula/CODE-CLEANUP.md).

This branch intentionally excludes the historical general R-vine/common-Q
experiments from the installed package. The supported publication method is
the fixed Gaussian-copula likelihood-score estimator. Convergence is to the
stationary set under the stated assumptions; no global-MLE guarantee is made.
