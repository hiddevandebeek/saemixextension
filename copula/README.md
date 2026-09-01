# Gaussian-copula FREM for saemix

This research fork adds one supported population estimator to the ordinary
`saemix()` workflow: fixed Gaussian-copula FREM with declared marginal
distributions, estimated by controlled-MCMC likelihood-score stochastic
approximation.

The estimator targets the observed-data likelihood of the declared model. It
is the stochastic-gradient route discussed in Delyon, Lavielle and Moulines
Section 8.2, not finite-statistic SAEM. The formal derivation is in
[`manuscript/score-sa/article.Rmd`](manuscript/score-sa/article.Rmd).

## Minimal workflow

```r
population <- gaussianCopulaFrem(
  etaSd = c(0.25, 0.30),
  covariates = subjectData[, c("WT", "eGFR")]
)

fit <- saemix(model, data, control, population = population)
```

`covariateMargins = "auto"` fits a support-aware starting family from Normal,
lognormal, Gamma, and Weibull. To isolate family selection from estimation,
declare them explicitly:

```r
population <- gaussianCopulaFrem(
  etaMargins = list(
    copulaMarginNormal(0.25),
    copulaMarginNormal(0.30)
  ),
  covariates = subjectData[, c("WT", "eGFR")],
  covariateMargins = list(
    copulaMarginCovariateLognormal(log(70), 0.20),
    copulaMarginCovariateGamma(7.5, 12)
  )
)
```

The initial incumbent uses declared eta margins, normally through `etaSd`.
Categorical covariates must be declared explicitly; non-Normal eta families may
either be declared or selected by the separate posterior screen below.

For optional automatic eta-family screening, first fit the Gaussian incumbent,
then use full posterior draws and start a fresh fixed-model fit:

```r
incumbent <- saemix(model, data, control, population = population)
selection <- copulaSelectEtaMargins(incumbent)
final <- saemix(model, data, control,
  population = selection$population)
```

`selection$table` reports the independent validation-BIC advantage, ESS,
Monte Carlo error, and eligibility of every candidate combination. If
`selection$diagnostics$selectionResolved` is false, increase
`posteriorDraws`; do not treat the current ranking as automatic evidence.

Omitting `population`, or passing `population = NULL`, uses unmodified stock
Gaussian saemix behavior.

## Supported statistical scope

- full fixed Gaussian copula;
- continuous fixed-support margins, continuous moving-support margins with a
  valid CDF/quantile map, and multiple dependent binary, ordinal, or explicitly
  ordered nominal categorical coordinates;
- joint likelihood `p(y, covariates)`;
- missing continuous Gaussian-copula covariates through exact conditional
  augmentation;
- missing categorical covariates and multivariate categorical dependence
  through exact truncated-Gaussian, fixed-support latent-uniform augmentation;
- fixed marginal families and copula structure during estimation;
- one likelihood-score recursion for population locations, marginal
  parameters, Gaussian dependence, and supported residual-error parameters.

Custom margins must declare either `support_fixed = TRUE` or, for a continuous
moving-support law, `support_fixed = FALSE`. The score route rejects undeclared
support. In the latter case, the runtime fixes its percentile coordinate and
propagates the candidate quantile through the response model. The score route
also rejects family/structure selection during fitting, non-Gaussian pair
copulas, conditional-likelihood fitting, indefinite proposal adaptation, and
gain powers at or below 0.75.

## Archived research paths

Historical experiments for general R-vines, direct natural-parameter margins,
adaptive selection, collapse prototypes and retained-particle common-Q fitting
are preserved under `legacy/package-prototypes`. They are not installed with
the package and are not alternative estimators supported by the manuscript.

## Main implementation path

```text
gaussianCopulaFrem()
  -> copulaPopulation()
  -> saemix(..., population = population)
  -> copulaScoreBatchUpdate()
  -> copulaScoreMstep()
  -> copulaGaussianFremPopulationScoreStep()
  -> llisCopula.saemix() / conddistCopula.saemix()
```

Core files:

- `../R/gaussianCopulaFremApi.R` -- minimal public constructor;
- `../R/copulaMargins.R` -- declared marginal distributions;
- `../R/gaussianCopulaFrem.R` -- Gaussian-copula density and conditioning;
- `../R/copulaScoreSa.R` -- current-batch score recursion;
- `../R/copulaLikelihood.R` -- observed likelihood and conditional sampling;
- `../R/main_estep.R`, `../R/main_mstep.R` -- small saemix hooks.

## Verification

```powershell
Rscript copula/tests/test-minimal-gaussian-copula-frem-api.R
Rscript copula/tests/test-score-sa-analytic-gradient.R
Rscript copula/tests/test-score-sa-fisher-oracle.R
Rscript copula/tests/test-score-sa-end-to-end.R
Rscript copula/tests/test-score-sa-gaussian-nested-null.R
Rscript copula/tests/test-gaussian-copula-frem-categorical-end-to-end.R
Rscript copula/tests/test-gaussian-copula-frem-multicategorical-end-to-end.R
Rscript copula/tests/test-score-sa-multicategorical-fisher-oracle.R
Rscript copula/tests/test-moving-support-score-sa.R
```

The systematic theory/code audit is in
[`audits/saem-literature-20260829/SAEM-LITERATURE-AUDIT.md`](audits/saem-literature-20260829/SAEM-LITERATURE-AUDIT.md).
For historical decisions and result provenance, see [`HANDOFF.md`](HANDOFF.md).
