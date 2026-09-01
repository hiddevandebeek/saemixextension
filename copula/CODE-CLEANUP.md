# Publication code boundary

Date: 31 August 2026

## Supported estimator

```text
gaussianCopulaFrem()
  -> copulaPopulation()
  -> saemix(..., population = population)
  -> controlled MCMC E-step
  -> copulaScoreBatchUpdate()
  -> copulaScoreMstep()
  -> copulaGaussianFremPopulationScoreStep()
```

This is the fixed Gaussian-copula FREM model with declared margins and the
likelihood-score stochastic-approximation estimator. Calling `saemix()` without
a population object retains the ordinary Gaussian saemix path.

## Production modules

- `R/gaussianCopulaFremApi.R`: public constructor.
- `R/copulaMargins.R`: marginal distributions and transformations.
- `R/gaussianCopulaCorrelation.R`: correlation/R-vine parameterization.
- `R/gaussianCopulaFrem.R`: density, conditioning, augmentation and sampling.
- `R/copulaScoreSa.R`: complete-data score and stochastic-approximation step.
- `R/copulaScoreRuntime.R`: integration with the saemix M-step.
- `R/copulaEtaSelection.R`: outer full-posterior eta-family screening.
- `R/copulaLikelihood.R`: observed likelihood and individual conditionals.
- `R/copulaFfem.R`: exact conditional full covariate representation.

## Removed from the installed package

The retained-particle common-Q estimator, collapse experiments, direct
natural-parameter margins, adaptive vine learning, Theorem-8 prototypes and
their C++ optimizer are preserved under `copula/legacy/package-prototypes`.
They are excluded from the source package and are not alternative public
estimators.

## Verification contract

Publication changes must pass the minimal API, analytical score, Fisher
identity, end-to-end score, Gaussian nested-null, categorical, moving-support,
custom-margin and FREM-to-FFEM tests. A clean staged source package must install
without compiled code.
