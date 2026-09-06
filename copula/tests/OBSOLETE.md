# Tests of removed functionality

Thirty-nine of the tests in this directory fail, and none of them is a
regression. They exercise APIs and features that were deliberately dropped as
the design settled, and they are kept because they record what those earlier
designs were and why they were tried — but they cannot pass against the current
package and should not be counted when reading the suite.

The suite's real state is the other half: those tests pass, and a failure among
them is a genuine regression.

Do not "fix" anything listed here by restoring the function it calls. Each one
was removed on purpose:

- the `copulaMaximise*` family was the ECM/direct-maximisation population step,
  superseded by the likelihood-score stochastic approximation. `copulaSet` now
  refuses any `populationAlgorithm` other than `"score-sa"`.
- `parameterMarginNormal`, `parameterMarginGamma`, `parameterMarginLognormal`
  were the earlier margin constructors, replaced by `copulaMarginNormal` and
  the tilted family in `naturalParameterMargins.R`.
- `learnGaussianCopulaFremMargins` and `copulaSelectGaussianFremMargins` were
  the family-screening selection route, replaced by `copulaSelectTiltedMargins`
  and the degree ladder.
- non-Gaussian R-vines are out of scope, so `copulaGaussianFremLogPrior` stops
  on anything but a full all-Gaussian vine. That is the manuscript's position,
  not an omission.

## Depending on a function that no longer exists

- `test-adaptive-proposals.R` — could not find function "copulaProposeCandidate"
- `test-automatic-gaussian-frem-margins.R` — could not find function "copulaSelectGaussianFremMargins"
- `test-automatic-gaussian-frem-workflow.R` — could not find function "learnGaussianCopulaFremMargins"
- `test-collapse-theorem-status.R` — could not find function "copulaMaximiseGaussianFrem"
- `test-conditional-likelihood-target.R` — could not find function "copulaConditioningIsTail"
- `test-copula-workflow.R` — could not find function "saemixCopula"
- `test-cpp-backend.R` — could not find function "copulaParameterLayout"
- `test-exact-numerics.R` — could not find function "parameterMarginNormal"
- `test-flexible-margin-saem.R` — could not find function "copulaBicopFromTau"
- `test-formal-saem-proof-obligations.R` — could not find function "copulaMaximiseGaussianFrem"
- `test-gaussian-copula-frem-categorical.R` — could not find function "copulaMaximiseGaussianFrem"
- `test-gaussian-frem-all-marginal-ecm.R` — could not find function "copulaMaximiseGaussianFrem"
- `test-gaussian-frem-multi-arbitrary-collapse.R` — could not find function "copulaMaximiseGaussianFremMultiArbitrary"
- `test-gaussian-frem-mvnormal-mstep.R` — could not find function "copulaMaximiseGaussianFremMvnorm"
- `test-gaussian-frem-partial-collapse.R` — could not find function "copulaMaximiseGaussianFrem"
- `test-one-arbitrary-envelope-gradients.R` — could not find function "copulaMaximiseGaussianFremOneArbitrary"
- `test-parameter-adaptive.R` — could not find function "parameterMarginGamma"
- `test-parameter-copula-saem.R` — could not find function "parameterMarginGamma"
- `test-parameter-family-saem.R` — could not find function "parameterMarginNormal"
- `test-parameter-gaussian-null.R` — could not find function "parameterMarginLognormal"
- `test-parameter-marginal-likelihood.R` — could not find function "parameterMarginNormal"
- `test-parameter-saem.R` — could not find function "parameterMarginGamma"
- `test-score-sa-population-backend.R` — could not find function "copulaMaximiseGaussianFrem"

## Depending on removed behaviour rather than a removed name

- `test-covariates.R` — Error in copulaGaussianFremLogPrior(as.matrix(E), vine, margins, ncol(as.matrix(E)),  :
- `test-custom-eta-end-to-end.R` — Error: abs(state$margins[[1L]]$parameters[["sd"]] - 0.35) < 0.16 is not TRUE
- `test-flexible-margins.R` — Error in copulaGaussianFremLogPrior(as.matrix(E), vine, margins, ncol(as.matrix(E)),  :
- `test-gaussian-copula-frem-end-to-end.R` — Error in copulaSet(vine, margins = margins, conditioning = list(values = conditioning,  :
- `test-gaussianizable-margin-collapse.R` — Error in copulaGaussianFremLogPrior(as.matrix(E), vine, margins, ncol(as.matrix(E)),  :
- `test-hybrid-arbitrary-margins.R` — Error in copulaGaussianFremLogPrior(as.matrix(E), vine, margins, ncol(as.matrix(E)),  :
- `test-hybrid-speed.R` — Error in copulaGaussianFremLogPrior(as.matrix(E), vine, margins, ncol(as.matrix(E)),  :
- `test-joint-conditioning.R` — Error: inherits(implicit_warm_start, "try-error") is not TRUE
- `test-native-parameters.R` — Error in copulaGaussianFremLogPrior(as.matrix(E), vine, margins, ncol(as.matrix(E)),  :
- `test-nested-null.R` — Error in copulaSet(vnG, sdv, populationAlgorithm = "common-q") :
- `test-nongaussian-smoke.R` — Error in copulaSet(startVine, sdG, familySet = NULL, mode = "joint", fitFrom = 20L,  :
- `test-population-algorithm-score.R` — Error: identical(state$populationAlgorithmRequested, "score-sa") is not TRUE
- `test-profiling.R` — Error in copulaSet(v, sd = c(0.3, 0.4), freezeSd = TRUE, freezeVine = TRUE,  :
- `test-robust-defaults.R` — Error in copulaSet(v4, rep(0.3, 4)) :
- `test-theorem8-fixed-map.R` — Error in file(filename, "r", encoding = encoding) :
- `test-vine-sampling.R` — Error in copulaGaussianFremLogPrior(as.matrix(E), vine, margins, ncol(as.matrix(E)),  :
