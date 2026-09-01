# Medium article: theoretical-basis and implementation audit

Date: 2026-08-31

## Scope and decision

The current medium article was checked against Yngman et al. (FREM), Comets et
al. (nonlinear mixed-effects notation), Delyon--Lavielle--Moulines equation
(74) and Fisher identity, Kuhn--Lavielle controlled-MCMC SAEM, and Fort et al.
controlled-Markov stochastic approximation. The article was also traced to the
current `score-sa` implementation and to the saved simulation objects.

Decision: **GO for the declared fixed-Gaussian-copula, prespecified-margin,
stationary-likelihood scope.** This is not a claim of ordinary finite-statistic
SAEM, global maximum-likelihood convergence, automatic family selection, or a
general non-Gaussian R-vine estimator.

## Traceability

| Article statement | Theoretical basis | Implementation evidence | Decision |
|---|---|---|---|
| Standard FREM block model, TPV, B and UPV | Yngman et al., equations 5--12 | Gaussian nested-null tests | Pass |
| Proper joint parameter-covariate distribution | Probability-integral transform and fixed mixed dominating measure | Margin normalization tests | Pass |
| Complete-score mean equals observed-data score | Delyon et al. equation 76 / Fisher identity | `test-score-sa-fisher-oracle.R` | Pass |
| Direct score recursion is not ordinary SAEM | Delyon et al. equation 74 and Section 8.2 | `populationAlgorithm="score-sa"` | Pass |
| Controlled MCMC preserves the conditional target | Kuhn--Lavielle; independence-MH minorization | Exact conditional-prior proposal and MCMC tests | Pass |
| Stationary-set convergence under stated assumptions | Fort et al. H1--H6, Proposition 4.6, Proposition 4.11 and Corollary 4.12 | Gain power 0.80, finite adaptation and runtime diagnostics | Conditional pass |
| Multiple ordered categorical margins | Fixed-support latent-uniform augmentation | Multicategorical Fisher and end-to-end tests | Pass |
| Parameter-dependent support | Fixed-percentile augmentation and pathwise response derivative | Moving-support Fisher test | Pass |
| Gaussian FREM is nested | Normal margins plus full Gaussian vine | `test-score-sa-gaussian-nested-null.R` | Pass |
| Conditional sampling with a covariate subset | Gaussian-score marginal conditioning | High-CRP VPC and individual conditional-density scripts | Pass |
| Translation to a conventional FFEM | FREM conditional Normal block algebra | `copulaFremToFfem`; algebraic and fitted-model tests | Pass when parameter margins are Normal and selected covariates are continuous |
| Translation with non-Normal parameter margins | Latent Gaussian conditional law plus fitted quantile maps | Conditional quantile and exact-simulation tests | Pass as a distributional FFEM |

## Corrections made during this audit

1. Replaced the likelihood values inherited from the older common-Q analysis
   by paired score-SA versus Gaussian-FREM evaluations using 5,000 importance
   samples and identical seeds. Mean differences are -0.01, 5.70, 4.14 and
   12.22 for the four scenarios.
2. Added the generating correlation matrix entries and exact covariate-margin
   parameters to the simulation methods.
3. Qualified categorical support as ordered categorical; permutation-invariant
   nominal models remain outside scope.
4. Disclosed that every estimated structural parameter must belong to the
   population-location block in the current implementation.
5. Corrected the all-Normal parameterization statement: marginal scales and
   correlation determine `Omega_FREM`; covariate locations remain separate.
6. Defined the population parameter subset `theta_p` explicitly and removed an
   undefined response-parameter symbol from the appendix.
7. Removed an artificial 99.95th-percentile cap from the VPC validation
   population. The exact upper-tail comparison gives log-scale RMSE 0.266 for
   Gaussian FREM and 0.123 for flexible-margin FREM.
8. Clarified that the displayed individual distributions condition on fixed
   fitted population parameters.
9. Updated an obsolete test that still expected multiple categorical margins
   to be rejected.
10. Added explicit failure enforcement to the replicated score-study runner
    for post-adaptation projection, backtracking, no-move, numerical-fallback,
    metric-freeze, proposal-freeze and block-schedule criteria.
11. Added an exact FREM-to-FFEM translation API. It returns transformed
    covariates, coefficient and UPV matrices for conventional FFEMs and exact
    conditional quantiles/simulation for distributional FFEMs. It rejects an
    invalid linear translation instead of silently moment-matching.

## Tests rerun

- analytical Gamma, Weibull and lognormal score comparisons;
- Fisher-identity Monte Carlo oracle;
- multicategorical Fisher oracle;
- parameter-dependent-support score test;
- full saemix score-SA end-to-end test;
- Gaussian nested-null test;
- population score backend comparison with common-Q;
- single- and multiple-categorical end-to-end tests; and
- likelihood-score population-algorithm declaration tests.
- FREM-to-FFEM block-algebra, conditional-moment, conditional-quantile,
  categorical, fitted-model and prediction-equivalence tests.

All theorem-relevant tests passed after removal of the obsolete categorical
expectation.

## Remaining empirical limitations

The reported 40-dataset analysis contains ten datasets per scenario and is
correctly labelled exploratory. It does not estimate coverage, cold-start
robustness, family-selection performance, or model-specific verification of
all asymptotic drift and regularity assumptions. A 100-replicate-per-scenario
study is running separately and should replace the exploratory table when
complete.
