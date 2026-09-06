# Handoff: Gaussian-copula FREM in SAEM

## Parameter-scale score: psi augmentation (2026-09-05; read first)

The 3000-iteration combined-study rerun (`output/rerun-combined-study` in the
bundle) showed a regression confined to the flexible arm: typical clearance
3.63 against a generating 3.5 (old estimator 3.54), between-dataset SD 0.25
against 0.17, and a smoothing-window score magnitude thirty times the Gaussian
arm's. It was present in the 34 datasets that retained a *lognormal*
clearance, whose model is identical to the Gaussian arm's, so it was not the
Gamma family. A one-dataset experiment fitting that identical lognormal model
both ways (scratch `exp1`) gave averaged score 0.005 and CL 3.45 on the
transformed-additive scale against 0.196 and 3.64 on the parameter scale.

Cause: `copulaScoreBatchUpdate` mapped every natural parameter to a fixed
percentile `u = F(psi; theta)` and the score step then differentiated through
`psi = Q(u; theta)` and the response likelihood (the "moving support"
construction of Appendix A). That augmentation is exact but its per-subject
score variance is set by the individual *response* information rather than
the population information, so the stochastic score is an order of magnitude
noisier and the recursion wanders along the weakly identified location/shape
ridge. For every family in the registry the support does not depend on the
parameters, so the construction was never needed there.

Fix: the parameter scale now holds psi fixed ("psi" route) and uses the
ordinary complete-data score, formed as centered differences of the per-row
complete objective in every coordinate (`copulaGaussianFremNaturalScoreInternal`
in `R/copulaScoreSa.R`; scoreMethod `natural-psi-centered-difference-score`).
The old behaviour is kept behind `copulaPopulation(..., scoreNaturalRoute =
"reference")` for margins whose support moves. Verified: on the natural
lognormal the psi-route score equals the transformed-additive analytic score
to 1e-9 in native coordinates and the Fisher covariance to 5e-5; on the Gamma
margin it equals the global centered difference to 1e-11
(`medium-paper-self-contained/tests/test-natural-psi-score-equivalence.R`).

The 500-dataset study must be rerun with this code before any of its numbers
are published; the bundle currently holds the reference-route rerun with its
regression documented in the Results text.

## Standard errors: terminal averaging pass (2026-09-05)

`copulaPopulation(..., scoreTerminal = K)` freezes the parameter at its
Polyak average for the last K iterations, keeps the chains running and
restarts the per-subject score averages with gain 1/k; the reported Fisher
information is the outer product of those averages at the end of the pass.
Calibrated on 105 datasets (N = 150, six observations, natural lognormal
margins, 2000 iterations with K = 500): SE/empirical-SD ratios 0.99-1.11 on
every parameter, Wald coverage 0.94-0.98. The post-hoc
`copulaScoreFisherInformation(..., thin = 10)` agrees; with unthinned draws
it understates typical-value errors by a third (Monte Carlo variance of the
per-subject averages). `copulaScoreSandwich` (Louis bread) is not calibrated
(clearance location ratio 0.68, sandwich 2.8) and should not be reported.
Default `scoreTerminal = 0` so that the running combined study is unchanged;
the pass must sit inside the smoothing phase (K <= 0.4 total).

## Publication-code cleanup (2026-08-31)

The installed extension now contains one estimator: fixed Gaussian-copula FREM
with declared margins, fitted by likelihood-score stochastic approximation.
General common-Q/R-vine fitting, collapse prototypes, direct natural-parameter
margins, adaptive learning, Theorem-8 experiments and their unused C++ backend
were moved to `copula/legacy/package-prototypes` and are excluded from package
builds. The active Gaussian density/conditioning module was reduced from 2,467
to 773 lines; the runtime was replaced by a compact score-only configuration.
All publication-path regression tests and the original saemix test suite pass,
and a clean staged source package installs without compilation. See
`copula/CODE-CLEANUP.md` and
`copula/audits/PUBLICATION-CODE-AUDIT-20260831.md`.

## Scope pivot (2026-08-25; current priority)

The active research target is now Gaussian-copula FREM versus standard Gaussian
FREM. General R-vine family/structure/truncation learning is deferred. The
Gaussian R-vine is retained only as an exact positive-definite parameterization
of an unrestricted Gaussian correlation matrix. The scientific extension is
flexible continuous margins plus exact collapsed handling of observed and
missing covariates. See `copula/GAUSSIAN-COPULA-FREM.md`.

The first focused implementation slice is complete:

- `R/gaussianCopulaFrem.R` evaluates the exact joint or conditional
  Gaussian-copula FREM density after marginalizing each subject's `NA`
  covariates;
- the ordinary SAEM independence kernel now draws exactly from
  `p(eta | C_obs)` on normal-score coordinates and transports through arbitrary
  continuous eta margins;
- the collapsed common-Q M-step accepts incomplete covariate rows;
- `copula/tests/test-gaussian-conditioning-collapse.R` passes density,
  Schur-moment, arbitrary-margin transport, missing-pattern, M-step, and public
  state gates;
- `copula/tests/test-gaussian-copula-frem-end-to-end.R` passes a nonlinear
  saemix fit with a lognormal covariate and subject-specific missingness.

The analytic Gaussian population update is now implemented. For Normal eta
margins and Normal/lognormal covariate margins, identity/log transformation
turns the entire population block into one weighted MVN. Marginal scales and
all Gaussian R-vine partial correlations are therefore recovered from one
covariance matrix, while missing covariates use pattern-wise Schur EM and the
structured saemix location block uses exact GLS. Unsupported or fixed-margin
models fall back to the literal optimizer, and every collapsed result is gated
against that literal common Q.

Verification and timings (2026-08-25):

- `copula/tests/test-gaussian-frem-mvnormal-mstep.R` matches the complete-data
  weighted MVN MLE to `2e-8` and checks missing-data literal-Q ascent;
- all Gaussian conditioning, end-to-end missing-covariate, conditional-target,
  and Gaussianizable-margin regression tests pass;
- `copula/benchmarkGaussianFremMstep.R` gives 23.42 s for the literal
  finite-difference M-step and 0.062 s for the collapse (378x), with Q differing
  by `3.7e-8`;
- the matched 120-iteration, every-iteration full fit fell from the earlier
  341.49 s to 9.64 s (35.4x); 6.47 s remains in 115 population updates and
  3.77 s in the E-step. See `out/collapsed_full_fit_benchmark.csv`.

Next: use the Gaussian finite-sufficient-statistic form incrementally across
SAEM iterations, then profile exact block/profile accelerators for arbitrary
margins. Do not silently apply the MVN covariance collapse to gamma/Weibull or
other margins whose Gaussianizing transformation depends nonlinearly on the
estimated marginal parameters.

### Automatic margin learning (2026-08-25)

Manual eta-margin declarations are no longer the intended user workflow.
`learnGaussianCopulaFremMargins()` now starts from an ordinary Gaussian saemix
fit, obtains conditional posterior draws, and selects eta and covariate margin
families on one shared weighted common-Q measure with a subject-level BIC
penalty. Default eta candidates are Normal, Student-t, and Laplace; default
continuous covariate candidates are Normal, lognormal, gamma, and Weibull.
This does not fit shrinkage-biased EBEs. Coordinate candidate sweeps jointly
reoptimize every continuous margin parameter and the full Gaussian dependence,
then return a frozen population for a fresh SAEM fit. Manual family lists remain
expert overrides.

Computational dispatch follows the selected statistical model automatically:
Normal/lognormal candidate vectors use
`gaussian-copula-frem-mvnormal-em`; other supported families use the literal
common-Q optimizer. `test-automatic-gaussian-frem-margins.R` checks family
recovery and backend dispatch, and
`test-automatic-gaussian-frem-workflow.R` checks the complete Gaussian-fit to
fixed-population learning path. Both pass. Independent-posterior and final
observed-likelihood validation are still required around the selected winner
before this becomes a one-call production default.

**Goal:** a robust, likelihood-correct Gaussian-copula FREM implementation with
flexible continuous margins and exact collapsed covariate handling.

**Status:** the fixed-model estimator now optimizes the correct common empirical
auxiliary surrogate over saemix population locations, all free parameters of
declared admissible continuous margins, and all native pair-copula parameters.
Valid block updates operate against that same common objective. This is the correct
fixed-specification observed-likelihood target; a finite diagnosed endpoint is
not automatically an MLE. Likelihood stationarity oracles and the
main copula-aware post-processing paths are implemented. Broader repeated
validation, accurate observed-likelihood Monte Carlo diagnostics, and
copula-aware uncertainty quantification are still required before the method
should be advertised unconditionally as production MLE software.

### Current likelihood-target and M-step continuation (2026-08-25)

- Covariate-augmented vines now expose `likelihoodTarget = c("joint",
  "conditional")`. The default preserves the existing proper joint target
  `g(y,c)`. The conditional route optimizes the distinct exact objective
  `log p(eta,c) - log p(c)` and targets `g(y|c)`.
- Exact conditional fitting currently requires the observed conditioning
  coordinates to form the terminal self-contained block of the R-vine order.
  `copulaConditioningSubvine()` then evaluates `p(c)` as a lower-dimensional
  vine. Incompatible orders and unconstrained structure reselection fail before
  fitting; they are not silently treated as conditional models.
- `copula/tests/test-conditional-likelihood-target.R` verifies sub-vine
  extraction, numerical normalization of `p(eta|c)`, common-Q ascent, API
  acceptance, and fail-closed ordering. Existing joint-conditioning,
  arbitrary-margin hybrid, and covariate tests still pass.
- The public defensive-importance likelihood now handles covariate-augmented
  snapshots: joint fits report `g(y,c)` and conditional fits report `g(y|c)`.
  MAP and posterior random-walk paths augment every latent proposal with the
  subject's fixed conditioning row. Conditional-prior independence proposals
  remain disabled until an exact conditional sampler is exposed there.
- The 400+100 one-cycle hybrid block-GEM phase took 418.0 seconds versus the
  724.1-second reference, but its paired conditional observed log likelihood
  was 0.153 lower (MCSE 0.019). An 800+200 continuation took 254.6 seconds only
  because its final simultaneous optimizer fell back (`conv=-1`, infinite
  gradient); its likelihood was 1.890 lower (MCSE 0.041). The current block
  schedule is rejected.

### Active publication run (2026-08-23 05:xx CEST; newest status)

- The long-chain calibration is complete and its authoritative direct checker
  passed all three N=100 Gumbel fixtures. The accepted calibration freeze is
  `copula/study_publication/calibration_v3/freeze.rds`; its freeze ID is
  `08572aa73140ae631b22e8135b86016b8ed21a0f2b586498dc9adf7be543573f`.
- The frozen balanced primary study is running with 10 workers into
  `copula/study_publication/out`. Its manifest contains all 660 prespecified
  cells (110 replicates for each Gaussian/Gumbel by N=50/100/300 stratum), with
  300+200 reported-fit iterations and 2,000 importance samples. Resume or audit
  progress by counting `out/primary/**/cell_r*.rds`; never patch a file listed
  by `study_code_hashes()` until summarization, figures, and generated Results
  have been completed under this freeze.
- On this Windows host, the `Start-Process` wrappers emit blank/nonzero exit
  codes even when `run_cell.R` writes and validates its atomic RDS. The balanced
  launcher therefore prints `FAIL` and continues. Treat the atomic files plus
  `study_cell_issues()` and the final source-hashed summarizer as authoritative,
  not the wrapper label. A direct validation of the first completed cells found
  zero schema/freeze/provenance issues.
- The user-owned `bigstudy.R` process (formerly PID 23000) has exited normally.
  The 10-worker primary run already saturates the available logical processors;
  do not add workers opportunistically even though physical-memory headroom has
  increased.
- Independent theory/literature/audience rounds 53--62 are in
  `copula/audits/`. They again passed the fixed-model observed-likelihood, EM/KL,
  Fisher, common-Q, Gaussian-nesting, MH, and convergence-boundary mathematics;
  added Duval and Robert-Grani\'e's 2007 heuristic non-curved-exponential
  SAEM-MCMC predecessor; and reorganized the paper for a first-time reader.
  The main Methods now carry the estimator and scientific design, while
  implementation/freeze edge cases are in Appendix B. The last remaining
  publication blockers are the validated primary results, the quantitative
  abstract/Discussion/Conclusion rewrite, rendered article/supplement QA, and
  author-supplied metadata/archive DOI.
- Round 63 found public-help inconsistencies that must be repaired only after
  the evidence bundle releases the source freeze: scope the inherited package-
  wide MLE/standard-error wording to ordinary Gaussian `saemix`; correct the
  `saemixCopula()` `criterion` help to describe the common-Q penalty rather than
  proposal AIC; document high-level versus legacy `familySet` behavior, both
  margin declaration routes, vector `truncLvl`, and warm-started fresh state;
  regenerate Rd files; and audit the broad namespace exports. The historical
  HANDOFF and README portions that were safe to change during the run have
  already been labelled/repaired. See
  `copula/audits/round63-public-doc-consistency.md`.
- Rounds 64--68 reconciled the four theory companions with the normative
  manuscript. They now distinguish the current continuous routes from
  prospective discrete/mixed or automatic margin-family learning, state the
  fixed parameter-independent transformation contract, and avoid conflating
  common-Q stationarity with MLE attainment. The final companion closure has no
  substantive open finding.
- A standalone base-R independent primary reanalysis was prepared in
  `copula/audits/independent_primary_reanalysis.R`. It does not source the
  official reducer, defaults to a non-writing dry run, and is undergoing a
  second independent static audit before it may be executed after all 660
  atomic cells are complete.
- Frozen `PROTOCOL.md` prose says “15 mapped result objects,” whereas the
  source-hashed `study_mapped_result_inventory()` contains 17 (`final_retention`
  and `runtime_environment` are separate objects). The executable inventory is
  authoritative. Do not edit the frozen protocol during the run; record this as
  a documentation erratum after validated Results and PDFs are finalized.

### Publication continuation (2026-08-23, supersedes the 2026-08-22 continuation)

- Runtime verification is now complete. Exact-numerics E1--E17, core margin
  and likelihood tests, hybrid/direct and adaptive-proposal tests, the
  three-dataset Gaussian nested-null comparison, non-Gaussian smoke tests,
  direct-parameter tests, and the source-to-paper mapping tests all pass.
  Round 51 independently re-audited the last repairs and returned FREEZE PASS.
- The last pre-freeze repairs make failed-but-timed oracle runtime denominators
  explicit, prevent the hybrid optimizer from exposing an internal sentinel as
  an exact common-Q endpoint, restrict sparse-refit beta holding to active
  population-location coefficients, and verify native reconstruction of a
  genuine four-dimensional truncated vine.

- The authoritative article source is `copula/manuscript/preprint.Rmd`; the
  companion source is `copula/manuscript/supplement.Rmd`. The title and abstract
  now state the exact claim boundary: the fixed-model observed-data likelihood
  and EM identity are established; classical MCMC-SAEM is used for finite-
  statistic models and finite-dimensional MCMC score stochastic approximation
  is the theorem-backed general route. Neither route guarantees the global MLE.
- Independent rounds 31--36 rechecked every heterogeneous-margin density,
  change-of-variables, EM/KL, Fisher, coupled-score, common-Q, Gaussian-nesting,
  and convergence-boundary claim. They repaired probit Jacobians, made support-
  transform saturation fail closed, replaced a too-strong substitution argument
  by the PIT push-forward proof, and fixed a malformed Proposition 1 math span.
  Round 36's final mathematical verdict is PASS with no open finding.
- Round 37 independently reread the original SAEM/SAEM-MCMC literature. It
  corrected one attribution: Delyon et al.'s general equation (6) already is a
  whole-function recursion equivalent to a weighted empirical measure; their
  *convergence proof* is what specializes to a curved-exponential finite-
  sufficient-statistic representation. Its independent post-repair verdict is
  PASS with no open P1/P2 finding.
- Rounds 33 and 35 closed the prospective publication pipeline. Atomic worker
  schema 7 and every freeze/manifest now bind `Rscript --vanilla`, the fixed RNG
  algorithms, R/platform and native dependency versions, and hashes of loaded
  `saemix` and `rvinecopulib` DLLs. The Results inventory has 17 mapped objects; the
  seven-output evidence manifest covers the article fragment, Supplementary
  Tables S1--S4, result RDS, and mapped-result-object inventory (not rendered
  display-cell lineage). Full-denominator family,
  exactness, likelihood, optimizer, and stage-specific oracle-runtime reporting
  passed a final static audit.
- Round 38's first-reader audit found a rigorous core but correctly classified
  the current draft as protocol-facing until frozen quantitative results exist.
  Arms, tail estimand, scope, declaration routes, and analysis sets now precede
  protocol detail; power and computational-contract algebra are in Appendix B.
  The final abstract, authored Results narrative, compact main display set,
  Discussion, and Conclusion must be written from the frozen outputs.
- The unrelated `ticbig.R` jobs have exited. The two user-owned `bigstudy.R`
  jobs remain protected and must not be stopped. Available memory was sufficient
  for the completed regression suite. The next computational step is the fresh
  three-cell 500+500 calibration under the final source hash, followed—only if
  its gate passes—by the balanced 660-cell 300+200 primary study.

### Publication continuation (2026-08-22, supersedes older study/status text below)

- The current paper source is `copula/manuscript/preprint.Rmd`; the earlier
  17-page PDF is an interim artifact and is not submission-ready. The draft now
  proves the fixed-specification observed-likelihood/EM target, gives the
  coupled margin--copula score and incumbent-relative cross-model identity,
  separates finite learning from fresh frozen-vine estimation, and moves the
  detailed convergence/numerical audit into appendices.
- Independent audit rounds 1--18 are indexed in `copula/audits/INDEX.md`.
  Rounds 9--18 again passed normalization, Jacobians, EM/KL, Fisher, common-Q,
  scaling, conditional stationarity, Gaussian nesting, and cross-model algebra.
  The primary-literature audit was performed afresh online and the manuscript
  now includes the closest asymmetric-NLME and pre-2017 copula mixed-model
  precedents while distinguishing pseudo-likelihood from observed-data
  likelihood. The most recent audits additionally repaired SAEM naming,
  claim-specific regularity, fixed-vine admissibility, change-of-variables
  conditions, posterior absolute continuity, bivariate evidence scope, and
  first-reader presentation. No mathematical P1 blocker remains.
- Round 15 found no published convergence theorem covering the generic
  function-valued empirical-Q implementation. That historical backend is no
  longer the proposed general estimator: Section 20 below records the finite-
  dimensional MCMC score route that closes the theory gap. The retained-pool
  code remains a development reference only.
- Round 16 stress-tested every displayed likelihood and convergence-roadmap
  identity after that expansion. Continuity at moving gradient limits, entropy
  finiteness, positive-definite Gaussian nesting, and the one-draw Delyon
  specialization were made explicit; closure found no remaining P1 or P2
  mathematical-validity issue.
- Round 17 found and incorporated the closest 2025 joint vine-margin
  observed-likelihood/selection precedent, plus the historical PIT-based NLME
  and longitudinal copula-random-effects precedents. Novelty is now explicitly
  limited to the incumbent-conditional mechanistic-NLME SAEM construction and
  separated screened-learning/fixed-refit workflow; `SAEM-MCMC` terminology is
  standardized. Closure found no remaining P1 or P2 literature-positioning
  issue.
- Round 18 repaired the finite-pool option so it retains complete newest batch
  identities rather than weight-sorted rows vulnerable to floating-point tie
  drift. It also made restart-burn boundaries and pool metadata/multiplicity
  contracts fail-closed and added targeted tests. The uncompressed publication
  settings are unchanged; closure found no remaining P1 or P2 alignment issue.
- Round 21 independently rederived every central likelihood and score identity
  and found no P1 defect. Its sole P2 was resolved by stating the topology,
  derivative-envelope/interchange, and differentiable maximizer-branch
  premises of the optional local measure-valued score-flow calculation. The
  calculation remains explicitly conditional, not a software convergence
  theorem.
- Round 20 independently aligned notation with the historic SAEM and `saemix`
  sources. Coordinate-specific covariate anchors, random/realized PITs, and the
  preceding-iterate Simulation--Stochastic approximation--Maximization cycle
  are now explicit; posterior-measure, EBE, Di{\ss}mann, entropy, and iterative-
  M terminology were also corrected. Its closure found no P1 or P2 issue.
- Round 24 independently reread the original SAEM and SAEM-MCMC papers online,
  together with the modern curved-exponential and particle-SAEM caveats. It
  again passed the likelihood theory and convergence boundary with no P1/P2.
  The manuscript now describes Delyon's LOC conditions by their actual
  isolation/definiteness/noise roles, uses Kuhn--Lavielle's uniform-ergodicity
  language, states the mixing and transient effects of finite MCMC composition,
  and says the fixed model defines the intended observed-likelihood objective
  rather than using a potentially ambiguous bare convergence ``target'' claim.
- The final round-21 closure restricted the weighted norm to the finite-$w$-
  moment signed-measure space and required both measures in the formal flow to
  have finite local $w$-moments. The referee then found no remaining P1 or P2
  mathematical issue.
- Round 19 found major simulation-estimand and reporting defects before the
  primary freeze. Round 22 repaired them prospectively, and round 23 strengthened
  the design to 110 replicates per
  cell, an explicit tail-effect planning calculation, a deployable qualifying-
  learned-or-Gaussian-fallback primary endpoint, secondary fresh-endpoint
  contrasts, stricter 10%-of-truth operational margins, fail-closed global
  conjunctions, full arm performance, tail robustness, direct fallback/failure
  provenance, and paired runtime ratios. See
  `copula/audits/round22-simulation-repair.md`.
- The round-10 study audit found and repaired two pre-freeze blockers:
  inference from incomplete output strata and loss of the completed learner
  record when a fresh refit throws. The subsequent from-scratch audit drove
  schema 6: it preserves learning/refit stages, reconstructs deterministic
  seeds and all frozen provenance, recomputes endpoint and likelihood status
  rather than trusting `valid`, verifies exact scenario tail truth, requires
  at least 99 of 110 qualifying Gaussian-reference pairs for the deployable
  composite workflow, and explicitly labels primary,
  oracle-secondary, and descriptive analyses. Freeze identity now binds the
  canonical plan and output root. Three post-repair static passes found no
  remaining ranked pipeline issue; runtime tests remain required.
- The publication study is `copula/study_publication/PROTOCOL.md`: 660 datasets
  across Gaussian/Gumbel truth and N=50/100/300, 110 replicates per cell,
  Gaussian and learned arms everywhere, plus a fixed-Gumbel oracle under
  Gumbel truth. Primary chains are 300+200; a three-dataset 500+500 calibration
  must pass under the identical source hash first.
- Round 23 independently reconstructed the deployable estimand and every
  denominator. Its prospective repairs distinguish per-size from joint-Holm
  power (idealized all-three power about 0.806 at 110 replicates), record and
  validate each worker execution environment in the atomic schema-6 cell,
  atomically bind derived summaries to ordered input/output hashes, and make
  figure generation consume the validated summary RDS rather than a mutable CSV.
  The source-hashed `test_summary_logic.R` exercises fallback, unavailable,
  99-pair, zero-variance, and malformed-conjunction branches. Static closure is
  in `copula/audits/round23-pipeline-referee.md`; R execution remains pending the
  memory guard.
- `R/parameterMargins.R` now separates literal exact direct-Q evaluation from
  the optimizer's finite internal penalty, rejects invalid starts/endpoints,
  exactly rebuilds the returned model, and reports a bound-projected gradient.
  Exact-numerics cases E13/E13b cover invalid and adjacent valid Normal tails.
- The current calibration freeze is intentionally absent. The preceding empty
  freeze was archived as
  `calibration_v3_stale_directsentinel_20260822_1900` after the source repair.
  Do not reuse any stale calibration directory.
- At this update, two unrelated `ticbig.R` processes (PIDs 30016 and 31568)
  consume almost all committed memory. Do not stop them. Run no R until they
  exit or committed-memory headroom is sustainably at least 3--5 GB. Then run
  `copula/tests/test-exact-numerics.R`, the relevant core suites and full nested
  Gaussian null, create a new calibration freeze, pass its gate, and launch the
  balanced primary manifest. Results placeholders may be replaced only from
  the frozen primary summary after independent reproduction.

### Direct natural-parameter continuation (2026-08-21, current)

The scientific target is now stated precisely as the **observed-data marginal
likelihood** obtained by integrating latent natural individual parameters
`psi`. It is not a “marginal distribution” target and it is never a likelihood
of point EBEs.

- `R/parameterMargins.R` implements support-aware Normal, Student, Laplace,
  lognormal, Gamma, Weibull, and beta margins for `psi`, each with an explicit
  anchor and link. `p_phi = p_psi(T(phi)) |det DT(phi)|` is used in every MCMC
  and likelihood path.
- `saemix(..., population = copulaPopulation(...))` is the package-native API.
  Fitted population state survives serialization. Simulation, MAP,
  conditional-distribution sampling, prediction/newdata, and defensive-IS
  likelihood use that state without mutable global configuration.
- Fixed-model Maximization jointly optimizes anchor/covariate coefficients,
  all free native margin parameters, and all free native pair-copula
  parameters on the same weighted absolute-phi particle approximation to the
  EM auxiliary function. `omega` is deliberately `NA` for direct populations;
  proposal covariance is stored separately and conditional population
  covariance is requested explicitly.
- Finite adaptation can now propose complete compatible margin-family vectors
  together with vine structure, family, rotation, and truncation. Proposals are
  reoptimized and compared on one common particle measure. Adaptation must end
  in a fresh, full-length fixed-fingerprint SAEM-MCMC run.
- Likelihood oracles pass for every built-in continuous margin, all internal
  transforms, population anchors/shapes, native copula parameters, and a
  residual parameter. The full-data likelihood engine agrees to
  `1.4e-14` under equivalent direct-lognormal and transformed-additive Gaussian
  declarations. End-to-end SAEM-MCMC recovery passes for all seven built-ins.
- The exact Gaussian-null end-to-end check matches stock fixed effects within
  `0.1%` and the working-scale covariance within `1.5%` on its prespecified
  metric. Heterogeneous Gamma/lognormal--Gumbel recovery, persisted likelihood,
  simulation, MAP, and conditional posterior paths pass.
- `copula/study_direct/` replaces the paused eta-scale study. Primary scientific
  arms use matched 500-iteration `c(300, 200)` schedules, with matched
  1,000-iteration `c(500, 500)` sensitivity runs, direct population activation
  from iteration 1, and a population Maximization every iteration. The older
  400+100 delayed-activation and sparse-update recipes are speed experiments,
  not correctness defaults.
- `copula/manuscript/preprint.pdf` is a rendered 17-page source-aligned draft.
  It distinguishes the exact likelihood target from stochastic-approximation
  convergence, a local stationary point, and a global maximum.

Remaining production gaps: likelihood ESS/Monte Carlo error is inadequate in
some quick nonlinear cells and must be improved before the large study;
copula-aware standard errors/FIM remain unimplemented; discrete direct margins
still need a label-aware latent state; and the general MCMC score backend must
be integrated and qualified. The function-valued particle approximation is not
the theorem-backed production route.

`llis.saemix()` now has a copula-aware defensive-importance-sampling path. Each
returned fit stores a compact immutable prior snapshot, so scoring is independent
of mutable global `.cop` state and survives serialization. Stock MAP, FIM,
Gaussian quadrature, and conditional-distribution post-processing are skipped or
refused for copula fits until their non-Gaussian versions exist.

### Continuation update (2026-08-21)

- The joint particle M-step now retains **absolute latent phi plus subject design
  rows** and jointly updates every estimated eta-location coefficient, including
  covariate effects. `tests/test-covariates.R` pins both empirical-Q recovery and
  an end-to-end Gaussian nested null with a covariate.
- Robustness is now split into an **adaptive learning phase after the initial
  Gaussian fit**, followed by a fresh fixed-model likelihood-targeting phase. Structure, families,
  rotations and truncation may jump while generating proposals, but the final
  run carries and asserts an invariant discrete-model fingerprint.
- Do **not** silently truncate or switch to Gaussian inside a production run.
  Truncation is a different model and the default full Gaussian vine must retain
  the arbitrary-MVN stock null. The default collapse guard records and warns
  without mutation; an outer controller may reject the candidate and return the
  separately fitted Gaussian baseline.
- `copula/R/adaptiveVine.R` implements the adaptive layer. Standard
  `rvinecopulib` parametric selection proposes structure, family, rotation and
  truncation candidates from subject-level conditional-draw slices. Complete
  marginal-family vectors can jump through the same interface. Every candidate
  is reoptimized on one common absolute-phi particle Q; failed candidates cannot
  win and the incumbent wins ties or gains below `minImprovement`. Slice-level
  BIC is proposal ranking only. Adaptation must stop after a finite phase, then
  restart a fresh fixed-fingerprint SAEM fit and compare observed likelihood.

---

## Historical development notes below; not current API or claim wording

The sections below preserve implementation history and old debugging lessons.
For current interfaces and scientific claims, use `copula/README.md`, the
generated help for `saemixCopula()` and `copulaPopulation()`, and
`copula/manuscript/preprint.Rmd`. Any older statement that a finite endpoint
``is the MLE'' or that rejected `pool`/`sa` modes are runnable is superseded.

## 0. Read this first

Three things will cost you a day each if you skip them.

1. **The ordering contract** (§3). A mismatch between the vine's variable order
   and the eta column order is silent -- you get a valid joint law, just not the
   one you specified. `tests/test-ordering.R` pins it exactly.
2. **Never call `set.seed()` inside the SAEM loop.** A seeded common-random-number
   block in the M-step reset R's global stream, so the E-step replayed identical
   MCMC draws and acceptance uniforms every iteration and the fit ran away
   monotonically. Use `withSeed()` (saves/restores `.Random.seed`). Test N7.
3. **Keep discrete selection OUT of the M-step.** Structure, families and
   truncation level must be chosen once and frozen (or in an outer loop).
   Discrete selection inside the M-step breaks Delyon et al. (1999) (M1) (Theta
   is no longer a fixed open subset of R^p), (M5) and (SAEM2): theta-hat jumps at
   selection boundaries so the mean field is discontinuous.

---

## 1. Environment

- Repo: `C:\package\saemix-copula`, branch **`copula-eta`**. Clone of
  `saemixdevelopment/saemixextension` (saemix 3.5). **Local only, never pushed.**
- R is **not on PATH**: `export PATH="/c/Program Files/R/R-4.5.3/bin/x64:$PATH"`
- Load with `devtools::load_all("C:/package/saemix-copula")`.
- In R use Windows paths (`C:/package/...`), not `/c/package/...`.
- Key packages: `rvinecopulib` 0.7.3.1.0 (canonical), `VineCopula` 2.6.1 (oracle
  only), `mvtnorm`, `numDeriv`, `combinat`.
- The box is CPU-contended (another Claude session runs jobs here). Keep runs
  small; check `Get-CimInstance Win32_Process -Filter "Name='Rscript.exe'"`
  before assuming a job is yours.

**Separate smoke tests from scientific runs.** N=40-60, short iteration schedules,
M=500 IS draws, and 1-2 replicates are useful only for software debugging. Any
convergence or method comparison must use at least the matched 500-iteration
schedule above, with the prespecified 1,000-iteration sensitivity and enough
replicates for Monte Carlo precision.

---

## 2. The model

```
eta_j ~ F_j(.; alpha_j)                   arbitrary centred proper margins
(u_1..u_d) ~ R-vine copula C(.; gamma)    flexible dependence
u_j = F_j(eta_j; alpha_j)
```

`transform.par` maps individual model parameters; it does not make the eta law
on that transformed scale arbitrarily non-Gaussian.  The current framework
therefore models both margins and dependence.  Built-ins cover Normal,
finite-variance Student, and Laplace margins; `copulaMargin()` and
`copulaMarginDistribution()` are the general finite-dimensional extension
boundary.  Population location remains `X beta`, so eta margins must be centred
and may not duplicate a free location parameter.

For continuous coordinates the prior log density is

```text
sum_j log f_j(eta_j; alpha_j) + log c_V(F_1(eta_1),...,F_d(eta_d); gamma).
```

Mixed continuous/discrete likelihood dispatch using CDF left limits is exact at
the distribution layer.  Discrete *latent eta* is not yet enabled end to end:
it needs a label-aware MCMC state and particle pool.  See `MARGINS-DESIGN.md`.

**Nested null.** An R-vine with Gaussian pair copulas at every edge *is*
MVN(0, Omega) exactly (Bedford-Cooke, partial correlations). So stock saemix is
the nested null, and any Gaussian-vine run must reproduce it. This is the single
most useful test you have -- use it after every change.

---

## 3. Ordering contract (load-bearing)

**rvinecopulib is canonical.** `dvine_structure(1:d)`;
`pair_copulas[[t]][[e]]` is the copula of `(e, e+t)` given the variables
strictly between them; and

```
vine variable j  ==  eta column j  ==  varList$ind.eta[j]
```

`VineCopula` uses a **transposed** array convention. It is an independent oracle
only -- never mix its objects with rvinecopulib's.

Two traps `tests/test-ordering.R` pins down exactly (density-based, not Monte
Carlo, so they pass to ~1e-14 or fail loudly):

- The inverse-Rosenblatt map `eta = T(z)` is **order-dependent even for a plain
  MVN**, so any non-centered parametrisation must pin the eta order.
- A vine containing a **rotated** pair copula does **not** survive relabelling
  (rotated copulas are not exchangeable): log-density shifts by ~1.1 versus
  ~1e-15 for Gaussian or unrotated Gumbel. With rotations in play the eta order
  is part of the model, not a free choice.

---

## 4. What is implemented

Package changes (only these four files differ from stock saemix):

| file | change |
|---|---|
| `R/copulaPrior.R` | new -- the whole copula prior, M-step, helpers |
| `R/main_estep.R` | null-gated copula path in the E-step |
| `R/main_mstep.R` | null-gated copula path + location-shift application |
| `DESCRIPTION` | `rvinecopulib` import, Collate entry |

The ordinary Gaussian path is retained when no flexible population is active,
subject to the separately documented fork-wide probit and extreme-tail
corrections. That gating makes the nested-null test meaningful: the same
`estep`/`mstep` code runs under two priors that are provably the same
distribution.

### Historical low-level API snapshot

This snapshot records an earlier experimental interface. Current recommended
entry points are `saemixCopula()` for finite dependence learning followed by a
new stochastic recursion, and `saemix(..., population=copulaPopulation(...))`
for a caller-declared fixed model. The current implementation rejects the
approximate `pool` and `sa` fitting modes shown below.

```r
copulaSet(vine, sd,
          familySet  = NULL,        # families the M-step may select from
          mode       = "joint",     # "joint" | "pool" | "sa"
          truncLvl   = Inf,         # 1L = truncate after tree 1
          refitEvery = 1L, fitFrom = 1L,
          poolMax    = Inf,
          freezeSd   = FALSE, freezeVine = FALSE,
          sdFromSS   = TRUE, verbose = FALSE)
copulaGet()      # snapshot: $vine, $sd, $delta, $famFixed, $trace, ...
copulaClear()    # revert to stock saemix
```

### Historical M-step experiments

- **`joint`** -- optimized the fixed model's common empirical auxiliary
  surrogate over location shifts, marginal scales, and native pair-copula
  parameters. This is not by itself an observed-score, convergence, or MLE
  certificate.
- **`pool`** -- rejected legacy experiment: IFM used the second moment for SD,
  then fit the copula on
  `pnorm(eta/sd)`. Q-hat carried as a weighted particle pool (Delyon eq. 6 with a
  weighted empirical measure). Consistent, not efficient.
- **`sa`** -- rejected legacy experiment: fit the vine to the current draws,
  then a Robbins-Monro
  step on tau. Cheapest. Carries a Jensen shift O(1/n_draws) -- RM on the
  parameter converges to `E[tau-hat]`, not the maximiser of the expected
  pseudo-loglik, and iteration-averaging does **not** remove it.

### How the E-step generalises

- **Kernel 1 is free.** It proposes *from the prior*, so the prior cancels in the
  MH ratio whatever the prior is -- only the proposal generator changes
  (`copulaRandEta`).
- **Kernels 2/3** need only the prior density (`copulaUeta`).
- **Kernel 4** (Laplace/MAP, the f-SAEM kernel) needs the copula Hessian and is
  **refused** rather than run with a mismatched Gaussian prior. Set
  `nbiter.mcmc[4] = 0`. Enabling it is on the TODO list.

### Simulated annealing is load-bearing

`copulaMstep` applies saemix's variance floor (`max(new, alpha1.sa * prev)`) to
the marginal SDs, anchored to the **initial** sd. Do not remove it. Without it a
weakly identified eta collapses: a tighter prior shrinks its posterior draws,
which shrinks sd, which tightens the prior -- and the copula absorbs the lost
marginal scale as near-perfect dependence (observed: tau -> 0.84, implied
rho 0.91 against a truth of 0.4). The vine parametrisation is *structurally* more
exposed to this than a covariance, because a copula fitted on normal scores is
scale-invariant, whereas a sample covariance cannot exceed sd_i * sd_j.

---

## 5. What is established

### Theory

- SAEM's exact-marginal-likelihood target **never depended on the prior being
  Gaussian**. (M1) is what makes the algorithm finite-dimensional and the proof
  work, not what makes the objective right.
- **The IFM M-step is consistent but not an MLE.** With
  `log p = sum_j [-log sd_j - eta_j^2/(2 sd_j^2)] + log c(u)`:
  ```
  d/d sd_j log p = [ -1/sd_j + eta_j^2/sd_j^3 ] + D_j,  D_j = -(1/sd_j) z_j phi(z_j) d_j log c
                   \____ what IFM solves ____/          \____ what it DROPS ____/
  ```
  `E[D_j] = 0` at the truth: with `g(u) = z phi(z)`, `g'(u) = 1 - z^2`, so
  integrating by parts `E_c[g d_j log c] = -E_c[1 - z^2] = 0` since `u ~ U(0,1)`
  implies `z ~ N(0,1)`. Verified at n=4e5 across Gaussian/Gumbel/Clayton, d=2
  and d=4 (all |values| < 0.006); the counterfactual fails as predicted (a wrong
  sd in the PIT moves it ~6 SE off zero). So IFM solves an unbiased *estimating
  equation*, not the score -- consistent, inefficient, no EM ascent guarantee,
  and invalid observed-information SEs/LRTs. This is an ES/IFM estimator
  (Elashoff & Ryan 2004, JCGS 13(1); Joe & Xu 1996).
- **`mode="joint"` repairs the common-surrogate score defect** -- see §5 numbers below.
- **1-truncation is established practice**, not a novelty. Nikoloulopoulos's vine
  copula mixed models (the one published line putting a vine on *latent* random
  effects) fix 1-truncation a priori from d>=4, at N=10-40. Licence: Joe, Li &
  Nikoloulopoulos (2010, JMVA) -- a vine has (tail) dependence in **all**
  bivariate margins provided the **tree-1** pair copulas do.
- **Parsimony argument.** A full vine at d=4 has 6 dependence parameters --
  *exactly* as many as an unstructured MVN correlation matrix. Shape flexibility,
  zero parsimony. 1-truncation halves it.

### Empirical (all numbers reproducible from `out/*.rds`)

**Joint common-Q update reduces the surrogate score residual** (d=2, N=40,
2 reps, `out/llcheck_{ifm,joint}_N40.rds`):

|          | IFM \|t\| | IFM shift | JOINT \|t\| | JOINT shift |
|----------|-----------|-----------|-------------|-------------|
| mu1      | 44.27     | 0.0058    | 4.39        | 0.0004      |
| mu2      | 32.67     | 0.0056    | 1.83        | 0.0002      |
| logSd1   | 17.03     | 0.0110    | 1.61        | 0.0010      |
| logSd2   | 22.86     | 0.0190    | 0.93        | 0.0007      |
| atanhTau | 1.67      | 0.0008    | 1.64        | 0.0008      |

Read the **shift**, not |t| -- a small gradient SE inflates |t| without meaning.
Residuals are 0.02-0.16% on the log scale vs 0.6-1.9% for IFM. Signs flip across
replicates in both arms, i.e. the residual is Monte Carlo noise, not bias.
`tau` was already approximately stationary for that empirical surrogate under
IFM because the copula block solved its own score equation. This is not an
observed-score or MLE certificate.

**Truncation fixes the failure mode and beats stock** (d=4, Gaussian truth,
5 reps, `out/trunc.rds`) -- mean |relative error| vs truth:

| arm | overall | sd(V2) |
|---|---|---|
| stock | 0.0437 | 0.114 |
| full vine | 0.0708 | 0.215 |
| **tree-1 truncated** | **0.0391** | **0.078** |
| frozen | 0.0455 | 0.124 |

**It breaks only when BOTH a weakly identified eta AND a full vine are present.**
Same d=4 model, only the sampling design changed (`out/decisive.rds`):

| | V2 shrinkage | copula/stock ratio |
|---|---|---|
| `iv2` | 0.213 | 1.17, and **intermittent** (per-rep 1.42, 0.90, 0.99, 1.57, 1.04) |
| `iv2rich` | 0.028 | 1.00, identical parameter by parameter |

Intermittency matters practically: you fit once and cannot tell which you got.
That is an argument for truncating by default, not for treating the full vine as
usually-fine.

**The threshold is free pair-copula parameters, not dimension** (`out/ladder.rds`):
d=2 (1 edge) and d=3 (3 edges) match stock replicate by replicate under a
Gaussian truth; d=4 full (6 edges) degrades; d=4 truncated (3 edges) behaves
like d=3.

**It estimates the TAIL better, which is the point** (d=4, Gumbel truth, 6 reps,
`out/rec2_sa.rds`, `out/recipe.rds`). Kendall's tau is the *wrong* metric -- the
twins are matched on tau by construction, so both models target the same tau and
the Gaussian estimates it more efficiently.

```
dLL = 2*(llCopula - llGauss), exact marginal likelihood
  Gaussian truth  mean  +1.40  range [-2.05, +3.80]
  Gumbel truth    mean +15.42  range [-0.02, +32.64]

joint tail probabilities, Gumbel truth
  pairHigh  true .4900   gauss .3469 (-29%)   copula .5020 (+2%)
  allHigh   true .0114   gauss .0036 (-3.2x)  copula .0097
  allLow    true .0017   gauss .0038 (+2.2x)  copula .0015
per-replicate |pairHigh error|: gauss .143 vs copula .050 -- copula better 6/6
```

Truncation keeps most of it: `pairHigh` error stock .153, full vine .059,
tree-1 .088, while tree-1 also has the **best** parameter error (.0424 vs stock
.0439, full .0461).

**Fitting a vine to EBEs/draws from a Gaussian fit is useless as a diagnostic**
(60+60 reps, N=100, best-case design, `out/power_N100_gumbel.rds`):
P(AIC selects >=1 non-Gaussian pair copula) = **1.000 under a Gaussian truth**,
at every level, every replicate. Calibrated power at 5%: 0.85 on true etas, 0.48
on conditional draws, 0.38 on EBEs. Do not use family selection as evidence.
(EBE caveat also backed by Savic & Karlsson 2009; Lavielle & Ribba 2016.)

### Refuted -- do not re-investigate

- **Pool truncation** as the cause of the nested-null gap: `poolMax` 40 vs 250
  gives 0.0735 vs 0.0729, i.e. no effect (`out/poolTrunc.rds`).
- **sd estimator source**: taking sd from saemix's exact SA sufficient statistic
  changes nothing (`pool-SS` 0.0711 vs `pool-v1` 0.0706, `out/arms.rds`).
- **M-step algorithm**: three unrelated machineries give an identical gap
  (pool 0.0706 / pool-SS 0.0711 / sa 0.0708 against stock 0.0437).
- **PIT misspecification** (posterior draws not Gaussian-marginal): real but
  negligible -- |skew| <= 0.14, |excess kurtosis| <= 0.22, sd ratio 1.00; a
  single slice shows nothing (Shapiro-Wilk p > 0.5). See `checkPIT.R`.

---

## 6. Gotchas

### rvinecopulib 0.7.3.1.0

- `selcrit="mbicv"` is a **silent no-op equal to `"bic"`** (fixed only on dev
  HEAD). Do not rely on it.
- mBICV's `psi0` sparsity dial is **inert at d<=6** (default 0.9 gives ~5-14% of
  the BIC term). Verified: at d=6, `trunc_lvl=NA` picks level 1 identically for
  `psi0` 0.9 and 0.2. You would need `psi0 < 0.5`, an unstudied regime.
- `trunc_lvl=NA` = automatic (always mBICV, ignores `selcrit`); `Inf` = none.
  **`VineCopula`'s `trunclevel=NA` means the opposite** (no truncation).
- Read the selected level from `dim(fit)["trunc_lvl"]`, **not**
  `fit$controls$trunc_lvl` (which keeps the user's input).
- Automatic selection **over-truncates at our n**: at d=4, n=100 a genuine
  tree-2 tau ~ 0.19 is missed 55-68% of the time. So do not claim "the data chose
  tree 1" -- it is partly the selector's floor.
- A truncated fit returns **fewer trees**; `copulaPadFlat()` pads the flat edge
  list with independence copulas to keep the canonical length and order.

### saemix internals

- `betas` must stay a **matrix** -- use `betas[] <- ...`, not `betas <- ...`.
- `t(COV) %*% COV` is **singular** (structurally zero columns). Apply a location
  shift to the **intercept row of `MCOV`**, not by least squares on `COV`. This
  is also the right semantics with covariates.
- Model scoring should use **n = number of SUBJECTS**, not observations
  (Delattre, Lavielle & Poursat 2014; Killiches & Czado 2018 Prop. 4.1).
- `varList$ind.eta` **changes at `kiter == nbiter.sa`** when some parameters have
  no IIV. The copula path asserts the dimension matches; if you support
  zero-IIV parameters you must handle this.

---

## 7. What to do next, in priority order

### P1 -- robustness by construction (this is the remaining core work)

1. **Test with covariates.** The top untested correctness risk. The
   Gaussian-surrogate GLS provably cannot move the `mu` fixed point for an
   *intercept-only* design (the GLS returns the column mean for any Omega), but
   **with covariates the weighting does not cancel**, and popPK models always
   have covariates. The intercept-shift path now exists, so this is testable.
2. **Treat `truncLvl = 1L` as a prespecified adaptive-phase candidate** for
   d >= 4, not a silent low-level default. It matches the literature, beat stock
   in the current experiment, and removed the collapse, but it also sets
   higher-tree Gaussian partial correlations to zero and therefore is not the
   unrestricted MVN nested null.
3. **Shrinkage gate.** After the stock warm-start fit, compute eta shrinkage and
   admit only etas below ~15% into the vine; leave the rest Gaussian. This
   targets the exact condition that breaks it (`out/decisive.rds`).
4. **Collapse guard.** Monitor the sd trajectory and |tau| during the run; if sd
   falls below a fraction of its warm start or |tau| pins near the cap, fall back
   to Gaussian and warn. Because the failure is *intermittent*, a fit-time
   detector matters more than a good average.
5. **Select the tree-1 structure once, then freeze.** At d=4 there are only
   d^(d-2) = 16 candidates (Cayley) -- enumerate, do not search. Measured loglik
   spread across the 16 is **3.4-6.7 units**, so the choice is not free.
   `runStructure.R` is **currently broken** (see §8).

### P2 -- usability

6. **Standard errors.** Now worth building: the joint fixed point solves the
   score, so observed-information SEs and LRTs are valid. saemix's own FIM
   machinery assumes Gaussian etas, so this needs a numerical observed FIM on
   the profiled likelihood, or a bootstrap.
7. **Fix the importance sampling.** `essMin` fell to 16/3000 for the copula
   marginal likelihood under a tail-dependent prior. `marginalLL()` already uses
   a defensive mixture (25% of draws from the prior), which fixed accuracy on a
   ground-truth case (1.0444 -> 0.9556 against saemix GQ 0.9512) but is still
   straining. Any likelihood-based model comparison depends on this.

### P3 -- breadth before believing it

8. Other truths: **Clayton** (lower tail -- the clinically dangerous corner for
   exposure, since AUC ~ 1/CL), **t** (both tails), strong-correlation Gaussian.
9. d=5, sparse designs, and eventually real data.
10. **Kernel 4** (f-SAEM MAP/IMH) under the copula path -- needs the copula
    Hessian. Currently refused. This is the E-step efficiency lever.

---

## 8. Known broken / incomplete

- **`runStructure.R` fails** with `object does not inherit from class
  vinecop_dist`, after printing the candidate count, inside the arm loop. It also
  builds **28** candidates instead of 16 (permutations of a D-vine path are
  counted twice -- reversing a path gives the same tree). Fix the dedup and trace
  the class error. The loglik-spread part of it works and is the source of the
  3.4-6.7 figure.
- `testJoint.R` reports Q on the **unshifted** draws, so `mode="joint"` looks
  slightly worse than it is once the location block is active. `/tmp/chk.R`-style
  direct checks (score at `E - delta`) are the correct comparison.
- `tests/test-nested-null.R` N5 splits by identifiability and is calibrated for
  the old `PK_*` model in `R/simData.R`, not the `simpleModels.R` ladder.
- `out/recovery.log` is empty (that run was killed as too slow; superseded by
  `runRecovery2.R`).
- `llcheck_*.R` were written by a subagent; `llcheck_core.R` holds `isPass`,
  `pack`/`unpack`, `vine2FromTau` and is the reusable part.

---

## 9. How to run things

```bash
export PATH="/c/Program Files/R/R-4.5.3/bin/x64:$PATH"
cd C:/package/saemix-copula/copula

Rscript tests/test-ordering.R        # exact, seconds -- run after ANY vine change
Rscript tests/test-nested-null.R     # end-to-end, minutes
Rscript testJoint.R                  # joint vs IFM M-step, unit level
Rscript testJointFit.R               # end-to-end smoke, all three modes
Rscript llcheck_joint.R 2 40 500 joint   # stationarity: reps, N, M, mode
Rscript llcheck_joint.R 2 40 500 ifm     # matched baseline
```

Experiment drivers are `run*.R`; results land in `out/`. Most take
`NREP N [family] [tag]` on the command line.

**The two tests to run after every change:** `test-ordering.R` (fast, exact) and
a Gaussian-vine nested-null fit (must reproduce stock saemix).

---

## 10. One-paragraph summary

An R-vine copula on the etas works. With `mode="joint"` it uses one common
likelihood-based population objective and the exploratory endpoints were
approximately stationary for that auxiliary surrogate (residual displacement
0.02-0.16%); this is not, by itself, an MLE certificate. Under a Gaussian truth it matches stock
saemix; under a non-Gaussian truth it corrects joint tail probabilities the
Gaussian model gets wrong by factors of 2-3, in every replicate. It has one
failure mode -- a weakly identified eta combined with a full vine -- which is
intermittent rather than systematic and is removed by truncating the vine after
tree 1, which also *beats* stock saemix on parameter error. The remaining work is
guarding that failure mode automatically, testing with covariates, and providing
standard errors.

---

## 11. Superseding publication state (22 August 2026)

The paragraph above records early exploratory work and is not the publication
claim. The current paper and frozen protocol are authoritative. The fixed
margin--vine construction has survived independent likelihood, Jacobian,
EM/Fisher, code-theory, historical-SAEM, and first-reader audits through round
30 without an open P1/P2 mathematical issue. Round 30 also closed the sole
generic custom-margin fail-closed defect by separating invalid optimizer trials
from the literal common-Q objective and exactly re-evaluating every accepted
endpoint. The high-level learner uses Normal
transformed-scale margins, makes at most one screened finite dependence move
from a Gaussian incumbent, then performs a newly initialized, warm-started
fixed-specification SAEM-MCMC fit; declared admissible continuous margins remain
lower-level fixed-model routes.

The prospective primary study contains 660 datasets (Gaussian/Gumbel truth;
N=50,100,300; 110 replicates per cell), with a 99-pair inference gate, exact
Gaussian fallback in the deployable workflow, a fixed-Gumbel oracle, and two
importance-sampling likelihood evaluations per returned arm. All 6,966 consumed
random streams are uniquely namespaced. The source-hashed Results pipeline now
maps validated `summary.rds` evidence to 12 main tables, 147-row Supplementary
Table S1, figures, lineage files, and atomic manifests; the manuscript has no
literal result placeholders and fails closed until those artifacts exist.

Two unrelated protected `ticbig.R` jobs (PIDs 30016 and 31568) still occupy
about 22 GB and have left less than 1 GB commit headroom. Do not stop them and
do not launch R until they exit or headroom remains safely above 3 GB
(preferably 5 GB). Once safe, run the source parsers and test suites, then the
three-cell 500+500 calibration gate, freeze the current source, execute the
660-cell 300+200 primary study, generate/audit Results, rewrite the empirical
Abstract/Discussion/Conclusion, and render/inspect the final PDF.

On 23 August, two additional protected `bigstudy.R` jobs (PIDs 23000 and
40004) were also active. Treat all four as user-owned and do not stop them.
Independent rounds 39--41 subsequently audited the implementation target, the
full publication-study analysis mapping, and saemix API compatibility. Their
post-repair verdicts are PASS with no open P1/P2 issue. In particular, all
common-Q optimizers now accept a stationary no-move without allowing an exact-Q
decrease; failed learner timing and likelihood-attempt counters are reported
with their correct denominators; and the public entry point distinguishes an
explicit isolated `population=NULL` from the historical omitted-argument
`copulaSet(); saemix(); copulaGet()` workflow. See `copula/audits/INDEX.md` and
the round-39--41 reports. These were static audits; R execution is still
required when the memory guard clears.

---

## 12. Runtime/profiling update (24 August 2026)

The user-requested long N=30, 500+500 selected-vine benchmark was interrupted
before it wrote results. No user-owned study process was stopped. The slowdown
was localized to the repeated joint population M-step over the growing
uncompressed empirical measure. In particular, `jointBackend="auto"` was
unnecessarily routing two-dimensional non-Gaussian vines to the literal R
objective even though the persistent compiled objective supports the same
fixed likelihood and L-BFGS-B parameterization.

Current changes:

- auto dispatch now uses the C++ objective for non-Gaussian `d >= 2`;
- `cores="auto"` resolves conservatively to at most four workers while the
  backward-compatible default remains `cores=1`;
- `profileFile`/`profileEvery` add an opt-in append-only CSV with current
  iteration, stage, empirical-pool size/ESS, backend, optimizer state, and
  E-step/joint-M-step/iteration wall times;
- repeated `do.call(rbind, rep(list(...)))` particle/design construction was
  replaced by equivalent row indexing;
- the conditional sampler now clamps round-off-negative online variances
  before `sqrt`, removing the anonymous `sqrt(varik)` NaN warning.

Verification: the R and C++ common-Q values agree to about `1e-15` in both 2D
and 4D tests; their L-BFGS-B fits reach the same Q; serial and four-core C++
evaluations agree exactly in the test. On a 6,000-row 2D pool, R/C++-1/C++-4
took 3.42/2.61/2.28 s. The complete 100-iteration N=30 quick workflow fell
from the earlier 127.6 s to 61.7--64.9 s and retained the same accepted Clayton
candidate and validation diagnostics. The final chain profiler identified the
joint M-step—not the E-step—as the dominant cost. See
`copula/out/profile-safe-speedups.csv`, `copula/profileSafeSpeedups.R`, and the
timestamped `copula/slides/data/saemix_profile_quick_*.csv` files.

---

## 13. Corrected covariate pitch benchmark and first 10x path (24 August 2026)

The covariate pitch benchmark is now rerun and publication-facing evidence no
longer uses the biased transition. The cause was a fresh `saemix()` process
incorrectly configured as `c(0,250)` with no burn-in: fresh chains, sufficient
statistics, and gain recursion were treated as a continuation, so an initial
residual-SD spike was retained by the decreasing gain. The corrected selected-
vine phase is a newly initialized warm-started `c(200,50)` fit with burn-in 5.
Its fail-closed Gaussian nesting gate passes (residual-SD relative error 0.12%;
conditional-omega relative error 4.4%).

Matched conditional marginal log likelihoods are -207.53 for prespecified
covariates, -207.32 for the exact automatic Gaussian collapse, -207.59 for the
independent direct Gaussian-vine validation, and -196.83 for the selected
R-vine. The selected-vine VPC quantile-curve RMSE is 0.084 versus 0.103 and
0.106. Arbitrary-subset conditional sampling was independently checked against
a logit-space Metropolis-Hastings implementation for weight only, eGFR only,
and both covariates; all checks pass and observed coordinates are exact.

Two exact Gaussian-collapse levels are now tested. For Normal eta margins,
Schur conditioning produces the ordinary-saemix random-effects distribution.
More generally, for any invertible continuous eta margins,
`z_e=Phi^-1(F_e(eta))` is conditionally Gaussian and
`eta_j=F_ej^-1(Phi(z_ej))` is a deterministic model transform. Thus arbitrary
eta and covariate margins retain the same marginal likelihood without generic
Gaussian-vine density evaluation. Student/Laplace eta density equality passes
to `3.38e-14`; the transport round trip passes to `2.11e-11` in
`copula/tests/test-gaussian-conditioning-collapse.R`.

For the non-Gaussian vine, Q-monotone generalized-EM intermediate blocks every
20 SAEM iterations plus the unchanged final simultaneous common-Q optimizer
reduce the fixed phase from 724.11 to 43.30 seconds (16.72x) and the full
workflow from 751.60 to 70.79 seconds (10.62x). Conditional marginal log
likelihood is -197.31 versus -196.83; the 0.48 difference is below the roughly
0.63 combined Monte Carlo uncertainty. VPC RMSE is 0.055 and residual SD is
0.10808 versus 0.10715. This is an accuracy-gated prototype, not yet a default.

The rebuilt 15-slide notebook/Reveal deck is
`copula/slides/copula-saem-pitch.ipynb` and
`copula/slides/copula-saem-pitch.slides.html`; all seven figures were rendered
and visually inspected, all image alt text is populated, and the discarded
-212.99 likelihood is absent. Reproduction and next collapse paths are tracked
in `copula/slides/TONIGHT-2026-08-24.md`. The main further exact routes are the
non-centred Rosenblatt transport, compatible Student-t scale mixtures,
Archimedean frailty augmentation, role-aware cancellation of covariate-only
vine factors, and partial Gaussian/linear-block collapse.

---

## 14. Exact Gaussian-FREM finite-statistic and conjugate collapses (25 August 2026)

Three additional exact paths now remove the remaining dominant costs for the
focused Gaussian-copula FREM model.

First, a frozen full-Gaussian model with Normal/lognormal margins is transformed
to its identity/log Gaussian scale. Within each missingness pattern, weighted
total, first and cross-moments of the transformed coordinates and linear design
row are complete sufficient statistics for the population common-Q. Thus even
continuous subject-specific designs do not create additional groups. In the
new varying-design regression gate, 600 rows reduce to 40 and the fitted beta,
correlation matrix and Q agree within `4.1e-14`. Those moments follow the same SA
recursion as the former unbounded empirical particle pool and are materialized
as at most `2p` symmetric pseudo-points for the existing Schur-EM M-step. Thus
the optimization sees exactly the same Q with constant memory. A 30,000-row
static check compressed to 24 rows with a Q difference of `1.03e-13`; a matched
saved 300-iteration fit changed fixed effects/correlations/log likelihood by
only about `1e-12` and fell from 211.18 to 9.56 seconds (22.1x).

The 30,000-row continuous-design benchmark reduces to 40 rows, takes 1.92
versus 0.22 seconds (8.73x for the standalone M-step), and differs by
`1.72e-13` in Q, `1.13e-14` in beta, and `7.75e-14` in the correlation matrix.
It is reproduced by
`copula/benchmarkGaussianFremAugmentedDesignCollapse.R`; results are in
`copula/out/gaussian_frem_augmented_design_collapse_benchmark.csv`.

Second, for an affine structural model with constant Gaussian error,
`eta | C_obs` is Normal by the existing Gaussian conditioning identity and the
individual posterior is conjugate Normal. `conddistCopula.saemix()` now returns
its exact moments and samples. The implementation verifies the affine identity
numerically and fails closed to Metropolis-Hastings for nonlinear or
heteroscedastic models. Exact posterior moments agree with an independent
analytic calculation within `2.1e-16`.

Third, the same conjugate identity analytically integrates eta in
`llisCopula.saemix()`. For the joint target it retains the separately evaluated
exact `f(C_obs)` term; diagnostics explicitly report
`method="exact-linear-gaussian"`, zero Monte Carlo error and infinite ESS.
Nonlinear PK models continue to use importance sampling.
The end-to-end regression test now fits both cases: it verifies the nonlinear
fallback and compares every subject's affine exact likelihood contribution to
a separately assembled Schur/conjugate Gaussian integral within `1e-10`.

Clean matched 1,000-iteration runtimes are 13.89/14.28 seconds for MCAR
Gaussian/automatic and 14.92/15.08 seconds for MAR. The corresponding old exact
raw-pool runs were about 1,458--1,469 seconds, a 97--102x reduction. Exact joint
observed log likelihood gains for automatic margins are +1.3263 (MCAR) and
+1.5309 (MAR). Automatic family learning itself remains about 215 seconds and
is now the main bottleneck: changing an arbitrary margin parameter changes all
normal scores, so the finite Gaussian moment identity does not apply. No lossy
screening shortcut was introduced.

Key implementation/tests:

- `R/gaussianCopulaFrem.R`: recursive moments, materialization, conditioning density;
- `R/copulaPrior.R`: SA-recursion integration and compressed M-step dispatch;
- `R/copulaLikelihood.R`: exact posterior and likelihood gate;
- `copula/benchmarkGaussianFremMomentCollapse.R`;
- `copula/benchmarkGaussianFremAugmentedDesignCollapse.R`;
- `copula/tests/test-gaussian-frem-mvnormal-mstep.R`;
- `copula/tests/test-gaussian-conditioning-collapse.R`;
- `copula/tests/test-gaussian-copula-frem-end-to-end.R`.

All listed regression tests plus `test-conditional-likelihood-target.R`,
`test-automatic-gaussian-frem-margins.R`, and
`test-automatic-gaussian-frem-workflow.R` pass on R 4.5.3, and a clean source
installation to `tmp/install-lib` succeeds. The
updated 19-slide pitch uses only the final 1,000-iteration exact
results and is in `copula/slides/copula-saem-pitch.ipynb` and
`copula/slides/copula-saem-pitch.slides.html`.

---

## 15. Matched 2,000-iteration check and FFEM translation (25 August 2026)

`run_nyberg_simple_1000.R` now accepts an optional total-iteration argument.
Matched 2,000-iteration runs use the same 400-iteration Gaussian start and an
uninterrupted 1,600-iteration fixed-model phase, preserving the nested random
stream through the first 1,000 iterations.

Exact observed log likelihoods at 2,000 iterations are:

- MCAR Gaussian -260.7718, automatic -259.4410, gain +1.3308;
- MAR Gaussian -264.4714, automatic -262.9882, gain +1.4832.

From 1,000 to 2,000 iterations, the automatic log likelihood improves by only
+0.0306 (MCAR) and +0.0810 (MAR); Gaussian FREM improves by +0.0261 and +0.1288.
Thus the flexible-margin advantage is essentially stable. Fixed/residual
parameters change by at most 0.00111, the largest correlation change is 0.0113,
and the range of any ordinary saemix trace column over the final 200 iterations
is below 0.00026. Recorded fixed-fit times are 36.27 seconds for the serial MCAR
Gaussian run and 45--47 seconds for three concurrent runs; the latter include
CPU contention, not additional statistical work.

Outputs are `*_2000.rds`, `*_2000.csv`, and `*_2000_parameter_trace.csv` in
`copula/examples/gaussian-copula-frem/out/nyberg_simple500/`, plus
`matched_2000_likelihood_comparison.csv`. The exact automatic length summary now
includes 2,000 iterations.

For communication with Nyberg, retain his distinction: GC-FREM is the
estimation model; an FFEM representation is derived afterward for GOF, VPC,
simulation and communication. Use the block notation
`Omega_FREM = (Omega_par, Omega_cov,par; Omega_par,cov, Omega_cov)`,
`beta_CP = Omega_par,cov Omega_cov^-1`, and conditional UPV
`Omega' = Omega_par - beta_CP Omega_cov,par`.
In GC-FREM these equations apply exactly to the latent Normal scores. With
Normal parameter margins they translate exactly to the parameter random-effect
scale; lognormal covariates enter through standardized log values. With
non-Normal parameter margins, report conditional medians/quantiles rather than
claiming an exact standard additive-Gaussian FFEM. The full derivation is now in
`copula/FREM-DESIGN.md`, section 2.1.

---

## 16. Exact partial collapse for one arbitrary margin (25 August 2026)

`R/gaussianCopulaFrem.R` now profiles a complete-data Gaussian-copula FREM
M-step exactly when precisely one margin is not Normal-Gaussianizable. If `a`
is that coordinate, its native density and normal score
`z_a = Phi^-1(F_a(x_a))` remain inside the nonlinear optimization. Conditional
on each trial value, all remaining transformed coordinates satisfy
`Y_G = alpha + B z_a + epsilon`, with Gaussian `epsilon`; weighted regression
therefore gives the exact MLEs of `alpha`, `B`, and `Psi`. The full covariance
is reconstructed as `Sigma_GG = Psi + B B'`, `Sigma_Ga = B`, so every Gaussian
scale and every copula dependence parameter is removed from the nonlinear
optimizer. The implementation independently re-evaluates the resulting model
through the original literal joint density and fails closed on disagreement.

The path currently supports no design and the ordinary full-rank constant
eta-intercept design. It deliberately falls back for missing data, general
covariate designs, conditional-likelihood targets, or two or more arbitrary
margins. For an arbitrary observed covariate, repeated particle rows are
compressed exactly to group masses, weighted Gaussian-coordinate sums, and
cross-products; cost therefore does not grow with the retained SA particle
history merely because the covariate repeats.

On the controlled four-dimensional Student-margin test, the all-Gaussian M-step
takes 0.14 seconds, the exact one-Student profile 0.11 seconds, and the legacy
literal optimizer 39.86 seconds. Profile and literal common-Q values are
-8.5476268953 and -8.5476286375 (difference 1.74e-6, with the profile slightly
higher); independent literal evaluation of the profiled result agrees within
1e-8. `copula/tests/test-gaussian-frem-partial-collapse.R` also proves that the
intercept-design and no-design parameterizations give the same objective and
locations.

An automatic two-eta example with a strongly skewed Gamma covariate selects
`normal/normal/gamma`. Timings for 350 iterations are 1.04 seconds for ordinary
saemix, 8.27 seconds for all-Normal Gaussian-copula FREM, and 18.11 seconds for
Gamma-margin Gaussian-copula FREM; the Gamma overhead relative to the same
copula model is therefore 2.19x, rather than the roughly 250--360x literal
M-step penalty. Automatic learning and validation take 4.73 seconds. The fresh
Gamma fit reports backend
`gaussian-copula-frem-one-arbitrary-profile`. Reproduction script and outputs:
`copula/examples/gaussian-copula-frem/run_heavytail_partial_collapse.R` and
`copula/examples/gaussian-copula-frem/out/heavytail_partial/`.

Next exact extensions are: condition the Gaussian block on multiple arbitrary
scores (leaving only their small internal correlation block nonlinear), extend
the pattern-wise Gaussian EM to missing coordinates, and add
transformation-Normal margins for which only a low-dimensional transform shape
must be profiled. Separately fitting then freezing arbitrary margins would be
IFM rather than the joint MLE and must not replace these exact paths.

---

## 17. Replicated saemix/NONMEM FREM validation (25 August 2026)

The complete reproducible study is in `copula/study_frem_validation/`. It uses
four frozen scenarios (`normal/normal`, `lognormal/normal`, `gamma/normal`, and
`weibull/lognormal` covariate margins), 10 independently simulated datasets per
scenario, N=100, and 1,500 SAEM iterations per primary fit (500 constant-gain +
1,000 decreasing-gain; two MCMC samples). Every dataset was fitted with:

1. standard full-block Gaussian FREM in NONMEM on klebsiella;
2. standard Gaussian FREM in saemix; and
3. Gaussian-copula FREM in saemix with the true covariate family specified and
   all continuous marginal parameters jointly estimated.

All 120 primary fits completed, as did 40 additional NONMEM EONLY evaluations
at the exact fitted saemix Gaussian-FREM parameter points. The raw/full
likelihood conversion is

`OFV_raw_full = OFV_NONMEM + n_DV*log(2*pi) - 2*log(J_raw_to_score)`.

Common AIC/BIC penalties are added only after this conversion. The same-point
NONMEM-minus-saemix deviance means are -0.094, +0.054, +0.060, and -0.015 for
the gamma, lognormal, normal, and Weibull scenarios; all scenario-level RMSEs
are 0.22--0.44 versus combined likelihood MCSEs of about 0.34. Thus the
independently implemented likelihoods agree within integration error. The
2.5--4.0 deviance differences between independently optimized NONMEM and
saemix fits reflect stochastic fitted endpoints, not an OFV conversion error.

Mean flexible-minus-Gaussian log-likelihood gains are +4.279 (gamma), +5.693
(lognormal), -0.062 (normal), and +12.334 (Weibull/lognormal). The flexible fit
wins all 10 non-Normal replicates and shows essentially no gain under the
nested Normal null. Structural and free-margin recovery are reported in
`RESULTS.md`. Current exact one-arbitrary-margin fits are roughly ten times
slower than the Gaussianizable collapse in this 1,500-iteration study; this is
the next implementation target, not a reason to weaken the likelihood target.

Primary artifacts include `out/all_results.csv`, `out/paired_results.csv`,
`out/recovery_summary.csv`, `out/margin_parameter_summary.csv`,
`out/same_parameter_likelihood.csv`, and `out/validation_summary.png`.
The archived primary NONMEM outputs are `out/nonmem_results_named.tar.gz`
(SHA256 `8A4318AA4DF08E9468B1201D26E15CA832315DF220B9D73CE12C48A6335E6F7F`);
same-point controls are `out/fixed_results_named.tar.gz`
(SHA256 `FB87150201FA4FA5D8FEC60413836917BA048FE989B5850DCF21A56CD961B5B9`).

The Jupyter/Reveal pitch deck was rebuilt from these results on 26 August 2026.
It now has 20 slides and separates the single-dataset MCAR/MAR workflow example
from the 40-dataset validation. New evidence slides show the replicated
flexible-margin likelihood gain, same-parameter NONMEM/saemix likelihood audit,
complete 1,500-iteration wall times, and replicated V/CL recovery. The timing
slide states explicitly that NONMEM ran on different hardware and that the
meaningful within-saemix comparison is approximately 46--47 seconds for
Gaussianizable flexible margins versus 438--450 seconds for the current exact
Gamma/Weibull profile. Sources are `copula/slides/copula-saem-pitch.ipynb` and
`copula/slides/gaussian_frem_pitch_content.py`; the ready-to-present export is
`copula/slides/copula-saem-pitch.slides.html`.

On 26 August, the pitch was extended to 22 slides with separate explanations
of the full Gaussian and one-arbitrary-margin partial Gaussian collapses. The
validation recovery slide now makes the result boundary explicit: mean
absolute relative V/CL bias is essentially identical (0.628% Gaussian versus
0.625% flexible), while the flexible model recovers lognormal, Gamma and
Weibull shape/scale parameters that do not exist in standard Gaussian FREM.

The local installation was also verified. PsN 5.4.0 was already installed at
`C:/PsN-5.4.0` but absent from PATH; NONMEM 7.6.0 is at `C:/nm76g64`. A complete
1,500-iteration `normal_normal_r01` run was executed through PsN inside the
NONMEM tree. The server used NONMEM 7.5.1, so independent SAEM endpoints differ
slightly (V 9.77960 vs 9.76285; CL 2.93194 vs 2.92990; server vs local).
At the exact fixed saemix point, however, all nine IMP evaluations agree to a
maximum 7.3e-9 and the mean OFV differs by 2.1e-9. Preserved results are
`out/local_nonmem_check/`, `out/local_nonmem_likelihood_comparison.csv`, and
`out/local_nonmem_primary_parameter_comparison.csv`; reproduction parsing is
in `summarize_local_nonmem_check.R`.

On 26 August the opening of the pitch was reorganized into a 24-slide
layman-first mathematical sequence: observations and latent quantities,
observed marginal likelihood, standard Gaussian FREM in Nyberg block notation,
Gaussian-copula FREM with arbitrary margins, and then full/partial Gaussian
collapse. The downstream robustness slide no longer uses paired RMSE dots. It
now shows conditional VPCs for the highest COV1 decile in all three non-Normal
validation scenarios, with the generating 10th/50th/90th percentiles overlaid
against standard and flexible FREM. Mean VPC-quantile RMSE decreases from 0.126
to 0.050 and improves in all three scenarios. This is explicitly described as
robustness to marginal-shape misspecification, not generic outlier robustness.
Reproduction is `study_frem_validation/make_tail_conditional_vpc.R`; outputs are
`out/tail_conditional_vpc.csv` and `out/tail_conditional_vpc_summary.csv`.

---

## 18. Envelope scores for non-Gaussian partial collapse (26 August 2026)

The exact one-arbitrary-margin profile is distribution-generic: it applies to
any declared continuous margin when all other coordinates have fixed
Gaussianizing transforms. The margin is not approximated or frozen. Its native
log density and `Phi^-1(F(x))` remain in the same joint common-Q objective;
only the conditional Gaussian block is maximized algebraically.

`R/gaussianCopulaFrem.R` now additionally accelerates the remaining small
optimizer with envelope scores for Gamma, Weibull, Laplace, and Student-t.
Weibull and Laplace scores are analytic. Gamma shape and Student-t degrees of
freedom each require one centered vectorized derivative because base R exposes
no corresponding CDF derivative; their other coordinates are analytic. The
optimizer automatically retries the identical numerical-gradient objective if
an analytic-score trial fails. Unsupported/custom continuous margins still use
the same exact partial collapse, with the numerical gradient.

`copula/tests/test-one-arbitrary-envelope-gradients.R` compares each accelerated
fit against the identical profile optimized numerically. Objective differences
are at most `6.1e-8`, arbitrary-margin parameter differences are at most
`6.3e-5`, and independently evaluated literal joint Q agrees within `1e-8`.
In the N=500 focused M-step test, observed speedups were approximately 2.1x
(Gamma), 2.2x (Weibull), 4.0x (Laplace), and 1.7x (Student-t). These are
optimizer-only timings, not claims about full SAEM wall time.

The exact boundary remains important: with two or more arbitrary margins the
current dispatcher uses the literal common-Q optimizer. The next algebraic
extension is to condition the Gaussian block on their vector of normal scores,
leaving only the small arbitrary-margin/dependence block nonlinear.

Matched complete-data 350-iteration smoke fits were also run through the full
saemix workflow (N=150, identical PK model and schedule). Gamma-margin GC-FREM
took 16.49 s versus 16.95 s for its Gaussian-FREM comparator; Weibull-margin
GC-FREM took 15.01 s versus 13.82 s. The flexible backend was
`gaussian-copula-frem-one-arbitrary-profile` in both cases. In the high-COV1
conditional VPC, RMSE changed from 0.212 to 0.111 for Gamma and from 0.234 to
0.131 for Weibull. These are implementation smoke checks, not a replacement
for the frozen 1,500-iteration replicated validation. Reproduction outputs are
`copula/slides/data/nonaffine_complete_vpc_gamma_350iter_benchmark.rds` and
`copula/slides/data/nonaffine_complete_vpc_weibull_350iter_benchmark.rds`.

---

## 19. Exact partial collapse for multiple arbitrary margins (26 August 2026)

`R/gaussianCopulaFrem.R` now generalizes the exact profile from one curved
coordinate to a curved block `A` of size two or larger, provided a nonempty
Gaussianizable block `G` remains. The nonlinear optimizer retains the native
margin parameters and eta locations in `A`, plus its internal Gaussian
correlation `R_AA`. For every trial point it forms `Z_A` and profiles

`Y_G = alpha + Z_A K + epsilon`, with `epsilon ~ N(0, Psi)`,

by weighted multivariate regression. It reconstructs
`Sigma_GA = K' R_AA` and
`Sigma_GG = Psi + K' R_AA K`, then maps the resulting full correlation matrix
back to the original Gaussian R-vine structure. Thus all Gaussian-block means,
scales, within-block dependence and cross-block dependence leave the nonlinear
optimizer without altering the joint common-Q objective.

`copula/tests/test-gaussian-frem-multi-arbitrary-collapse.R` verifies the public
automatic dispatcher, the ordinary eta-intercept design, and independent
literal-density identity for both two- and three-curved-coordinate examples.
For Student-t plus Laplace margins at N=450, the profile took 1.14 s and reached
Q=-8.6686736539; the former literal all-parameter optimizer took 24.41 s and
reached Q=-8.6687023272. The profile therefore ran 21.4x faster and attained a
slightly higher finite-iteration Q. Its independently evaluated literal Q
agreed within 1e-8; the three-curved example agreed exactly at printed
precision.

The implementation still fails closed for missing coordinates, general
covariate location designs, conditional-likelihood fitting, or when every
coordinate is curved (there is then no Gaussian block to profile). The next
speed step is a joint envelope score for the remaining curved block and exact
recursive compression where repeated observed curved-coordinate patterns make
that possible.

---

## 20. Closing the function-valued-recursion gap (26 August 2026)

The most general fixed margin--vine model no longer relies, theoretically, on a
changing weighted distribution. If finite sufficient statistics exist, use the
ordinary MCMC-SAEM recursion. Otherwise update the finite parameter vector with
the posterior average of the complete-data score (Gu--Kong/Gu--Li MCMC score
stochastic approximation). Fisher's identity makes its mean field exactly the
observed-data marginal-likelihood score. Under stable iterates, Robbins--Monro
gains, uniformly controlled parameter-indexed MCMC kernels, score moment and
Poisson-equation bounds, a bounded positive-definite metric, and summable
numerical bias, the limiting ODE is

`theta_dot = A(theta) gradient ell(theta)`.

Because `gradient ell' A gradient ell > 0` away from the stationary set, the
observed log likelihood is a strict Lyapunov function. Controlled-Markov
stochastic-approximation theory then gives almost-sure convergence in distance
to the stationary set of the exact observed likelihood. Isolated stable points
are local maxima; no global maximum is promised.

Authoritative details are in `copula/SCORE-SAEM-THEORY.md`,
`copula/MARGINAL-LIKELIHOOD-PROOF.md`, the article's fixed-specification section,
and Supplementary Note A. Fisher-identity oracle tests pass for a fixed effect,
Gamma shape, native copula parameter, and residual parameter. The important
implementation boundary is explicit: current collapsed validation routes use
finite-statistic MCMC-SAEM; the all-parameter score backend still has to be
integrated and qualified before the fully general package route is claimed.

---

## 21. Exact all-marginal Gaussian-copula ECM/profile (28 August 2026)

`R/gaussianCopulaFrem.R` now covers the case that defeated the earlier partial
collapse: every coordinate may have a non-affine continuous margin, leaving no
Normal/lognormal raw-coordinate block. For fixed Gaussian correlation
`R`, coordinate `j` is updated from

`log f_j + 0.5(1-Omega_jj) z_j^2 - z_j Z_-j Omega_-j,j`,

where `Omega=solve(R)` and `z_j=qnorm(F_j)`. Thus each marginal block uses only
its own density/CDF and one cached projection. With margins fixed, all particle
dependence needed for the correlation update collapses to `S=Z' W Z`. The
correlation block now uses a one-to-one hyperspherical Cholesky
parameterization and analytic gradient, then maps the fitted correlation back
to the frozen user R-vine structure once. Every accepted block is
nondecreasing on the literal common Q.

Observed non-Gaussian covariates are additionally grouped by exact value. Their
mass and weighted cached projection are sufficient for the marginal block, so
Gamma, Weibull, beta, or custom margin evaluations scale with distinct subjects
rather than retained conditional particles.

The new dispatcher backend is
`gaussian-copula-frem-marginal-ecm`. The oracle in
`copula/tests/test-gaussian-frem-all-marginal-ecm.R` combines Student and
Laplace eta margins, Gamma and Weibull covariates, and a generic beta margin.
It verifies the analytic correlation gradient, literal-density identity,
ordinary eta-intercept design, covariate grouping, and comparison with the old
literal all-parameter optimizer. On isolated N=320, d=5 repeated runs, runtime
fell from 50.24--73.47 s to 0.67--1.00 s (50--110x); common Q was slightly higher at the finite
optimizer limits (`-9.9877432` versus `-9.9880707`). Existing MVN, one-curved,
multi-curved, automatic-margin, missingness, exact-numerics, and hybrid tests
continue to pass.

Automatic order is now: full Gaussian moment collapse; one-curved exact
profile; multi-curved partial Gaussian profile; all-marginal ECM/profile; then
literal fallback. Ordinary SAEM iterations use one exact nondecreasing ECM
sweep; the exploration-to-smoothing transition and final polish use up to 12
sweeps. Remaining speed work is analytic native marginal-score
callbacks for the small multi/all-curved blocks and integration with the
finite-dimensional general score backend.

---

## 22. Presentation-fit conditional-kernel and sparse-update validation (28 August 2026)

The all-marginal ECM does not activate in the main presentation example because
that model has Normal eta margins, one Gamma/Weibull covariate, and one Normal
covariate; it uses `gaussian-copula-frem-one-arbitrary-profile`. Two exact
coding collapses were added for this route:

1. the E-step now constructs `p(eta | c_obs)` once and reuses it for both
   conditional-prior proposals and all random-walk MH ratios. Omitting
   `p(c_obs)` from a within-row ratio is exact cancellation. Tests GC12--GC13
   verify equality with joint-prior MH differences through `1.8e-15`, including
   arbitrary eta margins and missing-covariate patterns;
2. a validated full Gaussian R-vine is copied and its native partial-
   correlation matrices are replaced directly rather than reconstructing every
   `bicop_dist`. GR4b verifies correlation, native density and seeded sampling
   against the constructor through `2e-12`.

On the matched 350-iteration Gamma workflow these changes gave 8.20 s for the
flexible fit and 8.23 s for the Gaussian reference, with identical parameters
and VPCs. On the more variable 1,500-iteration run the every-iteration flexible
fit fell from the original 168.92 s to 156.94 s; its E-step fell from 90.67 to
47.96 s, leaving the repeated M-step dominant.

The existing public `refitEvery=5` control was then validated without changing
the package default. At 1,500 iterations:

- Gamma: 156.94 -> 56.18 s (2.79x); 299 versus 1,495 population updates;
  observed LL -2279.127 versus -2279.126 with MCSE about 0.22; maximum typical
  parameter difference 0.00161 and covariance difference 7.3e-5.
- Weibull: 119.72 -> 59.25 s (2.02x, parallel timing); observed LL -2306.860
  versus -2306.883 with MCSE about 0.22; maximum typical parameter difference
  0.00174 and covariance difference 1.4e-4.

Tail-VPC RMSE changes were about 0.001--0.0014. Two-thousand-iteration sparse
runs took about 84--86 s under concurrent execution and did not systematically
move closer to the finite every-iteration endpoints. Therefore the current
recommended accelerated presentation setting is 1,500 iterations with
`refitEvery=5`, a mandatory final full update, and likelihood/VPC gates.
`refitEvery=1` remains the standard compatibility/reference setting.

---

## 23. Formal arbitrary-margin likelihood/SAEM proof (28 August 2026)

`copula/FORMAL-GAUSSIAN-COPULA-SAEM-PROOF.md` is now the focused formal proof.
It follows the notation and claim boundaries of Delyon, Lavielle and Moulines
(1999) and Kuhn and Lavielle (2004):

- arbitrary proper absolutely continuous margins plus a positive-definite
  Gaussian copula normalize by the PIT/normal-score change of variables;
- the response density times that population density defines the ordinary
  complete law, whose latent integral is the observed-data marginal likelihood;
- the exact EM/KL and Fisher identities refer to this same likelihood;
- arbitrary *observed FREM covariate* margins retain a finite curved-
  exponential statistic because their normal scores are observed coefficient
  functions; missing-data dependence enters only through per-subject/pattern
  masses, first moments and cross-moments of Gaussianizable latent coordinates;
- exact recursive compression and conditional-prior MH cancellation do not
  alter the statistic, posterior, or likelihood;
- Delyon Theorem 5 plus Kuhn--Lavielle applies to that finite-statistic route
  under their M1--M5, SAEM1--SAEM4, stability, and Markov-kernel assumptions;
  Delyon Section 8.1/Theorem 8 covers a fixed iterative M map under MAX1--MAX2;
- arbitrary *latent individual-parameter* margins need not have finite
  sufficient statistics. Their theorem-backed general route is Delyon Section
  8.2's stochastic-gradient construction plus MCMC control, whose mean field is
  the exact observed score by Fisher's identity.

The proof explicitly corrects a potential overclaim about `refitEvery=5`.
Skipping only one block does not change the likelihood target, but the current
time-varying partial-delay schedule is not covered by Delyon Theorem 8 merely
because the delay is bounded. A fully macro-blocked construction is proved;
the implemented partial delay remains a numerically gated acceleration pending
an asynchronous-SA/augmented-state proof. `refitEvery=1` is the theorem-aligned
reference.

`copula/tests/test-formal-saem-proof-obligations.R` passes EM/KL identities for
Normal, Student, Laplace, lognormal, Gamma, Weibull, beta, and a custom logistic
margin (maximum error 2.84e-12), conditional/joint MH cancellation (5.55e-16),
recursive moment identity (4.44e-16), blocked-gain algebra (5.55e-17), and
profile-versus-literal common-Q equality (zero at printed precision).

The post-proof scope audit adds three explicit boundaries. First, Delyon
Theorem 8 is used only conditionally: MAX1 requires a unique global maximizer of
the finite-statistic auxiliary function and MAX2 requires geometric convergence
of one fixed iterative map; ordinary optimizer success or common-Q
nondecrease alone is insufficient. Second, population normalization and the
EM/Fisher identities hold for every fixed admissible R-vine, whereas the
quadratic finite-statistic collapse is Gaussian-copula-specific; general fixed
differentiable R-vines use the score route unless another finite statistic is
proved. Third, "arbitrary" means absolutely continuous fixed-support margins.
Conventional FREM represents dichotomous covariates as near-error-free numeric
responses and polychotomous covariates by dummy variables, which is not a true
Bernoulli/multinomial population law and can impute intermediate category
values when data are missing. Proper discrete copula-FREM remains a separate
rectangle-probability/latent-threshold implementation and proof task.

---

## 24. Proper categorical Gaussian-copula FREM (28 August 2026)

Observed categorical covariates no longer use the continuous Gaussian FREM
approximation. `R/copulaMargins.R` adds `copulaMarginBernoulli()`,
`copulaMarginOrdinal()`, and `copulaMarginCategorical()`. Binary/ordinal
categories are latent standard-Normal intervals with stick-breaking-logit
probabilities; nominal categories require an explicit latent order, which is
stored as part of the model fingerprint rather than silently inferred.

`R/gaussianCopulaFrem.R` now evaluates the proper mixed density/mass as the
continuous Gaussian-score density times a conditional Normal rectangle
probability. One discrete coordinate uses a Normal CDF difference; multiple
coordinates use deterministic `mvtnorm::Miwa` rectangle integration. Missing
categories are marginalized by omitting their latent coordinates. The exact
conditional kernel samples categorical latent scores from their truncated
Normal rectangle and then eta from the Gaussian conditional distribution.
Random-walk prior ratios use the corresponding selection-Normal mass, so the
category-pattern normalizer cancels exactly.

The Gaussian FREM dispatcher refuses continuous partial-collapse formulas when
any margin is discrete and uses the literal mixed common-Q objective. The
one-arbitrary continuous recursive statistic is also disabled for discrete
conditioning, preventing loss of stochastic history. `copulaBuildNative()` now
preserves mixed `var_types`, and copula likelihood importance sampling permits
discrete conditioning while continuing to reject discrete latent eta margins.

`copula/tests/test-gaussian-copula-frem-categorical.R` passes Bernoulli and
three-level ordinal rectangle oracles (8.88e-16), exact missing-category
marginalization, conditional/joint MH cancellation, selection-Normal moments,
public missing-category sampling, joint Bernoulli/dependence recovery, explicit-
order nominal normalization (2.22e-16), and multiple-category rectangle
normalization (2.22e-16). The nonlinear end-to-end SAEM and likelihood test in
`test-gaussian-copula-frem-categorical-end-to-end.R` also passes.

The theorem boundary remains explicit. The mixed observed likelihood is proper,
but after integrating categorical latent scores its log rectangle probability
is generally nonlinear in latent eta and need not have the continuous Gaussian
finite statistic. The current retained-Q categorical fit therefore has the
correct target but is not claimed under Delyon Theorem 5/8. A theorem-backed
production convergence result requires either carrying accepted categorical
latent scores with a threshold-update proof or using the general MCMC score
backend. Order-invariant nominal covariates additionally require a future
multinomial-probit utility block; the implemented scalar nominal route requires
an explicit scientifically accepted latent order.

---

## 25. Delyon Theorem 8 certificate and backend audit (28 August 2026)

`copula/THEOREM8-CERTIFICATE.md` now reproduces the operative content of MAX1
and MAX2 and maps every fitting backend to it. The main correction is that
Theorem 8 replaces exact M-step condition M5 with substantially stronger
conditions than generalized-EM ascent: every fixed auxiliary function must
have a unique differentiable global maximizer, and repeated applications of
one fixed continuous M-map must converge to it geometrically, uniformly on
compact sets.

A constructive sufficient certificate was added. If

`m I <= - Hessian_theta L(s,theta) <= M I`

uniformly on the declared statistic and convex parameter domains, the fixed
projected-gradient map with `0 < step < 2/M` satisfies MAX1--MAX2 with rate
`max(abs(1-step*m), abs(1-step*M))`. The Hessian bounds must be analytical or
interval-certified; a numerical grid, common-Q gate, optimizer status, or
small gradient is not a proof.

`R/theorem8.R` implements the contraction-rate algebra, fixed projected map,
and finite-trajectory falsification audit. The audit explicitly labels itself
as numerical rather than a proof. `copula/tests/test-theorem8-fixed-map.R`
passes determinism, common-objective monotonicity, the analytical geometric
bound, convergence, and rejection of a noncontractive step.

The resulting backend classification is deliberately conservative:

- a fully observed unconstrained Gaussian solve has an exact M-step;
- missing-pattern Gaussian EM, constrained Gaussian ECM, arbitrary-margin
  profiles, and all-marginal cyclic ECM need family/design-specific MAX1--MAX2
  certificates before an unconditional Theorem 8 claim;
- categorical retained-Q, arbitrary latent margins, and general R-vines remain
  augmented-state or score-SA problems;
- `refitEvery>1`, phase-dependent optimizer controls, fallback guards, and
  adaptive model changes are not one time-homogeneous Theorem 8 map.

The main paper, Online Resource, and standalone proof now use this narrower
statement. `refitEvery=1` with a frozen model is necessary for the reference
path, but is not by itself sufficient: the actual nonlinear M-map still needs
an exact solve or a MAX1--MAX2 certificate.

The collapse implementation now exposes this classification on every returned
Gaussian FREM M-step as `result$theorem8`, and fitted state retains it under
`copulaGet(fit)$lastJoint$theorem8`. Exact objective reduction, exact M-step,
fixed-map status, MAX1/MAX2 certification, and the applicable theorem route
are separate fields; an uncertified nonlinear optimizer is never labelled
theorem-aligned.

`augmentMissingGaussian=TRUE` is the new default for a frozen Gaussian copula
with identity/log-Gaussianizable continuous margins. Given an accepted eta
draw, `copulaGaussianFremImputeMissingConditioning()` draws missing covariate
normal scores from their exact conditional law and maps them back through the
declared margins. The SA statistic is therefore complete Gaussian moments.
With the saturated FREM intercept design the following mean/covariance M-step
is exact, returning this missing-covariate path to M5/Theorem 5. Setting the
option false retains the earlier analytic-marginalization/nested-EM route.

`copula/probeTheorem8OneMarginCurvature.R` falsified a tempting alternative for
non-Gaussian margins. On realistic 5x5 log-shape/log-scale boxes, the exact
Gamma profile had negative curvature at 10/25 points; Weibull had negative
curvature at 5/25 and four invalid boundary points. The smallest Weibull
curvature at its fitted optimum was only 3.6e-8. A robust global
strong-concavity certificate is therefore unavailable for these examples.
Their collapse remains exact, but an unconditional convergence result should
use another proof or the score-SA route.

`copula/tests/test-collapse-theorem-status.R` verifies exact versus conditional
backend labels, conditional-augmentation moments, automatic pool augmentation,
and restoration of the exact Gaussian M-step with both implicit and ordinary
intercept designs. Theorem status and end-to-end missing-covariate regression
tests pass.

---

## 26. Score-SA oracle experiments (28 August 2026)

The score route is no longer theory only. Three standalone experiments are
documented in `copula/SCORE-SAEM-THEORY.md` and write reproducible CSV/PNG
artifacts under `copula/out`.

`experimentProfileScoreMargins.R` first tested deterministic ascent on the
exact one-arbitrary-margin profiled auxiliary function. Sensible starts reached
the optimizer reference, but opposing shape/scale starts exposed flat or
inferior regions (one Gamma gap 0.0351; one invalid Weibull path). This confirms
that the nonconcavity is real and that profile-score ascent cannot manufacture
MAX1.

`experimentScoreSaMargins.R` implements Delyon equation (74) directly for one
latent eta and one observed Gamma/Weibull covariate. Eta is redrawn from its
exact conditional distribution at every iteration. Two starts per family
converged to independently optimized observed marginal likelihoods: Gamma gaps
6.5e-7 and 1.2e-8; Weibull gaps 3.6e-6 and 2.3e-7. Maximum relative parameter
difference was below 0.003.

`experimentScoreSaIncomplete.R` adds a genuine incomplete-data response model,
`y = eta + error`. Both the eta posterior and observed marginal likelihood are
analytic, providing a strict oracle for Fisher's identity. Two starts per
family again converged: Gamma observed-likelihood gaps 4.3e-7 and 4.5e-8;
Weibull gaps 4.1e-7 and 7.6e-8. Maximum relative parameter difference was
0.00110 (Gamma) and 0.000582 (Weibull). The convergence plot is
`copula/out/score_sa_incomplete_trace.png`.

The experiments use projected log-shape/log-scale coordinates, decreasing
Robbins--Monro gains, a bounded positive diagonal running score scale, exact
conditional simulation, and Polyak averaging. They establish feasibility and
the correct likelihood target in oracle models, not yet the full production
NLME proof. The next coding step is a fixed-model population score backend with
one stable unconstrained vector covering locations, marginal parameters,
Gaussian dependence, and scales; it should be activated after a clean restart
and compared with the current common-Q fit on the validation study. Section 27
records that implementation and its first comparison.

---

## 27. Package-level full score-SA backend (28 August 2026)

`R/copulaScoreSa.R` and `populationAlgorithm="score-sa"` now provide an
experimental production path through the real saemix loop. The default remains
`"common-q"`. The score coordinate contains every population-location
coefficient, every free continuous margin parameter, every Gaussian dependence
parameter, and every free residual-error parameter. Population and response
scores use the same imputed individual state. Estimated structural parameters
outside the population location block are rejected rather than silently mixed
with an unproved update.

The API enforces a theorem-oriented contract: full fixed Gaussian vine,
continuous transformed-additive margins, joint likelihood, `activeFrom=1`,
`fitFrom=1`, `refitEvery=1`, `warmStartOnActivate=FALSE`, no family/structure
learning, and no fallback fingerprint change. Model learning belongs before a
fresh score fit.

The diagonal positive metric adapts only during a finite initial segment and
is then frozen, with fixed positive eigenvalue
bounds. Numerical complete scores use a centered step
`h_k = h0 * sqrt(gamma_k)`. Under bounded third derivatives the bias is
`O(h_k^2)`, so its gain-weighted sum is bounded by a constant times
`sum(gamma_k^2)`. Polyak averaging starts only in the decreasing phase and is
used for terminal reporting. Projection, metric, score, gain, difference-step,
and theory-eligibility diagnostics are retained in the fitted copula state.

After the first four-scenario run showed that 500 full-size score steps could
drive a lognormal `sdlog` from 0.24 to about 1.05, the score gain was separated
from the stock SAEM statistic gain. The API now uses `scoreBurn` finite
metric-learning steps followed by
`scoreGainScale * (k-scoreBurn+1)^(-scoreGainPower)`, with power restricted to
`(0.5,1]`. This is the Robbins--Monro sequence used by the score proof and
prevents a long unit-gain random walk in weakly identified margin directions.

Fixed complete-data full-vector experiments reproduce exact-profile fits after
model-compatible initialization: Gamma/Weibull common-Q gaps 7.0e-8/1.2e-7,
margin differences below 0.00054, correlation differences below 0.00025, and
zero projections. Rough declarations can converge to inferior stationary
points, as expected for a local likelihood method; observed-margin and current
Gaussian-score initialization is therefore part of the workflow.

New passing tests:

- `test-score-sa-population-backend.R`: stateful dispatcher, coordinate
  transforms, frozen metric, no projection, and common-Q gap 1.9e-7;
- `test-score-sa-fisher-oracle.R`: posterior mean complete score agrees with
  the observed score (max error 0.00127; max MCSE 0.00243);
- `test-score-sa-end-to-end.R`: actual saemix E/M loop including population,
  residual, fixed effects, and terminal averaging.

The first nonlinear PK comparison used 500 score iterations versus 150
common-Q iterations. Gamma margins were 1.633/39.04 versus 1.640/38.88;
observed log likelihood differed by 0.42 with combined IS uncertainty about
0.78; neither route showed a meaningful likelihood disadvantage. Score-SA had
zero projections. Analytic scores were subsequently added for Gaussian
dependence, Normal, lognormal, Weibull, Laplace, population locations, and
residual errors. Gamma retains only its shape-CDF derivative on the summable
shrinking-difference schedule. `test-score-sa-analytic-gradient.R` agrees with
an independent all-coordinate numerical gradient to about 1.7e-10 for Gamma
and Weibull. Runtime fell from about 28 to 8.7 seconds versus 2.6 seconds for
common-Q.

---

## 28. Analytic and categorical score extensions (28 August 2026)

The full Gaussian dependence block now uses hyperspherical correlation angles
and its analytic complete-data gradient. Analytic score components also cover
Normal, lognormal, Weibull, Laplace, population-location, and residual-error
parameters. Gamma scale is analytic; only the Gamma shape CDF derivative uses
the shrinking centered-difference schedule. The analytic complete score agrees
with an independent global numerical derivative to 1.7e-10 for Gamma and
Weibull (`test-score-sa-analytic-gradient.R`). The 500-iteration nonlinear PK
runtime dropped to 8.7 seconds versus 2.6 seconds for common-Q, while its Gamma
margin was 1.639/38.88 versus 1.640/38.88 and likelihood differences remained
inside importance-sampling uncertainty.

The score dispatcher now also supports observed or missing binary/ordinal
conditioning margins. Eta remains the only missing state; the exact rectangle
term in `p(y,c_obs,eta)` is differentiated directly. This avoids the changing-
support issue that would arise from treating a threshold latent score as
complete data. The existing categorical retained-Q path remains outside
Theorem 5/8, but `populationAlgorithm="score-sa"` now supplies the Section 8.2
route. The combined observed/missing categorical end-to-end test passes with a
frozen metric, terminal averaging, and zero projection events.

---

## 29. Full score validation and robustness fixes (28 August 2026)

The first four-scenario production run exposed three implementation issues that
were fixed before the final validation:

1. The score recursion initially inherited 500 unit-size SAEM statistic gains,
   which destabilized weak lognormal-scale directions. Score gains are now
   separate and decrease from iteration one as
   `scoreGainScale * (k + scoreGainOffset)^(-scoreGainPower)`; `scoreBurn`
   controls only finite metric learning and the start of averaging.
2. The score study passed scalar `error.init=.12` for a proportional error
   model, causing the saemix constructor to substitute its default `b=1`.
   The corrected declaration is `c(0,.12)`. The score backend honors the
   declared residual start rather than silently applying a different M-step.
3. One extreme importance draw could make the vectorized mixed Gaussian density
   mark every row invalid. Margin validity is now rowwise, and individual
   non-finite importance weights are assigned zero weight. GC14 fixes this
   regression.

The Weibull observed-margin initializer also stalled at its hard-coded shape 2
when L-BFGS-B encountered a non-finite trial. A safe Nelder--Mead stage followed
by bounded polishing now recovers the native Weibull MLE (2.575 in the probe)
and has a regression test.

The theory-reference workflow is now: complete Gaussian FREM; use its parameter
values only as initialization; declare the fixed flexible model; start a fresh
score recursion with reset gain, metric, and averaging state. This is a finite
initialization, not an in-run model change.

All 40 validation datasets completed 1,500 score iterations. Metrics froze,
the runtime theory contract passed, and no projection occurred. With 5,000
shared importance base draws per score/common pair, score-minus-common log
likelihood averaged -0.0183 (SD 0.0490; range -0.139 to 0.112); the largest
absolute standardized difference was 0.573 using the conservative unpaired
combined MCSE. Maximum absolute differences were 0.0178 for V and 0.00430 for
CL. Score recovery remained unbiased at the same scale as common-Q.

Mean score/common runtime ratios were 1.37 (Normal), 1.37 (lognormal), 0.162
(Gamma), and 0.162 (Weibull). Thus score-SA is a somewhat slower theory oracle
for algebraically cheap cases and roughly six times faster for the expensive
nonlinear-margin common-Q paths. Primary artifacts are under
`study_frem_validation/out/score` and paired likelihood artifacts under
`study_frem_validation/out/score_pair`.

Manuscript scope decision: these score/common engineering comparisons and
runtime diagnostics are internal method validation, not paper results. The
article and Online Resource retain only the formal likelihood/convergence
argument needed to establish the method; do not fold the development history,
backend compliance audit, or score/common comparison figure back into them.

---

## 30. Single-estimator paper and optimized score implementation (28 August 2026)

The scientific method is now one estimator only: the Delyon Section 8.2
likelihood-score stochastic-approximation recursion coupled to the saemix
conditional MCMC step. It is not described as identical to ordinary
finite-statistic SAEM. Fisher's identity supplies the observed-likelihood mean
field. The older common-Q backend remains only as an internal validation
reference and is absent from the new paper; score-SA is the public default.

The score implementation was optimized without changing its gradient or
likelihood target:

1. Removed an unused generic pair-copula parameter layout that was rebuilt and
   compiled every score iteration.
2. Avoided duplicate Gaussian-vine materialization and repeated margin-object
   validation inside one iteration.
3. Added direct Normal/lognormal PIT, log-density, location-score, and
   scale-score algebra.
4. Score-SA now retains only the current controlled-MCMC batch; it no longer
   builds empirical-Q pools or nonlinear profile statistics that it never uses.
5. Analytical-score iterations use a finite-support validity check rather than
   evaluating a complete-data diagnostic objective every iteration. That
   objective is evaluated for numerical-score steps and the final averaged fit.
6. The ordinary Gaussian fixed-effect solve is bypassed when score-SA owns the
   complete fitted location block.

The N=200 four-dimensional mixed Gamma/Weibull microbenchmark fell from 22.79
to 2.41 seconds for 500 score updates (9.5-fold). Representative complete
1,500-iteration fits took 20.53 seconds (Normal/Normal), 23.30 seconds
(Gamma/Normal), and 23.18 seconds (Weibull/lognormal), with zero projections
and unchanged parameter/likelihood conclusions.

The replacement manuscript was written from a blank source at
`manuscript/score-sa/article.Rmd`. It contains one model, one estimator, one
Fisher/convergence proof, and full/partial score-collapse sections. Its current
PDF is `manuscript/score-sa/article.pdf`. The Jupyter pitch was reduced to 15
slides and now uses score-SA likelihood, recovery, timing, and conditional-VPC
results.

---

## 31. Independent Lavielle audit and corrected score validation (29 August 2026)

The independent primary-literature audit is
`audits/LAVIELLE-SAEM-SCORE-AUDIT.md`. It confirmed the central Delyon equation
(74)/Fisher-identity route but found two theorem-critical implementation gaps
and one score-wiring bug. All were corrected before the paper results were
accepted:

1. saemix random-walk proposal scales formerly adapted forever with constant
   gain. For score-SA they now tune only through `scoreBurn` and freeze.
2. Internal-coordinate projection, invalid-step backtracking, no-move events,
   and one-sided numerical fallbacks are now separate diagnostics. Post-freeze
   one-sided differences are rejected; the theorem requires all stabilization
   events eventually to cease.
3. The residual `etype` vector was mistakenly replicated even though saemix had
   already expanded it over chains (540 predictions versus 1,620 labels in the
   regression test). This invalidated the analytical score and triggered the
   former silent global fallback. `etype` is now passed once and response-block
   dimensions are enforced.

Therefore the earlier `out/score` and `out/score_pair` results are superseded
and must not be cited. The corrected 40 fits are in `out/score_corrected`; their
5,000-common-random-number likelihood audit is in
`out/score_pair_corrected`. Every corrected fit froze its metric and proposal
scale, and maximum post-freeze projection, backtracking, no-move, and numerical
fallback counts were zero. Corrected score-minus-reference likelihood
differences averaged -0.01868 (SD 0.05249; range -0.1481 to 0.09767); every fit
was within two combined Monte Carlo SEs and the largest standardized difference
was 0.609. Recovery remained essentially unchanged: V bias -0.50% to 1.14%,
CL bias -0.72% to 0.68%, and residual-SD bias -0.71% to 1.91%.

The new manuscript now uses Delyon's `f`/`g` augmented/observed likelihood
notation, the actual normalized internal-coordinate recursion, a fixed mixed
dominating measure, explicit controlled-Markov and critical-value assumptions,
finite proposal/metric tuning, iterate averaging, and analytical-score rather
than sufficient-statistic-collapse terminology. Its theorem-backed categorical
scope is restricted to at most one categorical coordinate until multivariate
rectangle quadrature has a summable score-error schedule.

---

## 32. Final Fort-compatible score route and referee GO (29 August 2026)

The final independent re-audit is appended to
`audits/LAVIELLE-SAEM-SCORE-AUDIT.md`. After several deliberately adversarial
rounds, its final verdict is **GO** for both the implementation and the formal
mathematics within the paper's declared Gaussian-copula FREM scope. This is a
conditional stationary-likelihood theorem, not a global-MLE theorem and not
Delyon Theorem 8.

The last theorem-critical correction changed the score gain default from
`0.62` to `0.80` and enforces `scoreGainPower` in `(0.75,1]`. With Lipschitz
kernel/score exponents and Fort drift moment `p >= 2`, this satisfies Fort et
al. H6. The paper now states H1--H6 precisely, including H4's integrated local
oscillation, an explicitly assumed coercive H5 extension and its global
critical-value condition, and the stopped-tail C-iv argument. It also treats
finite-difference error as a complete random absolutely summable perturbation,
uses the correct mixed Lebesgue/counting measure, and proves that exact
missing-covariate augmentation preserves the invariant minorized kernel.

The exact final 40 fits are in
`study_frem_validation/out/score_theory_final`; the shared 5,000-sample
likelihood comparisons are in
`study_frem_validation/out/score_pair_theory_final`. All 40 passed runtime
stabilization diagnostics; maximum post-freeze projection, backtracking, and
no-move counts were zero. Score-minus-reference likelihood differences had
mean -0.0541, SD 0.0730, range [-0.2830, 0.0536], maximum standardized
difference 1.166, and all 40 were within two combined Monte Carlo SEs.
Structural recovery remained close to truth: V bias -0.50% to 1.15%, CL bias
-0.71% to 0.69%, and residual-SD bias -0.68% to 1.89%.

The final standalone paper is `manuscript/score-sa/article.pdf`. The final
15-slide executable notebook and Reveal deck are
`slides/copula-saem-pitch.ipynb` and
`slides/copula-saem-pitch.slides.html`. Both use the `0.80` estimator and the
`score_theory_final` results. Core score, Fisher-identity, categorical,
end-to-end, and Gaussian-conditioning tests passed, and `R CMD INSTALL`
completed successfully.

---

## 33. Systematic all-relevant-SAEM literature audit (29 August 2026)

The systematic corpus, inclusion/exclusion rule, primary-source links and
full audit are in:

- `audits/saem-literature-20260829/CORPUS.md`
- `audits/saem-literature-20260829/SAEM-LITERATURE-AUDIT.md`

This round independently reread the Delyon SAEM and equation-74 score route,
Kuhn--Lavielle MCMC-SAEM, Gu--Kong score/MCMC lineage, Andrieu/Allassonnière/
Fort controlled-Markov theory, Polyak averaging, the curved-exponential-family
paper, saemix methodology, categorical SAEM, and the three main FREM papers.

It found and corrected three issues missed by the targeted referee rounds:

1. the manuscript now acknowledges standard Gaussian FREM's established
   robustness for empirical means/covariances under non-Normal error-free
   covariates and defines the proposed gain as full-density and conditional-
   tail modelling;
2. the theorem-backed score route now rejects custom margins unless they
   explicitly declare parameter-independent support with
   `support_fixed=TRUE`; all built-ins declare it;
3. the legacy general non-Gaussian `saemixCopula()` workflow now explicitly
   selects and labels the internal `common-q` research comparator, rather than
   being confused with the Gaussian-copula paper estimator.

Notation now follows Kuhn--Lavielle/saemix/FREM directly:
`phi_i = h(psi_i)`, `eta_i = phi_i - mu_i(beta)`. The all-Normal nested member
is displayed through the standard partitioned-Normal FREM coefficient and UPV
formulas. A new score-SA Gaussian nested-null test combines an exact MVN
density identity with the final shared-likelihood validation.

Final systematic-audit verdict: conditional GO within the declared fixed
Gaussian-copula, fixed-support, at-most-one-categorical-coordinate scope; no
global-MLE, automatic-family-selection, general R-vine, standard-error, or
coverage theorem is claimed.

---

## 34. Five-round implementation cleanup (29 August 2026)

The cleanup ledger is `CODE-CLEANUP.md`. The supported public workflow is now:

```r
population <- gaussianCopulaFrem(
  etaSd = c(0.25, 0.30),
  covariates = subjectData[, c("WT", "eGFR")])
fit <- saemix(model, data, control, population = population)
```

The constructor fixes the full Gaussian copula, joint likelihood, score
estimator, exact numerical policy, and theory-compatible gain. Stock saemix is
unchanged when `population` is absent or `NULL`.

Structural changes:

1. `copulaScoreBatchUpdate()` records only the current MCMC batch; the score
   route no longer enters the retained-particle pool.
2. `copulaScoreMstep()` now owns score-runtime state and diagnostics, while
   `copulaScoreSa.R` contains the pure score step.
3. Gaussian correlation/R-vine algebra moved to
   `gaussianCopulaCorrelation.R`; the duplicate load-order-dependent
   `copulaGaussianLogDensity()` definition was removed.
4. `gaussianCopulaFremApi.R` and `man/gaussianCopulaFrem.Rd` provide the small
   public interface; `copula/README.md` was rewritten around it.
5. Package collation and `.Rbuildignore` were cleaned. The source tarball now
   contains only package material.

All score, Fisher, Gaussian nested-null, conditioning/missingness, categorical,
legacy-workflow isolation, and direct-parameter nested-null tests passed.
`R CMD INSTALL --preclean --clean` passed. `R CMD check` on the cleaned tarball
completed with no errors; its two warnings/two notes are pre-existing
package-wide export/Rd issues recorded in `CODE-CLEANUP.md`.

The large common-Q/direct-parameter/adaptive-R-vine blocks remain loaded only
for historical scripts and tests. Their safe removal is the next cleanup
round; the ledger lists the required migration order rather than deleting them
blindly.
## 35. Second SAEM literature audit (2026-08-30)

A fresh audit extended the original corpus to Titterington (1984), Younes
(1989), Louis (1982), Gu and Zhu (2001), Liang (2010), Gruffaz et al. (2024),
and Nyberg and Jonsson (2025). Records are in
`audits/saem-literature-round2-20260830/`.

The likelihood/stationary-set proof remains valid conditional on its printed
assumptions. The manuscript was corrected to: credit the direct impute--score
recursion to Younes rather than Delyon; cite Louis for the conditional-score
identity; state that exact invariant or MH-corrected kernels are essential to
the present proof; limit iterate averaging to Cesaro consistency rather than
claiming Markov-noise efficiency; and distinguish the new full-density target
from FREM's established robustness to covariate-scope omission. The PDF was
rerendered successfully (13 pages; no unresolved citations) and visually
checked.

## 36. Multiple categorical coordinates (2026-08-30)

The former one-categorical-coordinate restriction is removed for the supported
score-SA Gaussian-copula FREM path. The mixed observed likelihood already used
joint Gaussian rectangle masses; the estimator now augments each categorical
coordinate by a unit-interval variable
`U=F(k-)+V{F(k)-F(k-)}`. This has parameter-independent support and integrates
exactly to the rectangle mass, so no `pmvnorm` derivative or quadrature error
enters the score recursion.

After every eta move, observed categorical latent scores, missing continuous
covariates, and missing categorical coordinates are sampled jointly from their
conditional Gaussian law. With multiple categorical coordinates, rectangle-
dependent random-walk eta kernels are disabled and the exact conditional-prior
independence kernel is used. Rejection caps fail closed. The theorem's Fisher,
minorization, and summable centered-difference arguments therefore apply to the
fixed augmented state.

Passing verification:

- CAT9--CAT12 in `test-gaussian-copula-frem-categorical.R`;
- `test-gaussian-copula-frem-multicategorical-end-to-end.R`;
- `test-score-sa-multicategorical-fisher-oracle.R` (maximum score discrepancy
  0.0062, maximum Monte Carlo SE 0.0113);
- original categorical, continuous Fisher, score end-to-end, and Gaussian
  nested-null regressions.

The derivation and scope certificate are in
`audits/MULTICATEGORICAL-SCORE-CERTIFICATE.md`. Ordered categorical and nominal
with an explicit latent order are supported. Permutation-invariant unordered
nominal variables still require a multinomial latent-utility extension.

## 37. Moving-support continuous margins (2026-08-30)

The score-SA route now accepts proper continuous margins declaring
`support_fixed=FALSE`. It maps latent eta and missing covariate draws to fixed
percentiles, holds those percentiles fixed during score evaluation, reconstructs
candidate native values through the quantile, and propagates moving eta values
through a fresh structural-model prediction. This includes the pathwise
response derivative that would be lost by holding the native eta fixed.

Built-in constructors are `copulaMarginUniformCentered()`,
`copulaMarginCenteredGamma()`, and `copulaMarginCovariateUniform()`. The generic
distribution adapter can use the same route when its continuous CDF/quantile
round trip is valid. Undeclared support and moving discrete labels still fail
closed.

`test-moving-support-score-sa.R` passes CDF/quantile probes, an independent
Fisher oracle (error `2.2e-10`), a nonlinear centered-Uniform/centered-Gamma eta
fit, and a missing Uniform-covariate fit. All existing score, Fisher,
categorical, multi-categorical and Gaussian-null regressions pass afterward.
The formal argument and the necessary nonregular-boundary qualification are in
`audits/MOVING-SUPPORT-SCORE-CERTIFICATE.md`.

The stationary-set theorem covers smooth compact interior trajectories. It
does not cover endpoint MLE asymptotics such as freely estimated Uniform bounds
when observed values determine the boundary.

## 38. Short-form article and predictive VPC figure (2026-08-30)

A separate five-page article was written at
`manuscript/short-form/article-short.Rmd` and rendered to
`manuscript/short-form/article-short.pdf`. It is not a patch of the long paper:
it keeps only the model, observed likelihood, score recursion, compact Fisher
proof, conditional stationary-set statement, main recovery/likelihood table,
tail-conditional VPC figure, and limitations.

`manuscript/score-sa/make_vpc_improvement_figure.R` turns the final
`score_theory_final` predictive outputs into `vpc_improvement.png` and
`vpc_improvement_metrics.csv`. Within each of the ten datasets per scenario,
subjects in the highest observed COV1 decile were conditioned on their observed
covariates and simulated 500 times using common random numbers. RMSE across the
10th, 50th and 90th response percentiles and seven time points was:

- Gamma/Normal: 0.1003 Gaussian versus 0.0313 flexible (68.8% lower);
- lognormal/Normal: 0.1802 versus 0.0756 (58.0% lower);
- Weibull/lognormal: 0.0952 versus 0.0261 (72.6% lower).

The mean RMSE decreased from 0.1252 to 0.0443 (64.6%). The short article states
that this is tail-conditional simulation-based predictive validation, where a
wrong marginal CDF changes the Gaussian score used in `p(eta | covariates)`.

An additional exploratory double-tail stress test is generated by
`study_frem_validation/make_extreme_tail_vpc_score.R` with 1,500 replicates per
subject. It evaluates the 1st, 5th, 95th and 99th response percentiles in the
same upper-COV1 population. `make_extreme_tail_vpc_figure.R` produces
`extreme_tail_vpc.png` and its metrics. Mean extreme-tail RMSE decreased from
0.1377 to 0.0551 (60.0%); scenario reductions were 67.1%, 56.3%, and 59.2%.
This is included as Figure 2 in the short article and explicitly labelled an
exploratory stress test rather than a replacement primary analysis.

## 39. Medium-form article (2026-08-31)

An 11-page medium manuscript now sits between the six-page short article and
the 14-page full theory paper:

- source: `manuscript/medium-form/article-medium.Rmd`;
- working PDF: `manuscript/medium-form/article-medium.pdf`;
- delivery PDF: `../output/pdf/flexible-margin-gaussian-copula-frem-medium.pdf`.

The main text retains the complete mathematical core: population density,
exact Gaussian-FREM nesting, categorical and moving-support fixed-reference
maps, observed marginal likelihood, score recursion, Fisher identity,
likelihood Lyapunov equation, and a grouped conditional convergence theorem.
Long arguments are moved to appendices:

- Appendix A: continuous normalization, multiple categorical coordinates,
  moving support, and the full Fisher calculation;
- Appendix B: exact conditional-proposal minorization, augmentation invariance,
  Markov-noise decomposition, and convergence mapping;
- Appendix C: almost-sure summability of centered numerical-score error and the
  limited iterate-averaging claim;
- Appendix D: the exploratory extreme-tail VPC.

The primary 10th/50th/90th tail-conditional VPC remains in the main results.
The PDF renders without unresolved citations and all 11 pages passed visual
inspection.

## 40. SAEM/FREM notation and layout alignment (2026-08-31)

The medium article was revised after a direct comparison with Delyon et al.,
Kuhn--Lavielle, Comets--Lavenu--Lavielle, Yngman et al., and Nyberg et al. It is
now 12 pages because the canonical standard-FREM construction is explicit.

Main notation changes:

- saemix NLMEM: `psi_i`, `phi_i = h(psi_i) = C_i mu + eta_par,i`, and residual
  parameters `sigma`;
- lowercase `c_i` is reserved for the observed FREM covariate vector so it does
  not collide with the saemix design matrix `C_i`;
- standard FREM: `Omega_FREM`, `Omega_par`, `Omega_cov`,
  `Omega_par,cov`, `B`, TPV, UPV, and `Omega'_par`;
- Delyon incomplete-data theory: complete likelihood `f_i`, observed likelihood
  `g_i`, log likelihood `ell`, augmented/missing state `Z_i`, and gain
  `gamma_k`;
- Gaussian copula scores are `xi_j`, avoiding collision with Delyon's `Z`;
- moving-support quantiles are written `F_j^{-1}`, leaving `Q` reserved for the
  EM auxiliary function.

The outer structure now follows the applied FREM papers: numbered Introduction,
Methods, Results, Discussion, Conclusion, and Declarations, followed by lettered
technical appendices. Simulation/VPC design is separated from Results. The
primary VPC now overlays empirical median and 10th/90th percentile marks on the
generating and fitted predictive trajectories. Between-dataset RMSE is reported
alongside mean bias.

The delivery PDF remains
`../output/pdf/flexible-margin-gaussian-copula-frem-medium.pdf`; it has no
unresolved citations and passed final 12-page visual inspection.

## 41. Second SAEM/FREM correction round (2026-08-31)

The remaining findings from the adversarial follow-up audit were implemented in
the medium article:

- all-Normal nesting is now stated as exact for the error-free Gaussian FREM
  target and as the zero-covariate-residual limit of the conventional
  NONMEM/PsN implementation;
- `Omega_par,cov` is explicitly `p x q` and
  `Omega_cov,par = t(Omega_par,cov)`, avoiding ambiguity from differing printed
  FREM subscript conventions;
- FREM response covariates `c_i` are explicitly excluded from the saemix design
  matrix `C_i` unless a separate conventional fixed-effect relation is intended;
- observed likelihood, likelihood-score estimation, convergence, analytical
  evaluation, and simulation design are all subsections of numbered Methods;
- the implemented update is identified adjacent to its equation as the
  internal-coordinate, projected, preconditioned, multiple-chain extension of
  Delyon equation (74);
- likelihood gains are labelled shared-importance Monte Carlo estimates, not a
  standalone non-nested selection test;
- predictive RMSE is explicitly against generating quantiles.

`make_tail_conditional_vpc_intervals.R` now generates 95% simulation intervals
for the 10th, 50th, and 90th percentiles using 500 simulated datasets over the
pooled upper-COV1 population. The primary VPC overlays these bands and empirical
percentile marks. The final medium PDF remains 12 pages, has no unresolved
citations, and passed visual inspection.
