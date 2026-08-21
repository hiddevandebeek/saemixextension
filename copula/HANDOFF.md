# Handoff: R-vine copula as the random-effects distribution in SAEM

**Goal:** a robust SAEM implementation with an R-vine copula on the etas.

**Status:** the estimator works and is a genuine MLE. What remains is robustness
engineering and validation breadth, not a question of whether the idea is sound.

---

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

**Work small and fast.** Default to N=40-60 subjects, M=500 IS draws, 2 reps, and
scale up only when an answer is ambiguous. The effects that matter here are not
subtle; runs sized for statistical comfort just block the loop.

---

## 2. The model

```
eta_j ~ N(0, sd_j^2)                      marginals stay Gaussian
(u_1..u_d) ~ R-vine copula C(.; theta_c)  dependence is the only new part
u_j = pnorm(eta_j / sd_j)
```

**Why Gaussian margins.** Non-Gaussian eta *marginals* are already free in
saemix (`transform.par`) and Monolix (`h(psi)` links) -- they are monotone maps
of a Gaussian and leave the curved exponential family (M1) intact. A Gaussian
copula with arbitrary margins is the same trick. The **only** new content in a
vine is non-Gaussian *dependence*, and that is exactly what breaks (M1).

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

**Stock saemix behaviour is unchanged unless `copulaSet()` has been called.**
That gating is what makes the nested-null test meaningful: the same
`estep`/`mstep` code runs under two priors that are provably the same
distribution.

### API

```r
copulaSet(vine, sd,
          familySet  = NULL,        # families the M-step may select from
          mode       = "joint",     # "joint" | "pool" | "sa"
          truncLvl   = Inf,         # 1L = truncate after tree 1
          refitEvery = 1L, fitFrom = 1L,
          poolMax    = 40L,
          freezeSd   = FALSE, freezeVine = FALSE,
          sdFromSS   = TRUE, verbose = FALSE)
copulaGet()      # snapshot: $vine, $sd, $delta, $famFixed, $trace, ...
copulaClear()    # revert to stock saemix
```

### The three M-step modes

- **`joint`** -- USE THIS. Maximises Q over `(mu-shift, sd, pair-copula taus)`
  together. Families and structure frozen after the first selection, so Q is a
  fixed function of a fixed-length parameter vector. This is the MLE.
- **`pool`** -- IFM: sd from the second moment, then the copula on
  `pnorm(eta/sd)`. Q-hat carried as a weighted particle pool (Delyon eq. 6 with a
  weighted empirical measure). Consistent, not efficient.
- **`sa`** -- no pool: fit the vine to the current draws, then a Robbins-Monro
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
- **`mode="joint"` fixes it** -- see §5 numbers below.
- **1-truncation is established practice**, not a novelty. Nikoloulopoulos's vine
  copula mixed models (the one published line putting a vine on *latent* random
  effects) fix 1-truncation a priori from d>=4, at N=10-40. Licence: Joe, Li &
  Nikoloulopoulos (2010, JMVA) -- a vine has (tail) dependence in **all**
  bivariate margins provided the **tree-1** pair copulas do.
- **Parsimony argument.** A full vine at d=4 has 6 dependence parameters --
  *exactly* as many as an unstructured MVN correlation matrix. Shape flexibility,
  zero parsimony. 1-truncation halves it.

### Empirical (all numbers reproducible from `out/*.rds`)

**Joint M-step is an MLE** (d=2, N=40, 2 reps, `out/llcheck_{ifm,joint}_N40.rds`):

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
`tau` was already stationary under IFM because the copula block always solved a
genuine score equation.

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
2. **Default `truncLvl = 1L`** for d >= 4 (currently `Inf`). Matches the
   literature, beats stock, removes the collapse.
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

An R-vine copula on the etas works. With `mode="joint"` it is a genuine MLE of
the copula model (all coordinates stationary for the exact marginal likelihood,
residual displacement 0.02-0.16%). Under a Gaussian truth it matches stock
saemix; under a non-Gaussian truth it corrects joint tail probabilities the
Gaussian model gets wrong by factors of 2-3, in every replicate. It has one
failure mode -- a weakly identified eta combined with a full vine -- which is
intermittent rather than systematic and is removed by truncating the vine after
tree 1, which also *beats* stock saemix on parameter error. The remaining work is
guarding that failure mode automatically, testing with covariates, and providing
standard errors.
