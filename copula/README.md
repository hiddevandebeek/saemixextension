# R-vine copulas as the random-effect distribution in SAEM

Research fork of saemix (branch `copula-eta`). Stock saemix behaviour is
unchanged unless `copulaSet()` has been called.

## The question

SAEM converges to a stationary point of the EXACT marginal likelihood. Does that
survive replacing the Gaussian eta distribution with an R-vine copula?

Short answer from the theory: the *target* survives untouched -- it never
depended on the prior being Gaussian. What breaks is assumption (M1) of Delyon,
Lavielle & Moulines (1999), the curved exponential family, which is what makes
the algorithm finite-dimensional and the convergence proof work.

And the scope of what is actually new is narrower than it looks:

| feature | already available? |
|---|---|
| non-Gaussian eta *marginals* | yes -- `transform.par`, Monolix `h(psi)`: monotone maps of a Gaussian, (M1) intact |
| Gaussian copula + arbitrary margins | yes -- same trick, (M1) intact, M-step unchanged |
| Gaussian pair copulas on a vine | **is** MVN exactly (Bedford-Cooke) -- the nested null |
| non-Gaussian *dependence* | NEW -- and precisely what breaks (M1) |

## Layout

    R/etaCopula.R    canonical vine representation, eta-scale density/sampling
    R/simEta.R       Gaussian null + tau-matched non-Gaussian twin
    R/simData.R      2-cmt IV bolus popPK model, 4 etas
    R/diagnostic.R   route-C statistic (vine vs Gaussian-vine on the same structure)
    R/marginalLL.R   independent IS marginal likelihood under an arbitrary prior
    R/replicate.R    one route-C replicate
    tests/           ordering + nested-null suites
    run*.R           experiment drivers, results land in out/

The copula prior itself lives in the package: `../R/copulaPrior.R`, with gated
hooks in `../R/main_estep.R` and `../R/main_mstep.R`.

## Ordering contract (load-bearing)

rvinecopulib is canonical: `dvine_structure(1:d)`, `pair_copulas[[t]][[e]]` is
the copula of `(e, e+t)` given the variables strictly between them, and

    vine variable j == eta column j == varList$ind.eta[j]

VineCopula uses a TRANSPOSED array convention; it is used only as an independent
oracle and its objects are never mixed with rvinecopulib's. `tests/test-ordering.R`
checks this exactly (density-based, not Monte Carlo).

Two traps it pins down:
- the inverse-Rosenblatt map `eta = T(z)` is order-dependent even for a plain
  MVN, so a non-centered parametrisation must pin the eta order;
- a vine containing a ROTATED pair copula is not invariant to relabelling
  (rotated copulas are not exchangeable), so with rotations the eta order is
  part of the model, not a free choice.

## Running

R is not on PATH:

    export PATH="/c/Program Files/R/R-4.5.3/bin/x64:$PATH"
    cd copula
    Rscript tests/test-ordering.R      # exact, seconds
    Rscript tests/test-nested-null.R   # end-to-end, minutes
    Rscript runPower.R 60 100 gumbel   # route-C power study
    Rscript analysePower.R
