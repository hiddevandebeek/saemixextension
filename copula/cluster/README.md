# Running studies on klebsiella

Everything here assumes the dev repo `C:\package\saemix-copula`. The bundle at
`C:\package\medium-paper-self-contained` is never edited by these scripts;
results are fetched into its `output/` directory by hand.

## The machine

`beekh@klebsiella.lacdr.leidenuniv.nl`, 64 cores, 125 GB. Password auth only --
put it in `NM_PASS`, never in a file. R there is 4.6.1 against 4.5.3 locally,
which is a reproducibility difference worth recording in anything published
from it.

Remote layout, all under `~/saemix_copula_study`:

```
lib/                 the private library the fork installs into
pkg/                 unpacked source
functions.R          combined-natural-frem-study helpers
ng_functions.R       non-gaussian-eta-demonstration helpers
study_*.R            chunk drivers
summarize.R          combined study summariser
summarize_replicates.R   eta demonstration summariser
out/                 per-chunk logs
out_<study>/         per-replicate rds
```

A private library rather than the system one, because a plain `saemix` may be
installed there for other work and silently replacing it would be rude.

## Deploy

```bash
cd C:/package/saemix-copula/copula/cluster
NM_PASS=... bash deploy.sh
```

Tars the working tree, not `git archive HEAD` -- the changes worth deploying are
usually uncommitted. Excludes `copula/`, which is 1.8 GB of experiment output
that `R CMD INSTALL` never looks at; the study scripts and the experiments'
`functions.R` are copied separately. The package is pure R, so install is just
untar and `R CMD INSTALL` -- no compiled code, unlike the admixr2 harness this
was adapted from.

**If install fails on `Collate`**, a new file in `R/` is missing from
`DESCRIPTION`. `devtools::load_all` ignores `Collate` so this never shows up
locally, but the package will not install anywhere until it is fixed. Every file
in `R/` must be listed.

## Launch

```bash
# on the cluster
cd ~/saemix_copula_study
SEED_BASE=0 bash launch.sh <study> <nchunk> <reps_per_chunk> KEY=VALUE...
```

```bash
# the combined study, 500 replicates at 3000 iterations
SEED_BASE=0 bash launch.sh margin_selection 50 10 \
  MODE=final ITERS=3000 FUNCTIONS=$HOME/saemix_copula_study/functions.R

# the eta demonstration, 200 replicates
SEED_BASE=0 bash launch.sh eta_demonstration 25 8 \
  FUNCTIONS=$HOME/saemix_copula_study/ng_functions.R
```

One chunk is one core; the package is single-threaded. Keep `nchunk` under the
free core count -- the box is shared, and 50 chunks against a background load of
7 sat at 53 of 64, which is about right.

Chunks never collide: each gets its own `TAG`, its own log, and a disjoint slice
of replicate indices. Replicates seed themselves from their index, so replicate
*i* is the same simulated data as published replicate *i* -- comparisons against
the frozen results are **paired**, not independent. Every replicate writes its
own rds and a rerun skips what exists, so an interrupted run costs one replicate
per chunk.

## Watch

```bash
NM_PASS=... bash watch.sh margin_selection
```

Three traps, all of which produce a confident and wrong "everything died":

- R appears as `--file=study_x.R`, so `pgrep -f "Rscript study"` matches nothing.
- **`pkill -f PATTERN` matches the shell running it**, because the pattern
  string is in that shell's own command line. Bracket the first character --
  `pkill -f "[f]ile=study_x.R"` -- or you will kill your own session before it
  launches anything. This silently ate one launch.
- A chunk that has not finished its first replicate has written nothing, which
  is not the same as having failed.

## Summarise and fetch

The published summarisers hardcode a Windows `root`. Do not edit them; give them
a tree of the right shape and rewrite the one line:

```bash
mkdir -p cs/replicates && cp out_margin_selection/replicate_*.rds cs/replicates/
sed "s|^root <- .*|root <- \"$HOME/saemix_copula_study/cs\"|" summarize.R > cs_summarize.R
R_LIBS_USER=$HOME/saemix_copula_study/lib Rscript cs_summarize.R final
```

`summarize.R` defaults to `mode = "pilot"` with no argument and will then claim
there are no replicate files. Pass `final`. It needs `dplyr`, `tidyr`,
`patchwork`, `ggplot2`, all present on the cluster.

Outputs land in `<root>/out`. Fetch with `scp` into
`medium-paper-self-contained/output/rerun-<study>/`.

## Study settings

`combined_run_replicate` (combined study) and `ng_run_replicate` (eta
demonstration) do the work; the `study_*.R` files only assign replicates to
chunks and fix paths. Final settings are 250 subjects, 50000 VPC draws, 700
screening draws, 2500 likelihood draws. `ITERS` overrides the iteration count;
nothing else about the gain needs touching, because preheat is `total/5` so the
ratio that matters stays at 0.2 whatever the total.

Timings: about 8-9 minutes per replicate at 3000 iterations, 250 subjects. 500
replicates on 50 chunks took roughly 90 minutes.

## The parameter-scale score route

`copulaPopulation(..., scoreNaturalRoute = "psi")` is the default from
2026-09-05 and is what the combined study must be run with: the natural
parameter is the latent variable and the population score is the ordinary
complete-data score. The earlier `"reference"` route (fixed percentile,
differentiated through the quantile map and the response) is what produced
the 3000-iteration rerun fetched into the bundle, and is the cause of the
flexible arm's clearance drift there; see `copula/HANDOFF.md`. A study rerun
needs a fresh `deploy.sh` (the fix is in the working tree), the previous
`out_margin_selection` moved aside (a rerun skips replicates that exist), and
then the same launch line as before.

## Standard errors

Fits that report standard errors must run with
`copulaPopulation(..., scoreTerminal = 500)` (or a fifth of the run, at most
two fifths): the last K iterations freeze the parameter at its Polyak
average and re-average the per-subject scores with gain 1/k, which is what
makes the Delattre-Kuhn information calibrated (see `copula/HANDOFF.md`).
The combined study does not report standard errors and runs with the
default `scoreTerminal = 0`.

## Things not to change

- **`scoreGainPower` stays in (0.75, 1]**, currently 0.8. The validator enforces
  it because Fort et al.'s H6 requires it; it is theorem-critical.
- **`scoreMetricRidge` stays at 1e-3.** At 1e-6 a two-compartment model with
  four random effects and three observations fails outright, and an undersized
  ridge silently degrades posterior draws enough to invert a shape comparison
  while the fit still looks converged.
- **Preheat length matters more than it looks.** 1000 of 1500 iterations gave
  53% coverage where 200 of 1500 gave 96%. Leave it at the default so it scales
  with the run.
- `gainScale` no longer exists; the three-phase schedule replaced it. An old
  script passing it fails with "unused argument", which is how the eta study
  turned out to be unrunnable against the current package.

## Traces

`combined_run_replicate` now stores a thinned iteration trace in
`result$trace$standard` and `result$trace$flexible` -- every tenth entry plus the
last, carrying `kiter`, gain, objective, `scoreMax`, residual, typical values,
margin and copula parameters. The two arms are stored separately because their
margin parameters differ by family.

The trace holds two entries per iteration; take the last per `kiter` or the
residual draws a sawtooth that is an artefact of logging.
`output/rerun-combined-study/make_convergence_plots.R` does this.

For convergence-to-truth plots, `result$truth` is stored per replicate. Plot
**derived** margin quantities via `copulaMarginDerived` rather than raw family
parameters: the flexible arm's family is selected, so a raw `margin_shape` is
not comparable across replicates and would drop every lognormal-selected one.
