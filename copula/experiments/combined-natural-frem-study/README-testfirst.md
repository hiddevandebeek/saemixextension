# The test-first workflow, evaluated end to end

This is the 500-replicate study of the workflow the manuscript proposes:
the score test of the latent margins runs first, and a coordinate is offered
to the candidate ranking only where a direction of the test rejects.

It exists because the published 500-replicate study
(`study_margin_selection.R`, summarised into `data/combined-study`) does not
do that. There, candidate clearance families are ranked straight from the
incumbent's posterior draws and volume is held lognormal, so the study
evaluates the ranking-and-refit route rather than the proposed procedure.
Both reviews of 7 September 2026 asked for that distinction to be made, and
the manuscript stated it as a limitation.

## What is the same, and what is not

The generator, the replicate indices and every seed derive from the
replicate number exactly as before, so replicate *i* here is the same
dataset as replicate *i* there. The standard (Gaussian FREM) arm, the
covariate screening, the refit, the importance-sampling likelihoods and
every reported endpoint are computed by the same code. Only the stage that
decides which parameter margins may change differs:

| stage | published study | this study |
| --- | --- | --- |
| standard fit | yes | yes |
| score test on the standard fit | not run | both parameters, both directions, bootstrap-calibrated |
| candidates offered to the ranking | clearance only, always | only coordinates where a direction rejects |
| volume | held lognormal | eligible when the test rejects on it |
| ranking | always run | run only if some direction rejects |
| refit | yes | yes |

If no direction rejects, no ranking is run and the flexible arm keeps
lognormal parameter margins, differing from the standard arm only in the
screened covariate margin.

Both arms of this study were fitted with the current package, so the study
is internally consistent; it is not directly comparable, fit for fit, with
the published study, which predates the re-initialisation change.

## Files

- `functions_testfirst.R` — `combined_fit_pair_test_first()` and
  `combined_run_replicate_test_first()`. The replicate driver of
  `functions.R` is reused, with its fitting stage substituted, so the record
  written per replicate has the published schema plus the score test
  (`shapeTest`), its decision (`shapeTestRejects`, `rankingRun`) and the
  per-stage timings.
- `study_margin_testfirst.R` (in `copula/cluster`) — chunk driver for
  klebsiella, resumable, one file per replicate.
- `summarize_testfirst.R` — runs the published summariser's body verbatim on
  these replicates, then adds the test-first counts:
  `shape_test_table.csv`, `shape_test_decisions.csv`,
  `shape_test_rejections.csv` and extra columns in `quality_summary.csv`.

## Reproducing

```bash
# on klebsiella, from ~/saemix_copula_study
SEED_BASE=0 bash launch.sh margin_testfirst 50 10 MODE=final ITERS=3000 BOOTSTRAP=59
```

Settings: 250 subjects, eight samples, 3,000 iterations per fit, 700
screening draws with a 2,000-draw retry, 2,500-draw likelihood samples, and
a score test on 400 posterior draws (thinning 4) with 59 bootstrap
replicates of 200 draws. One replicate takes about ten minutes on one core.

```powershell
Rscript summarize_testfirst.R replicates_testfirst out_testfirst
```
