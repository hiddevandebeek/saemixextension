# Combined natural-parameter and covariate-margin experiment

This experiment replaces the manuscript's separate natural-clearance and
conditional-covariate simulations with one paired nonlinear study. Each
dataset contains:

- lognormal natural volume;
- Gamma natural clearance (shape 2.5, geometric mean 3.5);
- Gamma CRP (shape 2, scale 20);
- Gaussian-score correlations 0.20 (V--CL), 0.15 (V--CRP), and 0.60
  (CL--CRP);
- 250 subjects, eight observations per subject, and proportional residual SD
  0.12 in the final study.

The Gaussian comparator uses Normal additive parameter effects and a Normal
CRP margin. The flexible workflow screens the observed CRP margin, fits a
fixed flexible-covariate incumbent, screens the natural clearance family from
two independent full-posterior pools, and starts a fresh fixed-family score-SA
fit. No family changes within a score recursion.

The 10-dataset pilot uses smaller Monte Carlo settings but the identical model
and workflow. Run and inspect it with:

```powershell
Rscript run.R 10 8 pilot
Rscript summarize.R pilot
Rscript inspect_pilot.R
```

All prespecified pilot gates passed. The final checkpointed run is:

```powershell
Rscript run.R 500 8 final
Rscript summarize.R final
```

Each replicate is written atomically. Restarting `run.R` skips completed
replicates. Raw replicate files and worker logs are intentionally ignored by
Git; the final summary tables and publication figure are copied into the
self-contained manuscript bundle after validation.
