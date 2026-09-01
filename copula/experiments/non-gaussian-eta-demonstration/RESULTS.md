# Free natural-parameter margin study

## Model and workflow

Individual clearance is generated directly as Gamma on its natural additive
scale, with shape 2.5 and geometric mean 3.5. No Gamma or Student variable is
exponentiated. The `saemix` parameter transform is used only as a reversible
internal coordinate map; the natural density and its exact Jacobian define the
same likelihood under identity or log storage.

Ten untouched data-generating seeds contain 180 subjects with eight IV
one-compartment observations each. Every dataset follows the realistic
screen-and-refit workflow:

1. fit ordinary `saemix` with lognormal individual clearance;
2. draw full individual conditional distributions from that fit;
3. screen all registered positive clearance families (lognormal, Gamma and
   Weibull) by posterior likelihood reweighting;
4. start a fresh SAEM fit with the selected family, resetting all SA/MCMC state
   while retaining the standard estimates only as numerical starting values;
5. estimate all continuous parameters for 1,000 exploratory plus 500
   decreasing-gain iterations.

Volume remains standard in this targeted clearance study. Family screening is
support-based and independent of `transform.par`.

## Principled score metric

The historical scalar `scoreScale = 0.05` is replaced by `scoreScale = "auto"`.
During the finite learning period, the implementation estimates a regularized
diagonal empirical score-information matrix and uses its inverse as
$A_k\approx I(\theta_k)^{-1}$. Eigenvalues are bounded away from zero and
infinity. The metric and MCMC proposal tuning stop adapting before the
decreasing-gain phase; population parameters remain free throughout.

## Ten-dataset results

- The free-margin fit had higher observed log likelihood in 10/10 datasets;
  every gain exceeded twice the combined likelihood Monte Carlo error.
- Median likelihood gain was 7.16 (2.5th--97.5th percentiles 1.65--10.45).
- Clearance selection chose Gamma in 7 datasets and Weibull in 3; all selected
  volume margins remained lognormal.
- The generating probability `P(CL > 8)` is 0.100. Median fitted probability
  was 0.111 for the free model and 0.119 for lognormal clearance.
- The free fit had lower overall and extreme VPC log-RMSE in all 10 datasets.
- Median VPC log-RMSE across the 1st, 10th, 50th, 90th and 99th percentiles was
  0.219 for free margins and 0.811 for lognormal clearance.
- Median extreme-quantile VPC log-RMSE was 0.266 versus 1.266.
- Median final-fit runtime ratio was 2.02 (48.6 seconds for the selected-family
  fit); median screening time was 59.5 seconds.

## Verification

- Natural Gamma, lognormal and Weibull densities and quantiles agree under
  identity and log `saemix` coordinates after the exact Jacobian.
- Natural Gamma normalizes to one under both mappings.
- The direct natural-parameter E-step, fixed-reference score, fresh fit and
  observed-likelihood scorer pass an end-to-end PK test.
- The analytic score agrees with the global numerical score to the stated
  regression tolerance.
- All 20 final validation fits passed post-learning projection, backtracking,
  no-move and metric-freeze checks.

The earlier exponentiated centred-Gamma experiment is not used as evidence.

Run `run_replicates.R 10 5 --force` to regenerate the fits and
`summarize_replicates.R` to regenerate the summaries and figure.
