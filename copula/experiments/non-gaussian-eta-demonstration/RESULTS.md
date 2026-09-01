# Free natural-parameter margin study

## Model and workflow

Individual clearance is generated directly as Gamma on its natural additive
scale, with shape 2.5 and geometric mean 3.5. No Gamma or Student variable is
exponentiated. The `saemix` parameter transform is used only as a reversible
internal coordinate map; the natural density and its exact Jacobian define the
same likelihood under identity or log storage.

Two hundred untouched data-generating seeds contain 180 subjects with eight IV
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

If the initial ranking is unresolved, screening is repeated with 2,000 posterior
draws. If it remains unresolved, the standard family is retained.

## Two-hundred-dataset results

- The free-margin procedure had higher observed log likelihood in 163/200
  datasets; 148 gains exceeded twice the combined likelihood Monte Carlo error.
  Twenty-seven unresolved datasets conservatively retained the standard fit.
- Median likelihood gain was 4.69 (2.5th--97.5th percentiles -0.56--13.00).
- Clearance selection chose Gamma in 149 datasets, Weibull in 17 and lognormal
  in 34 (including the 27 retained-standard cases). Volume remained lognormal.
- The generating probability `P(CL > 8)` is 0.100. Median fitted probability
  was 0.106 for the free model and 0.108 for lognormal clearance.
- The free procedure had lower overall and extreme VPC log-RMSE in 172/200
  datasets; nearly all ties were conservative standard retention.
- Median VPC log-RMSE across the 1st, 10th, 50th, 90th and 99th percentiles was
  0.219 for free margins and 0.689 for lognormal clearance.
- Median extreme-quantile VPC log-RMSE was 0.293 versus 1.078.
- Median final-fit runtime ratio was 2.05; median screening time was 95.7
  seconds.

## Verification

- Natural Gamma, lognormal and Weibull densities and quantiles agree under
  identity and log `saemix` coordinates after the exact Jacobian.
- Natural Gamma normalizes to one under both mappings.
- The direct natural-parameter E-step, fixed-reference score, fresh fit and
  observed-likelihood scorer pass an end-to-end PK test.
- The analytic score agrees with the global numerical score to the stated
  regression tolerance.
- All completed validation fits passed post-learning projection, backtracking,
  no-move and metric-freeze checks.

The earlier exponentiated centred-Gamma experiment is not used as evidence.

Run `run_replicates.R 200 10 --force` to regenerate the fits and
`summarize_replicates.R 200` to regenerate the summaries and figure.
