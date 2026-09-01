# Non-Gaussian eta demonstration

## Question

Does a fitted non-Gaussian eta margin recover population tails and predictive
tails that a Gaussian eta model misses, while retaining the same structural,
residual and dependence models?

## Design

One illustrative dataset contains 180 subjects and eight IV one-compartment PK
observations per subject. Volume has a Normal eta. Clearance has a centred Gamma
eta with shape 2.5 and standard deviation 0.45:

\[
\eta_{CL}=G-E(G),\qquad G\sim\mathrm{Gamma}(2.5,\;0.45/\sqrt{2.5}).
\]

Thus `E(eta_CL)=0`; on `CL_i=TVCL exp(eta_CL,i)` the distribution is bounded
below and right-skewed. The generating correlation between the latent Normal
scores for volume and clearance is 0.35. The proportional residual standard
deviation is 0.12.

Two fixed Gaussian-copula score-SA models were fitted for 1,000 plus 500
iterations from separate fresh starts:

1. Normal margins for both etas;
2. Normal volume eta and centred-Gamma clearance eta.

This is a fixed-family demonstration. It does not claim that the current
automatic Normal/Student/Laplace screen selects Gamma.

## Results

- Observed log likelihood: -265.45 (Gaussian) versus -247.46 (centred Gamma),
  a gain of 17.99. The separate likelihood Monte Carlo standard errors were
  0.44 and 0.44.
- The generating probability `P(|eta_CL| > 1)` was 0.0342. The fitted values
  were 0.0133 for the Gaussian margin and 0.0212 for the centred-Gamma margin.
- VPC log-RMSE across the 1st, 10th, 50th, 90th and 99th percentiles was 0.939
  for the Gaussian fit and 0.499 for the centred-Gamma fit.
- Extreme-quantile VPC log-RMSE (1st and 99th percentiles) was 1.477 versus
  0.767.

## Shared score and E-step acceleration

The original moving-support implementation used global centered differences
for all seven fitted coordinates and took 97.69 seconds for this same
1,500-iteration Gamma fit. The shared implementation now:

1. retains analytic scores for fixed-support margins, Gaussian dependence and
   residual error;
2. differentiates the structural response once per referenced eta coordinate
   and propagates it to location and moving-margin parameters by the chain rule;
3. caches the Gaussian correlation and Cholesky factor once per E-step
   iteration.

The final centred-Gamma fit took 57.59 seconds, a 41.0% reduction, versus 38.72
seconds for the Gaussian fit. The optimized and original implementations
reached the same printed estimates, likelihood and predictive metrics. The
hybrid score matched the former global finite-difference oracle within
`3.35e-10`; the cached E-step prior matched the literal density within `2e-12`.

The centred-Gamma fit therefore recovers the asymmetric population shape and
substantially improves the predictive tails. It does not perfectly recover the
generating shape from one finite dataset (estimated shape 3.96 versus 2.5), so
this figure should remain an illustration until supported by repeated-dataset
results.

Run `Rscript --vanilla run.R --force` to refit. Without `--force`, the script
uses the cached fits and regenerates the analysis and figure.
