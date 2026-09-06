## Generalized-Pareto tail index of the importance weights. It estimates how
## many moments the weight distribution has, which effective sample size cannot
## do: ESS reports the concentration of the weights that were drawn, k reports
## the tail that was not. Validated against synthetic draws of known shape,
## because a silently wrong tail estimate is worse than none.

suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet=TRUE)})

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-52s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}
set.seed(6501L)

## known shapes must be recovered
for (k in c(.1, .3, .5, .8, 1.5)) {
  n <- 40000L
  x <- sort((1 / k) * (runif(n)^(-k) - 1))       # GPD(shape k, scale 1)
  tail <- x[(n - 3000L):n] - x[n - 3000L]
  est <- saemix:::copulaGpdShape(tail)
  ok(sprintf("GPD shape %.1f recovered", k), abs(est - k) < .15,
    sprintf("k-hat %.3f", est))
}

## degenerate inputs must not produce a spurious number
ok("constant weights give shape zero",
  identical(saemix:::copulaParetoK(matrix(0, 3L, 500L)), rep(0, 3L)))
ok("too few draws give NA",
  all(is.na(saemix:::copulaParetoK(matrix(rnorm(30L), 3L, 10L)))))

## light and heavy tails are separated
light <- matrix(log(runif(20L * 4000L)), 20L)             # bounded weights
heavy <- matrix(log(runif(20L * 4000L)^(-1.2)), 20L)      # Pareto-ish tail
kLight <- saemix:::copulaParetoK(light)
kHeavy <- saemix:::copulaParetoK(heavy)
ok("bounded weights stay below the 0.7 gate",
  max(kLight, na.rm = TRUE) < .7, sprintf("max %.2f", max(kLight, na.rm = TRUE)))
ok("heavy tail exceeds the 0.7 gate",
  median(kHeavy, na.rm = TRUE) > .7,
  sprintf("median %.2f", median(kHeavy, na.rm = TRUE)))

## the metric surfaces it without disturbing the existing quantities
logw <- matrix(rnorm(40L * 2000L, 0, .4), 40L)
m <- saemix:::copulaPosteriorBridgeMetricFromLogWeight(logw)
ok("metric reports a finite maximum k",
  is.finite(m$maximumParetoK), sprintf("%.3f", m$maximumParetoK))
ok("metric reports a finite median k", is.finite(m$medianParetoK))
ok("delta log-likelihood is unchanged in form",
  is.finite(m$deltaLogLik) && is.finite(m$mcse) &&
    is.finite(m$minimumEssFraction))
ok("maximum k is at least the median", m$maximumParetoK >= m$medianParetoK)

if (nFail) stop(sprintf("Pareto-k checks failed: %d", nFail))
cat("Pareto-k diagnostic checks passed\n")
