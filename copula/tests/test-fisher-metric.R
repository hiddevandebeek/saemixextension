## The Fisher metric: is it a metric, and is it the right one?
##
## The recursion can be preconditioned by an estimate of the Fisher
## information, which is Gu and Li's B_k. Their Theorem 1 makes the payoff
## precise -- if B_k converges to the information then the recursion is
## asymptotically efficient -- but it also imposes the condition that every
## B_k be positive definite, and that condition is where three earlier
## attempts died.
##
## The one that looked most principled was Louis' identity, which writes the
## observed information as the posterior mean complete-data curvature less the
## posterior variance of the complete-data score. It is exactly true and
## numerically hopeless here: measured at a converged point those two terms
## were 1421 and 1632 on the first coordinate, for a difference near 20. The
## missing-information fraction is 99%, so the answer is a one percent residue
## of the things it is made of, and the fourteen percent error of a differenced
## Monte Carlo gradient flips its sign. The smallest eigenvalue came out at
## -544 and eleven fits in twelve died.
##
## What works is the empirical covariance of the per-subject scores, centred,
##
##   B_k = (1/N) sum_i (Delta_i - Delta_bar) (Delta_i - Delta_bar)',
##
## which is positive semi-definite by construction and needs no curvature at
## all. Delta_i is subject i's posterior mean score, accumulated along the
## recursion. The distinction from E[s s' | y] -- the outer product of the
## mean, not the mean of the outer product -- is not a detail: on this problem
## those differ by a factor of eighty, and the version this code used to
## accumulate divided every step by eighty and stranded the correlation at
## 0.34 against a true 0.45.
##
## So this test asserts the two things that failed before: that the matrix
## stays a metric at every update, and that it is the information and not some
## other matrix eighty times its size.

suppressPackageStartupMessages({
  library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)
  library(rvinecopulib)
})
source("C:/package/saemix-copula/copula/experiments/combined-natural-frem-study/functions.R")

nFail <- 0L
ok <- function(label, pass, detail) {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-52s %s  %s\n", label,
    if (isTRUE(pass)) "PASS" else "**FAIL**", detail))
}

set.seed(6120000L)
truth <- combined_truth(); n <- 250L
rho <- .45; sdV <- truth$sdlogV; sdCL <- truth$sdlogV * 1.4
eta <- sweep(matrix(rnorm(n * 2L), n, 2L) %*%
  chol(matrix(c(1, rho, rho, 1), 2L, 2L)), 2L, c(sdV, sdCL), "*")
psi <- sweep(exp(eta), 2L, c(truth$V, truth$CL), "*")
dd <- data.frame(id = rep(seq_len(n), each = length(truth$times)),
  dose = truth$dose, time = rep(truth$times, n))
dd$y <- pmax(combined_pk(psi, dd$id, cbind(dd$dose, dd$time)) *
  (1 + truth$residual * rnorm(nrow(dd))), 1e-8)

vine <- copulaGaussianRvineFromCor(matrix(c(1, .15, .15, 1), 2L, 2L),
  rvinecopulib::dvine_structure(1:2))
options(saemix.fisherTrace = TRUE)
fit <- saemix(combined_model(), combined_data(dd),
  combined_control(6120060L, 800L),
  population = copulaPopulation(vine,
    margins = list(copulaNaturalMarginLognormal(.25),
      copulaNaturalMarginLognormal(.30)),
    scale = "parameter", populationAlgorithm = "score-sa",
    scoreScale = "auto", scoreBurn = 50L))
options(saemix.fisherTrace = FALSE)

state <- copulaGet(fit)
pieces <- state$scoreState$louisPieces
trace <- state$scoreState$louisTrace
ok("the metric was actually built", !is.null(pieces) && !is.null(trace),
  if (is.null(trace)) "no trace" else sprintf("%d updates", nrow(trace)))
if (is.null(pieces)) { cat("\n1 FAILURES\n"); quit(status = 1L) }

## Gu and Li's standing condition on B_k, checked at every update rather than
## only at the end, because it is early on that a curvature-based metric loses
## definiteness.
worst <- min(trace[, "metric"], na.rm = TRUE)
ok("positive definite at every update", worst > 0,
  sprintf("smallest eigenvalue seen %.4g", worst))

## It has to be the information. The empirical FIM is the independently
## validated quantity, and centring is the only difference between them, which
## vanishes as the mean score goes to zero.
ratio <- diag(pieces$metric) / diag(pieces$meanOuter)
ok("agrees with the empirical Fisher information",
  all(ratio > .9 & ratio < 1.1),
  sprintf("diagonal ratios %.3f to %.3f", min(ratio), max(ratio)))

## And it has to be the outer product of the mean score, not the mean of the
## outer product. Those are the same matrix only if the posterior score has no
## spread at all; the metric must sit with the mean-score outer product, not
## between the two.
inflation <- median(diag(pieces$outerProduct) / diag(pieces$metric))
ok("is not the posterior mean of the outer product",
  inflation > 1.05 && max(abs(ratio - 1)) < (inflation - 1) / 2,
  sprintf("E[s s'|y] is %.3f times larger on the diagonal", inflation))

## Finally the estimates themselves, since a metric that satisfies every
## structural property and still lands in the wrong place is no use.
estimate <- c(as.numeric(fit@results@fixed.effects)[1:2],
  unname(state$margins[[1L]]$parameters[1L]),
  unname(state$margins[[2L]]$parameters[1L]),
  copulaGaussianRvineCor(state$vine, 2L)[1L, 2L])
target <- c(truth$V, truth$CL, sdV, sdCL, rho)
relative <- abs(estimate - target) / abs(target)
ok("estimates are close to the truth", all(relative < .12),
  sprintf("largest relative error %.1f%% (correlation %.3f)",
    100 * max(relative), estimate[5L]))

cat(if (nFail) sprintf("\n%d FAILURES\n", nFail) else
  "\nthe metric is the information and stays positive definite\n")
if (nFail) quit(status = 1L)
