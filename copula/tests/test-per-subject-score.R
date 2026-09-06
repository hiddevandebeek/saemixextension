## The Fisher information is the covariance of the per-subject score, so it
## cannot be formed from the aggregate the recursion previously kept. Splitting
## the score by latent row is only safe if the split is exact, and three of the
## four blocks are already written as weighted sums over rows so they are.
##
## The fourth is not: the correlation block goes through the scatter matrix S,
## and its angle gradient is affine in S rather than linear -- the term in
## Omega alone survives at S = 0. Both mistakes that makes possible were made
## and caught here: double-counting the off-diagonal of a symmetric bump, and
## dropping the constant. So the identity is checked on every iteration of a
## real fit rather than once on a constructed example, which is what
## `saemix.scorePerRow` turns on.

suppressPackageStartupMessages({
  library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)
  library(rvinecopulib)
})
source("C:/package/saemix-copula/copula/experiments/combined-natural-frem-study/functions.R")

options(saemix.scorePerRow = TRUE)
set.seed(1501000L)
truth <- combined_truth(); R <- combined_correlation(truth); n <- 250L
z <- matrix(rnorm(n * 3L), n, 3L) %*% chol(R)
margins <- combined_truth_margins(truth)
typical <- matrix(rep(c(truth$V, truth$CL), each = n), nrow = n)
psi <- copulaNaturalMarginsQuantile(pnorm(z[, 1:2, drop = FALSE]), typical,
  margins[1:2])
crp <- margins[[3L]]$quantile(pnorm(z[, 3L]), margins[[3L]]$parameters)
dd <- data.frame(id = rep(seq_len(n), each = length(truth$times)),
  dose = truth$dose, time = rep(truth$times, n))
dd$y <- pmax(combined_pk(psi, dd$id, cbind(dd$dose, dd$time)) *
  (1 + truth$residual * rnorm(nrow(dd))), 1e-8)
sxdata <- combined_data(dd)

standard <- saemix(combined_model(), sxdata, combined_control(1501011L, 800L),
  population = combined_population(combined_initial_vine(),
    combined_standard_margins(crp), crp, "transformed-additive"))
start <- combined_fit_start(standard); state <- copulaGet(standard)
scales <- vapply(state$margins[1:2], function(margin)
  unname(margin$parameters[["sd"]]), numeric(1))
flat <- lapply(1:2, function(j)
  copulaNaturalMarginTilted(0L, scales[j], NULL, "positive"))

## Any disagreement stops the fit from inside the score step, so reaching the
## end is the check.
fit <- try(saemix(combined_model(start$fixed, start$residual, start$omega),
  sxdata, combined_control(1501060L, 800L),
  population = combined_population(state$vine,
    c(flat, list(combined_flexible_covariate_margin(crp))), crp,
    "parameter")), silent = TRUE)

if (inherits(fit, "try-error")) {
  cat("**FAIL** ", as.character(fit), "\n")
  quit(status = 1L)
}
cat(sprintf("%-54s %s\n",
  "per-row score sums to the aggregate, every iteration", "PASS"))
cat("\nper-subject score check passed\n")
