## The centred margin evaluator remembers a few results addressed by their
## content, because a FREM fit evaluates the same conditioning column against
## the same fixed covariate margin on every objective and gradient call of the
## score step. Being a cache, it has to be invisible: the same call must give
## the same answer, a changed column must not be served a stale one, and
## neither must a changed margin on an unchanged column.

suppressPackageStartupMessages({
  library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)
})

nFail <- 0L
ok <- function(label, pass) {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-52s %s\n", label, if (isTRUE(pass)) "PASS" else "**FAIL**"))
}
fresh <- function(x, margins) {
  .copulaMarginColumnCache$entries <- list()
  copulaGaussianFremEvaluateMargins(x, margins)
}

set.seed(4)
n <- 40L
x <- cbind(stats::rnorm(n, 0, .3), stats::rgamma(n, 5, 1))
margins <- list(copulaMarginNormal(.3),
  copulaFitCovariateMargin(x[, 2L], "gamma"))

cold <- fresh(x, margins)
ok("a repeated call gives the same answer",
  identical(copulaGaussianFremEvaluateMargins(x, margins), cold))

moved <- x; moved[, 1L] <- moved[, 1L] + .1
ok("a moved column is not served a stale answer",
  identical(copulaGaussianFremEvaluateMargins(moved, margins),
    fresh(moved, margins)))
ok("returning to the first values still agrees",
  identical(copulaGaussianFremEvaluateMargins(x, margins), cold))

wider <- list(copulaMarginNormal(.9), margins[[2L]])
warm <- copulaGaussianFremEvaluateMargins(x, wider)
ok("a changed margin on unchanged values is respected",
  identical(warm, fresh(x, wider)) && !identical(warm$z[, 1L], cold$z[, 1L]))

## The augmented evaluator re-enters the same function one column at a time,
## where every column is position one. Keying the cache on position rather than
## content made those calls collide, which is why it is content-addressed.
single <- copulaGaussianFremEvaluateMargins(x[, 2L, drop = FALSE],
  margins[2L])
ok("a single-column call agrees with the two-column one",
  isTRUE(all.equal(single$z[, 1L], cold$z[, 2L])) &&
    isTRUE(all.equal(single$logMargin[, 1L], cold$logMargin[, 2L])))
ok("and the two-column call is unchanged afterwards",
  identical(copulaGaussianFremEvaluateMargins(x, margins), cold))

cat(if (nFail) sprintf("\n%d FAILURES\n", nFail) else
  "\nall margin column cache checks passed\n")
if (nFail) quit(status = 1L)
