## Sampling/order regression tests for a genuinely non-canonical R-vine.
##
## A column permutation can leave a mathematically valid vine while attaching
## the dependence to the wrong model parameters.  These checks therefore use
## unequal margins, asymmetric/rotated pair copulas, and a non-natural R-vine
## order, then compare the package density with an independent factorization.

suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
  library(rvinecopulib)
})

n_fail <- 0L
ok <- function(label, pass, detail = "") {
  if (!isTRUE(pass)) n_fail <<- n_fail + 1L
  cat(sprintf("%-66s %s %s\n", label,
              if (isTRUE(pass)) "PASS" else "**FAIL**", detail))
}

set.seed(61)
structure <- rvine_structure_sim(4, natural_order = FALSE)
pair_copulas <- list(
  list(
    bicop_dist("gumbel", parameters = 1.8),
    bicop_dist("clayton", rotation = 90, parameters = 1.2),
    bicop_dist("frank", parameters = -2.4)
  ),
  list(
    bicop_dist("gaussian", parameters = 0.42),
    bicop_dist("joe", rotation = 180, parameters = 1.55)
  ),
  list(bicop_dist("t", parameters = matrix(c(-0.35, 7), nrow = 2)))
)
vine <- vinecop_dist(pair_copulas, structure)
margins <- lapply(c(0.17, 0.31, 0.52, 0.83), copulaMarginNormal)

set.seed(62)
u <- rvinecop(6000, vine)
eta <- copulaMarginsQuantile(u, margins)
evaluated <- copulaMarginsEvaluate(eta, margins, numericalPolicy = "exact")
roundtrip_error <- max(abs(evaluated$vine_u - u))
ok("V1 PIT/quantile round trip preserves every R-vine coordinate",
   roundtrip_error < 2e-12, sprintf("maxabs=%.2e", roundtrip_error))

reference <- log(dvinecop(u, vine)) + evaluated$log_margin
package <- copulaLogPrior(eta, vine, margins = margins,
                          numericalPolicy = "exact")
density_error <- max(abs(reference - package))
ok("V2 package prior density uses the identical non-canonical ordering",
   density_error < 2e-10, sprintf("maxabs=%.2e", density_error))

permuted <- eta[, c(2, 1, 3, 4)]
permuted_density <- copulaLogPrior(permuted, vine, margins = margins,
                                   numericalPolicy = "exact")
ok("V3 a column permutation is detected as a different fitted law",
   mean(abs(package - permuted_density)) > 0.1,
   sprintf("meanabs=%.3f", mean(abs(package - permuted_density))))

copulaClear()
copulaSet(vine, margins = margins, familySet = NULL, mode = "joint")
set.seed(63)
draw <- copulaRandEta(6000)
draw_u <- copulaMarginsEvaluate(draw, margins,
                                numericalPolicy = "exact")$vine_u
tau_reference <- cor(u, method = "kendall")
tau_draw <- cor(draw_u, method = "kendall")
tau_error <- max(abs(tau_reference - tau_draw))
ok("V4 kernel-1 prior draws reproduce the declared R-vine dependence",
   tau_error < 0.035, sprintf("max Kendall-tau diff=%.3f", tau_error))
copulaClear()

cat(sprintf("\n%d failure(s)\n", n_fail))
if (n_fail > 0L) quit(status = 1L)
