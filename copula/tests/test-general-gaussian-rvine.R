## Exact Gaussian covariance mapping on arbitrary full R-vine structures.
suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
  library(rvinecopulib)
})

n_fail <- 0L
ok <- function(label, pass, detail = "") {
  if (!isTRUE(pass)) n_fail <<- n_fail + 1L
  cat(sprintf("%-70s %s %s\n", label,
              if (isTRUE(pass)) "PASS" else "**FAIL**", detail))
}

dmvn0 <- function(z, R) {
  ch <- chol(R)
  q <- rowSums((z %*% solve(ch))^2)
  -0.5 * (ncol(z) * log(2 * pi) + 2 * sum(log(diag(ch))) + q)
}

set.seed(82431)
roundtrip_error <- density_error <- native_error <- numeric()
orders <- character()
for (d in 2:6) for (replicate in 1:5) {
  A <- matrix(rnorm(d * d), d, d)
  R <- cov2cor(crossprod(A) + diag(d) * 0.5)
  structure <- rvine_structure_sim(d, natural_order = FALSE)
  vine <- copulaGaussianRvineFromCor(R, structure)
  rebuilt <- copulaGaussianRvineCor(vine)
  roundtrip_error <- c(roundtrip_error, max(abs(rebuilt - R)))
  u <- matrix(runif(400 * d, 1e-6, 1 - 1e-6), ncol = d)
  z <- qnorm(u)
  reference <- dmvn0(z, R) - rowSums(dnorm(z, log = TRUE))
  exact <- copulaNestedGaussianLogDensity(z, vine,
    lapply(rep(1, d), copulaMarginNormal)) - rowSums(dnorm(z, log = TRUE))
  native <- log(dvinecop(u, vine))
  density_error <- c(density_error, max(abs(reference - exact)))
  native_error <- c(native_error, abs(reference - native))
  orders <- c(orders, paste(structure$order, collapse = "-"))
}
ok("GR1 correlation -> arbitrary R-vine -> correlation is exact",
   max(roundtrip_error) < 2e-11,
   sprintf("maxabs=%.3g", max(roundtrip_error)))
ok("GR2 general-R-vine Gaussian specialization equals the MVN density",
   max(density_error) < 2e-10,
   sprintf("maxabs=%.3g", max(density_error)))
ok("GR2b native recursion agrees away from its rare extreme-tail instability",
   unname(quantile(native_error, .99)) < 2e-8,
   sprintf("99th percentile=%.3g; exact specialization is used",
     unname(quantile(native_error, .99))))
ok("GR3 test exercised non-natural general R-vine orders",
   any(vapply(strsplit(orders, "-", fixed = TRUE), function(x)
     !identical(as.integer(x), seq_along(x)), logical(1))),
   sprintf("unique orders=%d", length(unique(orders))))

R4 <- matrix(c(1, .3, -.2, .1,
               .3, 1, .25, -.15,
               -.2, .25, 1, .4,
               .1, -.15, .4, 1), 4, 4, byrow = TRUE)
v4 <- copulaGaussianRvineFromCor(R4, cvine_structure(c(3, 1, 4, 2)))
ok("GR4 an explicit non-D-vine C-vine also round-trips",
   max(abs(copulaGaussianRvineCor(v4) - R4)) < 2e-12)
R4b <- matrix(c(1, .18, -.12, .21,
                .18, 1, .31, -.09,
                -.12, .31, 1, .27,
                .21, -.09, .27, 1), 4, 4, byrow = TRUE)
fast4 <- copulaGaussianRvineUpdateCor(v4, R4b)
reference4 <- copulaGaussianRvineFromCor(R4b, v4$structure)
u4 <- matrix(runif(2000, 1e-5, 1 - 1e-5), ncol = 4L)
set.seed(82432); drawFast <- rvinecop(500L, fast4)
set.seed(82432); drawReference <- rvinecop(500L, reference4)
ok("GR4b fast Gaussian R-vine update equals constructor and native behavior",
   max(abs(copulaGaussianRvineCor(fast4) - R4b)) < 2e-12 &&
     max(abs(log(dvinecop(u4, fast4)) - log(dvinecop(u4, reference4)))) < 2e-12 &&
     max(abs(drawFast - drawReference)) < 2e-12)

truncated <- truncate_model(v4, 1L)
bad <- try(copulaGaussianRvineCor(truncated), silent = TRUE)
ok("GR5 exact unrestricted mapping fails closed for truncated vines",
   inherits(bad, "try-error"))

cat(sprintf("\n%d failure(s)\n", n_fail))
if (n_fail > 0L) quit(status = 1L)
