## Acceptance tests for the direct natural-parameter margin prototype.
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R")
source("R/parameterMargins.R")

fail <- 0L
ok <- function(label, condition, detail = "") {
  if (isTRUE(condition)) cat("PASS", label, detail, "\n") else {
    cat("FAIL", label, detail, "\n"); fail <<- fail + 1L
  }
}

## The standard log-Normal saemix model is exactly a direct lognormal margin.
eta_normal <- copulaMarginNormal(.30)
adapter <- parameterMarginTransformedAdditive(eta_normal, transform = 1L)
direct <- parameterMarginLognormal(.30, anchor = "median", link = "log")
psi <- exp(seq(log(4), log(22), length.out = 101L))
predictor <- rep(log(10), length(psi))
ld_adapter <- adapter$log_density(psi, predictor, adapter$parameters)
ld_direct <- direct$log_density(psi, predictor, direct$parameters)
ok("P1 transformed Normal eta equals direct lognormal density",
   max(abs(ld_adapter - ld_direct)) < 1e-12)
ok("P2 transformed Normal eta equals direct lognormal CDF",
   max(abs(adapter$cdf(psi, predictor, adapter$parameters) -
           direct$cdf(psi, predictor, direct$parameters))) < 1e-12)

## Exponentiated Student eta has no finite positive moments.
log_student <- parameterMarginTransformedAdditive(copulaMarginStudent(.3, 5), 1L)
student_moments <- log_student$moments(log(10), log_student$parameters)
ok("P3 exponentiated Student moment guard", all(is.infinite(student_moments)))

## Direct families have the declared natural-scale anchors.
gamma_mean <- parameterMarginGamma(shape = 8, anchor = "mean")
gamma_median <- parameterMarginGamma(shape = 8, anchor = "median")
weibull_median <- parameterMarginWeibull(shape = 2.5, anchor = "median")
beta_mean <- parameterMarginBeta(precision = 25)
ok("P4 direct gamma mean anchor",
   abs(gamma_mean$moments(log(10), gamma_mean$parameters)[, "mean"] - 10) < 1e-12)
ok("P5 direct gamma median anchor",
   abs(gamma_median$quantile(.5, log(10), gamma_median$parameters) - 10) < 1e-12)
ok("P6 direct Weibull median anchor",
   abs(weibull_median$quantile(.5, log(10), weibull_median$parameters) - 10) < 1e-12)
ok("P7 direct beta mean anchor",
   abs(beta_mean$moments(stats::qlogis(.35), beta_mean$parameters)[, "mean"] - .35) < 1e-12)

## Normalization on real, positive, and unit supports.
integral <- function(m, predictor, lo, hi)
  stats::integrate(function(x) exp(m$log_density(
    x, rep(predictor, length(x)), m$parameters)), lo, hi,
    rel.tol = 1e-9)$value
ok("P8 direct Student normalizes",
   abs(integral(parameterMarginStudent(.4, 6), 2, -Inf, Inf) - 1) < 1e-8)
ok("P9 direct gamma normalizes",
   abs(integral(gamma_mean, log(10), 0, Inf) - 1) < 1e-8)
ok("P10 direct beta normalizes",
   abs(integral(beta_mean, stats::qlogis(.35), 0, 1) - 1) < 1e-8)

## Joint Gaussian-copula/lognormal likelihood equals the old eta likelihood plus
## the exact log-Jacobian of psi = exp(mu + eta).
vine <- etaVineGaussian(matrix(c(1, .55, .55, 1), 2))
margins_direct <- list(direct, direct)
margins_eta <- list(eta_normal, eta_normal)
u <- withSeed(912L, rvinecopulib::rvinecop(300L, vine))
pred <- cbind(rep(log(10), 300L), rep(log(3), 300L))
psi2 <- parameterMarginsQuantile(u, pred, margins_direct)
eta <- log(psi2) - pred
new_log <- parameterLogPrior(psi2, pred, vine, margins_direct)
old_log <- copulaLogPrior(eta, vine, margins = margins_eta) - rowSums(log(psi2))
ok("P11 joint transformed-Gaussian null includes exact Jacobian",
   max(abs(new_log - old_log)) < 2e-12)

if (fail) stop(fail, " direct parameter-margin tests failed")
cat("All direct parameter-margin acceptance tests passed.\n")
