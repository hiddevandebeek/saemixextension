source("copula/tests/helper-load.R")

set.seed(280842L)
n <- 90L; d <- 3L; dEta <- 1L
R <- matrix(c(1, .4, -.2, .4, 1, .25, -.2, .25, 1), d, d)
structure <- rvinecopulib::cvine_structure(c(2, 1, 3))
vine <- copulaGaussianRvineFromCor(R, structure)
X <- matrix(1, n, 1L); locMap <- matrix(c(1, 0, 0), 1L, d)
beta <- .12; w <- rep(1 / n, n)

for (family in c("gamma", "weibull", "lognormal")) {
  margins <- list(copulaMarginNormal(.4),
    switch(family,
      gamma = copulaMarginCovariateGamma(1.5, 45),
      weibull = copulaMarginCovariateWeibull(.9, 80),
      lognormal = copulaMarginCovariateLognormal(log(75), .3)),
    copulaMarginCovariateNormal(90, 16))
  E <- copulaMarginsQuantile(rvinecopulib::rvinecop(n, vine), margins)
  E[, 1L] <- E[, 1L] + beta
  f <- rep(seq(.8, 2.2, length.out = n), each = 2L)
  y <- f * (1 + .1 * rnorm(length(f)))
  response <- list(y = y, f = f, etype = rep(1L, length(f)),
    pres = c(0, .12), free = 2L)
  analytic <- copulaGaussianFremPopulationScoreStep(E, w, margins, vine,
    d, dEta, .4, X, locMap, beta, 1L,
    scoreScale = .01, finiteDifference = 2e-5,
    response = response, analyticScore = TRUE)
  numeric <- copulaGaussianFremPopulationScoreStep(E, w, margins, vine,
    d, dEta, .4, X, locMap, beta, 1L,
    scoreScale = .01, finiteDifference = 2e-5,
    response = response, analyticScore = FALSE)
  error <- max(abs(analytic$score - numeric$score))
  if (error >= 3e-4) print(rbind(analytic = analytic$score,
    numeric = numeric$score, difference = analytic$score - numeric$score))
  cat(family, "analytic score max error", signif(error, 4), "\n")
  stopifnot(error < 3e-4,
    identical(analytic$scoreMethod,
      "analytic-with-declared-local-numerical-components"))
}
cat("score-SA analytic gradient checks passed\n")
