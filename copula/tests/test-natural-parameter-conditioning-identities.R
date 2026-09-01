suppressPackageStartupMessages({
  library(devtools); devtools::load_all("C:/package/saemix-copula", quiet = TRUE)
  library(rvinecopulib)
})
set.seed(902201)
n <- 500L; dEta <- 2L
R <- matrix(c(1, .25, .15, .25, 1, .55, .15, .55, 1), 3L, 3L)
vine <- copulaGaussianRvineFromCor(R, cvine_structure(c(3, 1, 2)))
covMargin <- copulaMarginCovariateGamma(2, 20)
etaMargins <- list(copulaMarginNormal(.22), copulaMarginNormal(.4))
naturalMargins <- list(copulaNaturalMarginLognormal(.22),
  copulaNaturalMarginLognormal(.4))
z <- matrix(rnorm(n * 3L), n, 3L) %*% chol(R)
eta <- sweep(z[, 1:2], 2L, c(.22, .4), "*")
crp <- matrix(covMargin$quantile(pnorm(z[, 3L]), covMargin$parameters), ncol = 1L)
predictor <- matrix(rep(log(c(20, 3.5)), each = n), nrow = n)

standard <- copulaGaussianFremLogPrior(cbind(eta, crp), vine,
  c(etaMargins, list(covMargin)), dEta, "joint")
natural <- copulaNaturalFremLogPrior(eta, crp, vine,
  c(naturalMargins, list(covMargin)), dEta, predictor, c(1L, 1L), "joint")
stopifnot(max(abs(standard - natural)) < 2e-10)

kernel <- copulaNaturalFremConditionalKernel(crp, vine,
  c(naturalMargins, list(covMargin)), dEta, predictor, c(1L, 1L))
conditional <- copulaNaturalFremLogPrior(eta, crp, vine,
  c(naturalMargins, list(covMargin)), dEta, predictor, c(1L, 1L),
  "conditional")
stopifnot(max(abs(kernel$negative(eta) + conditional)) < 2e-10)

draw <- kernel$random()
phi <- predictor + draw
psi <- copulaWorkingToNatural(phi, c(1L, 1L))
typical <- copulaWorkingToNatural(predictor, c(1L, 1L))
evaluated <- copulaNaturalMarginsEvaluate(psi, typical, naturalMargins)
zc <- qnorm(covMargin$cdf(crp[, 1L], covMargin$parameters))
expected <- outer(zc, R[1:2, 3L])
stopifnot(max(abs(colMeans(evaluated$z - expected))) < .08)
cat("natural-parameter conditioning identities passed\n")
