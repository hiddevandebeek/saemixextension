source("copula/tests/helper-load.R")

set.seed(280841L)
N <- 40L; etaSd <- .4; residualSd <- .25; rho <- .55
R <- matrix(c(1, rho, rho, 1), 2L, 2L)
vine <- copulaGaussianRvineFromCor(R,
  rvinecopulib::cvine_structure(c(2, 1)))
etaMargin <- copulaMarginNormal(etaSd)
truthMargin <- copulaMarginCovariateGamma(1.4, 48)
joint <- copulaMarginsQuantile(rvinecopulib::rvinecop(N, vine),
  list(etaMargin, truthMargin))
etaTruth <- joint[, 1L]; covariate <- joint[, 2L]
y <- etaTruth + residualSd * rnorm(N)
template <- copulaMarginCovariateGamma(1.3, 50)
theta <- log(template$parameters)

makeMargin <- function(logParameters) {
  parameters <- template$parameters; parameters[] <- exp(logParameters)
  copulaMarginWithParameters(template, parameters)
}
observed <- function(logParameters) {
  margin <- makeMargin(logParameters)
  zc <- qnorm(margin$cdf(covariate, margin$parameters))
  meanResponse <- etaSd * rho * zc
  sdResponse <- sqrt(etaSd^2 * (1 - rho^2) + residualSd^2)
  mean(margin$log_density(covariate, margin$parameters) +
    dnorm(y, meanResponse, sdResponse, log = TRUE))
}
complete <- function(logParameters, eta) {
  margin <- makeMargin(logParameters)
  mean(copulaGaussianFremLogPrior(cbind(eta, covariate), vine,
    list(etaMargin, margin), 1L, "joint"))
}
gradient <- function(fun, x, ..., h = 2e-4) vapply(seq_along(x), function(j) {
  plus <- minus <- x; plus[j] <- plus[j] + h; minus[j] <- minus[j] - h
  (fun(plus, ...) - fun(minus, ...)) / (2 * h)
}, numeric(1))

observedScore <- gradient(observed, theta)
margin <- makeMargin(theta)
zc <- qnorm(margin$cdf(covariate, margin$parameters))
priorVariance <- etaSd^2 * (1 - rho^2)
posteriorVariance <- 1 / (1 / priorVariance + 1 / residualSd^2)
posteriorMean <- posteriorVariance *
  (etaSd * rho * zc / priorVariance + y / residualSd^2)
B <- 1200L
scores <- matrix(NA_real_, B, length(theta))
for (b in seq_len(B)) {
  eta <- posteriorMean + sqrt(posteriorVariance) * rnorm(N)
  scores[b, ] <- gradient(complete, theta, eta = eta)
}
fisher <- colMeans(scores); mcse <- apply(scores, 2L, sd) / sqrt(B)
error <- abs(fisher - observedScore)
stopifnot(all(error <= 4 * mcse + 5e-4))
cat(sprintf(paste0("score-SA Fisher oracle passed; max error %.3g, ",
  "max MCSE %.3g\n"), max(error), max(mcse)))
