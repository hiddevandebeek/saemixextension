suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})

p <- .25
rawScale <- function(par)
  par["sd"] * p * (1 - p) / sqrt(1 - 2 * p + 2 * p^2)
rawCenter <- function(par) {
  b <- rawScale(par); b * (1 - 2 * p) / (p * (1 - p))
}
rawLogDensity <- function(x, par) {
  b <- rawScale(par); u <- x / b
  log(p * (1 - p) / b) - u * (p - as.numeric(u < 0))
}
rawCdf <- function(x, par) {
  b <- rawScale(par)
  ifelse(x < 0, p * exp((1 - p) * x / b),
    1 - (1 - p) * exp(-p * x / b))
}
rawQuantile <- function(u, par) {
  b <- rawScale(par)
  ifelse(u < p, b / (1 - p) * log(u / p),
    -b / p * log((1 - u) / (1 - p)))
}
margin <- copulaMarginCenteredCustom("asymmetric-laplace", c(sd = .4),
  1e-6, 10, rawLogDensity, rawCdf, rawQuantile,
  rawCenter, function(par) par["sd"], support_fixed = TRUE,
  set_scale = function(par, value) { par["sd"] <- value; par },
  roles = "scale")
copulaMarginValidate(margin)
kink <- -rawCenter(margin$parameters)
integral_piecewise <- function(fun)
  integrate(fun, -Inf, kink, rel.tol = 1e-10, subdivisions = 500L)$value +
  integrate(fun, kink, Inf, rel.tol = 1e-10, subdivisions = 500L)$value
stopifnot(isTRUE(margin$centered), isTRUE(margin$scale_is_sd),
  identical(margin$metadata$parameter_independent_support, TRUE),
  abs(integral_piecewise(function(x)
    exp(margin$log_density(x, margin$parameters))) - 1) < 1e-8,
  abs(integral_piecewise(function(x)
    x * exp(margin$log_density(x, margin$parameters)))) < 1e-8,
  abs(integral_piecewise(function(x)
    x^2 * exp(margin$log_density(x, margin$parameters))) - .4^2) < 1e-7)

bad <- try(copulaMarginCenteredCustom("bad", c(sd = .4), 1e-6, 10,
  rawLogDensity, rawCdf, rawQuantile, rawCenter,
  function(par) par["sd"], support_fixed = NA), silent = TRUE)
stopifnot(inherits(bad, "try-error"))

## A smooth bimodal distribution verifies that the adapter is not restricted
## to one of the built-in unimodal family shapes.
mixWeight <- c(.65, .35)
mixMean <- function(par) par["scale"] * c(-.6, 1.0)
mixSd <- function(par) par["scale"] * c(.18, .28)
mixCenter <- function(par) sum(mixWeight * mixMean(par))
mixScale <- function(par) {
  mu <- mixMean(par); centre <- sum(mixWeight * mu)
  sqrt(sum(mixWeight * (mixSd(par)^2 + (mu - centre)^2)))
}
mixLogDensity <- function(x, par) {
  mu <- mixMean(par); s <- mixSd(par)
  log(mixWeight[1L] * dnorm(x, mu[1L], s[1L]) +
    mixWeight[2L] * dnorm(x, mu[2L], s[2L]))
}
mixCdf <- function(x, par) {
  mu <- mixMean(par); s <- mixSd(par)
  mixWeight[1L] * pnorm(x, mu[1L], s[1L]) +
    mixWeight[2L] * pnorm(x, mu[2L], s[2L])
}
mixQuantile <- function(u, par) {
  mu <- mixMean(par); s <- mixSd(par)
  vapply(u, function(probability) uniroot(function(x)
    mixCdf(x, par) - probability,
    lower = min(mu - 12 * s), upper = max(mu + 12 * s),
    tol = 1e-12)$root, numeric(1))
}
baseMixSd <- mixScale(c(scale = 1))
mixture <- copulaMarginCenteredCustom("user-bimodal-normal-mixture",
  c(scale = .4 / baseMixSd), 1e-6, 20,
  mixLogDensity, mixCdf, mixQuantile, mixCenter, mixScale,
  support_fixed = TRUE,
  set_scale = function(par, value) {
    par["scale"] <- value / baseMixSd; par
  }, roles = "scale")
copulaMarginValidate(mixture)
q <- mixture$quantile(c(.01, .1, .5, .9, .99), mixture$parameters)
stopifnot(max(abs(mixture$cdf(q, mixture$parameters) -
  c(.01, .1, .5, .9, .99))) < 1e-8,
  abs(mixture$scale(mixture$parameters) - .4) < 1e-12)

R <- matrix(c(1, .55, .55, 1), 2L)
state <- list(vine = copulaGaussianRvineFromCor(R,
    rvinecopulib::cvine_structure(c(2, 1))),
  margins = list(mixture, copulaMarginCovariateNormal(0, 1)),
  d = 2L, dEta = 1L, dConditioning = 1L,
  conditioningName = "COV", variableName = "eta_mix")
translated <- copulaFremToFfem(state, "COV")
analytic <- copulaFfemQuantile(translated,
  matrix(qnorm(.9), ncol = 1L, dimnames = list(NULL, "COV")),
  c(.05, .5, .95))
direct <- mixture$quantile(pnorm(.55 * qnorm(.9) +
  sqrt(1 - .55^2) * qnorm(c(.05, .5, .95))), mixture$parameters)
stopifnot(max(abs(analytic$value - direct)) < 1e-10)
cat("custom centered eta margin checks passed\n")
