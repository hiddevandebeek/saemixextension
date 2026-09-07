## Every registered natural-margin family must work as a population margin in
## an actual fit, on every support the registry declares -- not merely as a
## density in isolation. Real-support families in particular were never
## exercised end to end: the study model is all-positive, so Normal, Student
## and Laplace etas only ever appear in unit tests.
##
## One compact model per support class, each with the transform that makes that
## support natural: identity for the real line, log for the positive line,
## logit for the unit interval.

suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet=TRUE)})
suppressPackageStartupMessages(library(rvinecopulib))

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-52s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}
control <- function(seed) list(seed = seed, save = FALSE, save.graphs = FALSE,
  print = FALSE, displayProgress = FALSE, warnings = FALSE,
  nbiter.saemix = c(200L, 150L), nbiter.mcmc = c(2, 2, 2, 0),
  ll.is = FALSE, fim = FALSE, map = FALSE)

## y = level_i + slope * time, level_i real  -> identity transform
linear <- function(psi, id, x) psi[id, 1L] + psi[id, 2L] * x[, 1L]
## y = level_i * exp(-rate * time), level_i positive -> log transform
decay <- function(psi, id, x) psi[id, 1L] * exp(-psi[id, 2L] * x[, 1L])
## y = fraction_i * plateau, fraction_i in (0,1) -> logit transform
share <- function(psi, id, x) psi[id, 1L] * psi[id, 2L] * (1 - exp(-x[, 1L]))

build <- function(kind, n = 120L, seed = 4001L) {
  set.seed(seed)
  times <- c(.5, 1, 2, 4, 6, 8)
  spec <- switch(kind,
    real = list(model = linear, transform = c(0L, 0L), psi0 = c(level = 2, slope = .6),
      draw = function(n) rnorm(n, 2, .8)),
    positive = list(model = decay, transform = c(1L, 1L),
      psi0 = c(level = 10, rate = .25), draw = function(n) rlnorm(n, log(10), .3)),
    unit = list(model = share, transform = c(3L, 1L),
      psi0 = c(fraction = .4, plateau = 8),
      draw = function(n) plogis(rnorm(n, qlogis(.4), .6))))
  first <- spec$draw(n)
  second <- rep(spec$psi0[2L], n)
  psi <- cbind(first, second)
  dd <- data.frame(id = rep(seq_len(n), each = length(times)),
    time = rep(times, n))
  dd$y <- spec$model(psi, dd$id, cbind(dd$time)) + rnorm(nrow(dd), 0, .15)
  list(data = saemixData(name.data = dd, header = TRUE, name.group = "id",
      name.predictors = "time", name.response = "y", verbose = FALSE),
    model = saemixModel(model = spec$model, modeltype = "structural",
      description = kind,
      psi0 = matrix(spec$psi0, 1L, dimnames = list(NULL, names(spec$psi0))),
      transform.par = spec$transform,
      covariance.model = diag(c(1L, 1L)), omega.init = diag(c(.3^2, .2^2)),
      error.model = "constant", error.init = c(.15, 0), verbose = FALSE))
}

## a two-coordinate Gaussian copula; the family under test occupies the first
## coordinate and a fixed base family the second
vine2 <- copulaGaussianRvineFromCor(matrix(c(1, .25, .25, 1), 2L, 2L),
  rvinecopulib::dvine_structure(1:2))
partner <- list(real = copulaNaturalMarginNormal(1),
  positive = copulaNaturalMarginLognormal(.25),
  unit = copulaNaturalMarginLognormal(.25))

families <- list(
  real = list(normal = copulaNaturalMarginNormal(.8),
    student = copulaNaturalMarginStudent(.8, 6),
    laplace = copulaNaturalMarginLaplace(.8)),
  positive = list(lognormal = copulaNaturalMarginLognormal(.3),
    gamma = copulaNaturalMarginGamma(4),
    weibull = copulaNaturalMarginWeibull(3),
    generalizedgamma = copulaNaturalMarginGeneralizedGamma(.3, 1)),
  unit = list(logitnormal = copulaNaturalMarginLogitNormal(.6),
    beta = copulaNaturalMarginBeta(8)))

for (kind in names(families)) {
  spec <- build(kind)
  for (nm in names(families[[kind]])) {
    margin <- families[[kind]][[nm]]
    fit <- try(suppressWarnings(saemix(spec$model, spec$data,
      control(4100L + nchar(nm)),
      population = copulaPopulation(vine2,
        margins = list(margin, partner[[kind]]),
        scale = "parameter", populationAlgorithm = "score-sa",
        scoreScale = "auto", scoreBurn = 50L, scoreGainScale = .15,
        scoreGainPower = .8, scoreGainOffset = 30,
        scoreFiniteDifference = 1e-4, scoreProjection = 24))), silent = TRUE)
    if (inherits(fit, "try-error")) {
      ok(sprintf("%-8s %-17s fits", kind, nm), FALSE,
        sub("\n.*", "", as.character(fit)))
      next
    }
    state <- copulaGet(fit)
    fitted <- state$margins[[1L]]
    finite <- all(is.finite(fitted$parameters)) &&
      all(is.finite(as.numeric(fit["results"]["fixed.effects"])))
    ## the fitted margin must still be a proper density
    u <- c(.05, .5, .95)
    typical <- rep(switch(kind, real = 2, positive = 10, unit = .4), 3L)
    q <- try(fitted$quantile(u, typical, fitted$parameters), silent = TRUE)
    proper <- !inherits(q, "try-error") && all(is.finite(q)) &&
      max(abs(fitted$cdf(q, typical, fitted$parameters) - u)) < 1e-6
    ok(sprintf("%-8s %-17s fits and stays proper", kind, nm),
      finite && proper,
      paste(sprintf("%.3f", fitted$parameters), collapse = " "))
  }
}

if (nFail) stop(sprintf("family end-to-end checks failed: %d", nFail))
cat("all registered families fit end to end\n")
