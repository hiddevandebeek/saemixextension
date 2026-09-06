## Executable proof obligations accompanying
## copula/FORMAL-GAUSSIAN-COPULA-SAEM-PROOF.md.
suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})

nfail <- 0L
ok <- function(label, pass, detail = "") {
  if (!isTRUE(pass)) nfail <<- nfail + 1L
  cat(sprintf("%-76s %s %s\n", label,
    if (isTRUE(pass)) "PASS" else "**FAIL**", detail))
}

## P1: EM/KL identity for every built-in continuous shape and one custom law.
y <- .7; sigma <- .35
families <- list(
  normal = list(lower = -Inf, upper = Inf, p0 = c(0, .8), p1 = c(.15, 1.05),
    ld = function(z, p) dnorm(z, p[1], p[2], log = TRUE)),
  student = list(lower = -Inf, upper = Inf, p0 = c(.7, 5), p1 = c(.9, 9),
    ld = function(z, p) dt(z / p[1], p[2], log = TRUE) - log(p[1])),
  laplace = list(lower = -Inf, upper = Inf, p0 = c(0, .55), p1 = c(.1, .8),
    ld = function(z, p) -log(2 * p[2]) - abs(z - p[1]) / p[2]),
  lognormal = list(lower = 0, upper = Inf,
    p0 = c(-.1, .55), p1 = c(.12, .72),
    ld = function(z, p) dlnorm(z, p[1], p[2], log = TRUE)),
  gamma = list(lower = 0, upper = Inf, p0 = c(2.1, .45), p1 = c(3, .31),
    ld = function(z, p) dgamma(z, p[1], scale = p[2], log = TRUE)),
  weibull = list(lower = 0, upper = Inf, p0 = c(1.3, .9), p1 = c(2.1, .72),
    ld = function(z, p) dweibull(z, p[1], scale = p[2], log = TRUE)),
  beta = list(lower = 0, upper = 1, p0 = c(2.2, 4.5), p1 = c(3.1, 3.8),
    ld = function(z, p) dbeta(z, p[1], p[2], log = TRUE)),
  custom_logistic = list(lower = -Inf, upper = Inf,
    p0 = c(-.1, .7), p1 = c(.2, 1.0),
    ld = function(z, p) dlogis(z, p[1], p[2], log = TRUE)))

emErrors <- vapply(families, function(case) {
  logh <- function(z, par) dnorm(y, mean = z, sd = sigma, log = TRUE) +
    case$ld(z, par)
  likelihood <- function(par) integrate(function(z) exp(logh(z, par)),
    case$lower, case$upper, rel.tol = 1e-10)$value
  L0 <- likelihood(case$p0); L1 <- likelihood(case$p1)
  posterior0 <- function(z) exp(logh(z, case$p0) - log(L0))
  Q <- function(par) integrate(function(z) {
    density <- posterior0(z); logValue <- logh(z, par)
    out <- numeric(length(z)); keep <- density > 0 & is.finite(logValue)
    out[keep] <- density[keep] * logValue[keep]; out
  }, case$lower, case$upper, rel.tol = 1e-9)$value
  KL <- integrate(function(z) {
    logp0 <- logh(z, case$p0) - log(L0)
    logp1 <- logh(z, case$p1) - log(L1)
    density <- exp(logp0)
    out <- numeric(length(z))
    keep <- density > 0 & is.finite(logp0) & is.finite(logp1)
    out[keep] <- density[keep] * (logp0[keep] - logp1[keep]); out
  }, case$lower, case$upper, rel.tol = 1e-9)$value
  (log(L1) - log(L0)) - (Q(case$p1) - Q(case$p0) + KL)
}, numeric(1))
ok("P1 every continuous built-in/custom margin satisfies the EM/KL identity",
  max(abs(emErrors)) < 2e-8,
  sprintf("max error=%.3g (%s)", max(abs(emErrors)),
    names(which.max(abs(emErrors)))))

## P2: conditional p(eta|c) and joint p(eta,c) give identical MH ratios.
R <- matrix(c(1, .62, .62, 1), 2L)
vine <- copulaGaussianRvineFromCor(R,
  rvinecopulib::cvine_structure(c(2, 1)))
margins <- list(copulaMarginNormal(.4),
  copulaMarginCovariateGamma(1.4, 50))
conditioning <- matrix(c(30, 90, 160), ncol = 1L)
eta0 <- matrix(c(-.1, .05, .2), ncol = 1L)
eta1 <- matrix(c(.03, -.08, .26), ncol = 1L)
conditionalU <- copulaGaussianFremConditionalUetaEvaluator(
  conditioning, vine, margins, 1L)
conditionalDifference <- conditionalU(eta1) - conditionalU(eta0)
jointDifference <- -copulaGaussianFremLogPrior(
  cbind(eta1, conditioning), vine, margins, 1L, "joint") +
  copulaGaussianFremLogPrior(
    cbind(eta0, conditioning), vine, margins, 1L, "joint")
cancelError <- max(abs(conditionalDifference - jointDifference))
ok("P2 conditional-prior collapse preserves every rowwise MH ratio",
  cancelError < 2e-12, sprintf("max error=%.3g", cancelError))

## P3: recursively updated moments equal the expanded weighted history.
set.seed(280828L)
p <- 3L; state <- list(mass = 0, sum = numeric(p), cross = matrix(0, p, p))
values <- matrix(numeric(), 0L, p); weights <- numeric()
for (k in seq_len(25L)) {
  batch <- matrix(rnorm(4L * p), 4L, p)
  gain <- if (k <= 4L) 1 else 1 / (k - 3L)
  state$mass <- (1 - gain) * state$mass + gain
  state$sum <- (1 - gain) * state$sum + gain * colMeans(batch)
  state$cross <- (1 - gain) * state$cross +
    gain * crossprod(batch) / nrow(batch)
  weights <- (1 - gain) * weights
  values <- rbind(values, batch)
  weights <- c(weights, rep(gain / nrow(batch), nrow(batch)))
}
momentError <- max(abs(c(state$mass - sum(weights),
  state$sum - colSums(weights * values),
  state$cross - crossprod(values, weights * values))))
ok("P3 recursive mass/first/cross moments equal the expanded SA history",
  momentError < 2e-12, sprintf("max error=%.3g", momentError))

## P4: q micro-gain updates equal one blocked update with effective Gamma.
old <- -.37; innovations <- c(.2, -1.1, .8, .4, -.2)
gains <- c(.08, .07, .06, .05, .04)
micro <- old
for (r in seq_along(gains))
  micro <- (1 - gains[r]) * micro + gains[r] * innovations[r]
effective <- 1 - prod(1 - gains)
blockWeight <- vapply(seq_along(gains), function(r)
  gains[r] * if (r == length(gains)) 1 else
    prod(1 - gains[(r + 1L):length(gains)]), numeric(1))
blockInnovation <- sum(blockWeight * innovations) / effective
blocked <- (1 - effective) * old + effective * blockInnovation
ok("P4 fixed blocked gain is algebraically identical to its micro-updates",
  abs(micro - blocked) < 2e-15,
  sprintf("error=%.3g", micro - blocked))

## P5: the one-arbitrary profile evaluates the package's literal common Q.
set.seed(280829L)
d <- 3L; dEta <- 1L; n <- 220L
R3 <- matrix(c(1, .35, .55, .35, 1, -.12, .55, -.12, 1), 3L)
structure <- rvinecopulib::cvine_structure(c(2, 1, 3))
truthVine <- copulaGaussianRvineFromCor(R3, structure)
truthMargins <- list(copulaMarginNormal(.42),
  copulaMarginCovariateGamma(1.5, 42),
  copulaMarginCovariateNormal(90, 18))
E <- copulaMarginsQuantile(
  rvinecopulib::rvinecop(n, truthVine), truthMargins)
startMargins <- list(copulaMarginNormal(.36),
  copulaMarginCovariateGamma(2.0, 34),
  copulaMarginCovariateNormal(86, 22))
startVine <- copulaGaussianRvineFromCor(diag(d), structure)
w <- rep(1 / n, n)
fit <- copulaMaximiseGaussianFrem(E, w, startMargins, startVine,
  d, dEta, maxit = 150L, diagnostics = TRUE)
literal <- sum(w * copulaGaussianFremLogPrior(
  sweep(E, 2L, fit$delta, "-"), fit$vine, fit$margins, dEta, "joint"))
ok("P5 one-arbitrary collapse equals the literal common-Q density",
  identical(fit$backend, "gaussian-copula-frem-one-arbitrary-profile") &&
    abs(fit$value - literal) < 2e-9,
  sprintf("error=%.3g", fit$value - literal))

cat(sprintf("\n%d failure(s)\n", nfail))
if (nfail) quit(status = 1L)
