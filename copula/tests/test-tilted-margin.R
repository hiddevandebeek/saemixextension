## The exponentially tilted margin: a base family multiplied by exp of an
## orthonormal Legendre expansion in its own probability integral transform.
## What has to hold is that degree 0 is exactly the base family, that every
## degree is a proper anchored density with a consistent CDF and quantile,
## that the family really is the exponential family it claims to be -- score
## L_k - E[L_k], Hessian -Var(L), both matching numerical derivatives -- and
## that the Newton fit reaches shapes no shipped family can reach.

suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet=TRUE)})

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-54s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}
set.seed(4021L)

## --- degree 0 is exactly the base family ---
x <- seq(.5, 14, length.out = 400L); tx <- rep(3.5, length(x))
base <- copulaNaturalMarginTilted(0L, .3, support = "positive")
lognormal <- copulaNaturalMarginLognormal(.3, anchor = "geometric-mean")
ok("degree 0 reproduces lognormal exactly",
  max(abs(base$log_density(x, tx, base$parameters) -
          lognormal$log_density(x, tx, lognormal$parameters))) < 1e-10)
r <- seq(-3, 3, length.out = 300L); tr <- rep(.4, length(r))
ok("degree 0 reproduces Normal exactly",
  max(abs(copulaNaturalMarginTilted(0L, .8, support = "real")$log_density(
    r, tr, c(scale = .8)) -
    copulaNaturalMarginNormal(.8)$log_density(r, tr, c(sd = .8)))) < 1e-10)

## --- proper, invertible and anchored at every degree ---
for (beta in list(.4, c(.4, -.3), c(.3, .2, -.2, .1))) {
  margin <- copulaNaturalMarginTilted(length(beta), .3, beta, "positive")
  copulaNaturalMarginValidate(margin)
  u <- c(.02, .25, .5, .75, .98); typical <- rep(3.5, 5L)
  q <- margin$quantile(u, typical, margin$parameters)
  draws <- margin$random(2e5L, rep(3.5, 2e5L), margin$parameters)
  span <- range(log(draws))
  total <- stats::integrate(function(t)
    exp(margin$log_density(exp(t), rep(3.5, length(t)), margin$parameters) + t),
    span[1L] - 4, span[2L] + 4, rel.tol = 1e-10, subdivisions = 6000L)$value
  ok(sprintf("K=%d proper, invertible, anchored", length(beta)),
    abs(total - 1) < 1e-5 &&
      max(abs(margin$cdf(q, typical, margin$parameters) - u)) < 1e-4 &&
      abs(exp(mean(log(draws))) - 3.5) < .05,
    sprintf("integral %.6f anchor %.4f", total, exp(mean(log(draws)))))
}

## --- it is the exponential family it claims to be ---
## The location is held fixed while beta moves, by shifting the typical value
## with the tilt centre, so this is the natural parameterisation.
degree <- 3L
truth <- copulaNaturalMarginTilted(degree, .3, c(.2, -.1, .1), "positive")
sample <- truth$random(4000L, rep(3.5, 4000L), truth$parameters)
tiltCentre <- function(beta) {
  grid <- copulaTiltGrid(length(beta))
  lw <- drop(grid$basis %*% beta)
  w <- grid$weight * exp(lw - max(lw))
  sum(w / sum(w) * grid$z)
}
reference <- tiltCentre(c(.2, -.1, .1))
logLik <- function(beta) {
  margin <- copulaNaturalMarginTilted(degree, .3, beta, "positive")
  sum(margin$log_density(sample,
    rep(exp(log(3.5) + .3 * (tiltCentre(beta) - reference)), length(sample)),
    margin$parameters))
}
analytic <- function(beta) {
  grid <- copulaTiltGrid(degree)
  lw <- drop(grid$basis %*% beta)
  w <- grid$weight * exp(lw - max(lw)); w <- w / sum(w)
  z <- (log(sample) - (log(3.5) - .3 * reference)) / .3
  colSums(copulaTiltBasis(z, degree)) - length(sample) * drop(w %*% grid$basis)
}
difference <- function(f, beta, h = 1e-5) vapply(seq_along(beta), function(k) {
  up <- down <- beta; up[k] <- up[k] + h; down[k] <- down[k] - h
  (f(up) - f(down)) / (2 * h) }, numeric(1))
hessian <- function(f, beta, h = 1e-3) {
  H <- outer(seq_along(beta), seq_along(beta), Vectorize(function(i, j) {
    at <- function(si, sj) { v <- beta; v[i] <- v[i] + si * h
      v[j] <- v[j] + sj * h; f(v) }
    (at(1, 1) - at(1, -1) - at(-1, 1) + at(-1, -1)) / (4 * h^2) }))
  (H + t(H)) / 2
}
for (point in list(c(0, 0, 0), c(.4, -.3, .2), c(-.5, .4, -.15))) {
  exact <- analytic(point); numerical <- difference(logLik, point)
  relative <- max(abs(exact - numerical)) / max(1, max(abs(numerical)))
  eigenvalues <- eigen(hessian(logLik, point), only.values = TRUE)$values
  ok(sprintf("score and Hessian at (%s)", paste(point, collapse = ",")),
    relative < 1e-4 && all(eigenvalues < 0),
    sprintf("score error %.1e, largest eigenvalue %.3g", relative,
      max(eigenvalues)))
}

## --- the Newton fit reaches shapes the registry cannot ---
laplaceQuantile <- function(u) ifelse(u < .5, log(2 * u),
  -log(2 * (1 - u))) / 4
grid <- seq(.002, .998, length.out = 4000L)
targets <- list(
  Gamma = list(margin = copulaNaturalMarginGamma(2.5), bound = .006),
  bimodal = list(margin = copulaNaturalMarginLognormalMixture(.35, 1.2, .22),
    bound = .012),
  Laplace = list(margin = NULL, bound = .008))
for (name in names(targets)) {
  target <- targets[[name]]$margin
  if (is.null(target)) {
    sample <- 3.5 * exp(laplaceQuantile(stats::runif(20000L)))
    q <- 3.5 * exp(laplaceQuantile(grid))
    true <- log(2) - 4 * abs(log(q) - log(3.5)) - log(q)
  } else {
    sample <- target$random(20000L, rep(3.5, 20000L), target$parameters)
    q <- target$quantile(grid, rep(3.5, length(grid)), target$parameters)
    true <- target$log_density(q, rep(3.5, length(q)), target$parameters)
  }
  fitted <- copulaNaturalMarginTiltedFit(sample, rep(3.5, length(sample)), 4L,
    "positive")
  divergence <- mean(true - fitted$margin$log_density(q, rep(3.5, length(q)),
    fitted$margin$parameters))
  ok(sprintf("degree 4 fit approximates %s", name),
    is.finite(divergence) && divergence < targets[[name]]$bound,
    sprintf("KL %.5f", divergence))
}

## --- the histogram used inside the Newton step changes nothing ---
sample <- copulaNaturalMarginGamma(2.5)$random(30000L, rep(3.5, 30000L),
  copulaNaturalMarginGamma(2.5)$parameters)
residual <- log(sample) - log(3.5)
binned <- copulaNaturalMarginTiltedNewton(residual, .3, 4L)
exact <- copulaNaturalMarginTiltedNewton(residual, .3, 4L, bins = 1e6L)
ok("histogram Newton matches the full-sample Newton",
  max(abs(binned$beta - exact$beta)) < 1e-3,
  sprintf("largest coefficient difference %.2e",
    max(abs(binned$beta - exact$beta))))

cat(if (nFail) sprintf("\n%d FAILURES\n", nFail) else "\nall tilted-margin checks passed\n")
if (nFail) quit(status = 1L)
