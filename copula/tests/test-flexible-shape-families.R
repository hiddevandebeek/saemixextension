## Two families added to widen the natural-margin registry beyond a single
## fixed shape: a generalized gamma that contains Gamma, Weibull and lognormal,
## and a two-component lognormal mixture that can be bimodal. Both are checked
## for exact nesting of the families they generalise, internal consistency of
## density/CDF/quantile, and that the anchor still pins the population location.

suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet=TRUE)})

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-54s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}
set.seed(7731L)
x <- seq(.4, 15, length.out = 500L); tx <- rep(3.5, length(x))
compare <- function(a, b) max(abs(a$log_density(x, tx, a$parameters) -
                                 b$log_density(x, tx, b$parameters)))

## --- generalized gamma nests Gamma and Weibull exactly ---
shape <- 2.5; sg <- 1 / sqrt(shape)
ok("generalized gamma at Q = sigma is Gamma",
  compare(copulaNaturalMarginGeneralizedGamma(sg, sg),
    copulaNaturalMarginGamma(shape)) < 1e-9)
wp <- 2.2
ok("generalized gamma at Q = 1 is Weibull",
  compare(copulaNaturalMarginGeneralizedGamma(1 / wp, 1),
    copulaNaturalMarginWeibull(wp)) < 1e-9)

## --- mixture reduces to lognormal when the components coincide ---
ok("mixture at zero separation is lognormal",
  compare(copulaNaturalMarginLognormalMixture(.3, 0, .3),
    copulaNaturalMarginLognormal(.3, anchor = "geometric-mean")) < 1e-9)

## --- consistency and anchoring for both ---
probe <- function(margin, label) {
  u <- c(.01, .1, .5, .9, .99); tv <- rep(3.5, 5L)
  q <- margin$quantile(u, tv, margin$parameters)
  inverse <- max(abs(margin$cdf(q, tv, margin$parameters) - u))
  h <- 1e-6
  numeric <- (margin$cdf(q + h, tv, margin$parameters) -
    margin$cdf(q - h, tv, margin$parameters)) / (2 * h)
  analytic <- exp(margin$log_density(q, tv, margin$parameters))
  density <- max(abs(numeric - analytic) / analytic)
  n <- 200000L
  draws <- margin$random(n, rep(3.5, n), margin$parameters)
  anchor <- abs(exp(mean(log(draws))) - 3.5)
  ok(label, inverse < 1e-7 && density < 1e-4 && anchor < .03 &&
    all(draws > 0), sprintf("inv %.0e dens %.0e anchor %.3f",
      inverse, density, anchor))
}
for (p in list(c(.3, .05), c(.63, .63), c(.45, 1), c(.8, 2)))
  probe(copulaNaturalMarginGeneralizedGamma(p[1], p[2]),
    sprintf("generalized gamma sigma=%.2f Q=%.2f", p[1], p[2]))
for (p in list(c(.3, .8, .25), c(.5, 1.6, .2), c(.1, 2.5, .35)))
  probe(copulaNaturalMarginLognormalMixture(p[1], p[2], p[3]),
    sprintf("mixture w=%.1f sep=%.1f s=%.2f", p[1], p[2], p[3]))

## --- the mixture really is bimodal when it should be ---
m <- copulaNaturalMarginLognormalMixture(.4, 1.8, .18)
gridPoints <- seq(.5, 25, length.out = 3000L)
dens <- exp(m$log_density(gridPoints, rep(3.5, length(gridPoints)), m$parameters))
turning <- sum(diff(sign(diff(dens))) < 0)
ok("well-separated mixture has two modes", turning == 2L,
  sprintf("%d modes", turning))

## --- both are registrable and recoverable by direct MLE ---
truth <- copulaNaturalMarginLognormalMixture(.35, 1.2, .22)
sample <- truth$random(4000L, rep(3.5, 4000L), truth$parameters)
nll <- function(v) {
  g <- try(copulaNaturalMarginLognormalMixture(stats::plogis(v[1]),
    exp(v[2]), exp(v[3])), silent = TRUE)
  if (inherits(g, "try-error")) return(1e10)
  -sum(g$log_density(sample, rep(3.5, length(sample)), g$parameters))
}
est <- stats::optim(c(stats::qlogis(.3), log(.8), log(.25)), nll,
  control = list(maxit = 2000L))$par
ok("mixture parameters recovered from clean draws",
  abs(stats::plogis(est[1]) - .35) < .08 && abs(exp(est[2]) - 1.2) < .2,
  sprintf("w %.2f sep %.2f", stats::plogis(est[1]), exp(est[2])))

if (nFail) stop(sprintf("flexible shape family checks failed: %d", nFail))
cat("flexible shape family checks passed\n")
