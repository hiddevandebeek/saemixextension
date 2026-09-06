## The quantities a reader reports, and whether they mean what they say.
##
## The classical parameterisation writes psi = theta * exp(eta) with eta
## normal. There theta is simultaneously the median and the geometric mean, and
## omega is simultaneously the standard deviation of log psi and, near enough,
## a coefficient of variation. Those coincidences are properties of the
## lognormal, not of the parameterisation, and a tilted margin separates them.
##
## What survives exactly is the geometric mean: the location is anchored so
## that E[log psi] = log theta for every tilt coefficient, which is what makes
## the rungs of the degree ladder comparable to one another. What does not
## survive is reading the scale parameter as an inter-individual spread.
##
## So `copulaStandardErrors` reports the median, the spread of log psi and a
## coefficient of variation as derived quantities, computed rather than read
## off a coefficient. This checks the three things that could go wrong.

suppressPackageStartupMessages({
  library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)
})

nFail <- 0L
ok <- function(label, pass, detail) {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-56s %s  %s\n", label,
    if (isTRUE(pass)) "PASS" else "**FAIL**", detail))
}

theta <- 3.5; sigma <- .3

## 1. Against closed forms, where they exist. A degree-0 tilted margin is
## exactly lognormal, so all three are known in advance.
base <- copulaNaturalMarginTilted(0L, sigma, NULL, "positive")
got <- copulaMarginDerived(base, theta, base$parameters)
exact <- c(median = theta, spread = sigma, cv = sqrt(exp(sigma^2) - 1))
ok("degree 0 matches the lognormal closed forms",
  max(abs(got - exact)) < 1e-8,
  sprintf("largest error %.2g", max(abs(got - exact))))

## 2. Against a different computation. The quadrature integrates against the
## Gaussian score; a fine quantile grid inverts the CDF instead, and the two
## share no code path beyond the margin itself.
u <- seq(1e-7, 1 - 1e-7, length.out = 200001L)
worst <- 0
for (beta in list(c(.4), c(0, .5), c(-.6, .4), c(.3, -.5))) {
  margin <- copulaNaturalMarginTilted(length(beta), sigma, beta, "positive")
  got <- copulaMarginDerived(margin, theta, margin$parameters)
  x <- margin$quantile(u, theta, margin$parameters); lx <- log(x)
  reference <- c(median = stats::median(x), spread = stats::sd(lx),
    cv = stats::sd(x) / mean(x))
  worst <- max(worst, max(abs(got - reference) / abs(reference)))
}
ok("tilted margins match a fine quantile grid", worst < 1e-3,
  sprintf("largest relative error %.2g over four tilts", worst))

## 3. The anchor. This is the claim the whole parameterisation rests on, so it
## is checked across tilts rather than argued: the geometric mean is theta for
## every coefficient vector, while the median is not.
geometric <- medians <- numeric(0)
for (beta in list(numeric(0), c(.4), c(0, .5), c(-.6, .4), c(.3, -.5))) {
  margin <- if (!length(beta))
    copulaNaturalMarginTilted(0L, sigma, NULL, "positive") else
    copulaNaturalMarginTilted(length(beta), sigma, beta, "positive")
  quadrature <- copulaNormalQuadrature(160L)
  value <- margin$inverse_score(quadrature$node, theta, margin$parameters)
  geometric <- c(geometric, exp(sum(quadrature$weight * log(value))))
  medians <- c(medians, margin$inverse_score(0, theta, margin$parameters))
}
## The identity E[log psi] = log theta is exact in the definition; the anchor
## realises it through E_beta[Z], which is computed on the 2048-node tilt grid,
## so it holds to that grid's accuracy rather than to machine precision. The
## departure is about 2e-6 on a typical value of 3.5, and it does not shrink
## when this test integrates more finely, which is how it was identified as
## the grid's error and not this test's. That is five parts in ten million,
## against a standard error on the same quantity of around two percent.
ok("theta is the geometric mean under every tilt",
  max(abs(geometric - theta)) < 1e-4,
  sprintf("largest departure %.2g, i.e. %.1e relative",
    max(abs(geometric - theta)), max(abs(geometric - theta)) / theta))
ok("theta is not the median once tilted",
  max(abs(medians - theta)) > .1,
  sprintf("median moves by up to %.1f%%",
    100 * max(abs(medians - theta)) / theta))

## 4. The scale parameter is not the spread. If it were, none of this would be
## needed, so the test asserts the gap that motivates the reporting.
ratio <- numeric(0)
for (beta in list(c(0, .5), c(.3, -.5))) {
  margin <- copulaNaturalMarginTilted(2L, sigma, beta, "positive")
  ratio <- c(ratio,
    copulaMarginDerived(margin, theta, margin$parameters)[["spread"]] / sigma)
}
ok("the scale parameter misstates the spread when tilted",
  any(abs(ratio - 1) > .2),
  sprintf("sd(log psi)/scale ranges %.2f to %.2f", min(ratio), max(ratio)))

cat(if (nFail) sprintf("\n%d FAILURES\n", nFail) else
  "\nthe reported quantities are what they claim to be\n")
if (nFail) quit(status = 1L)
