## Unit-support natural-parameter margins: logit-normal and Beta.
## Checks the density/CDF/quantile triple, the anchor, registry wiring and
## that the screening grid can be built for "unit" declared support.

suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet=TRUE)})

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-58s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}

set.seed(4021L)
families <- copulaNaturalMarginFamilies("unit")
ok("unit support has registered families", length(families) >= 2L,
  paste(families, collapse = "/"))
ok("logitnormal registered", "logitnormal" %in% families)
ok("beta registered", "beta" %in% families)

probe <- function(margin, typical) {
  p <- margin$parameters
  u <- c(.005, .05, .25, .5, .75, .95, .995)
  tv <- rep(typical, length(u))
  q <- margin$quantile(u, tv, p)
  inside <- all(q > 0 & q < 1)
  inverse <- max(abs(margin$cdf(q, tv, p) - u))
  h <- 1e-6
  numeric <- (margin$cdf(q + h, tv, p) - margin$cdf(q - h, tv, p)) / (2 * h)
  analytic <- exp(margin$log_density(q, tv, p))
  density <- max(abs(numeric - analytic) / pmax(analytic, 1e-8))
  n <- 200000L
  draws <- margin$random(n, rep(typical, n), p)
  centre <- if (margin$anchor == "mean") mean(draws) else stats::median(draws)
  list(inside = inside, inverse = inverse, density = density,
    anchor = abs(centre - typical), range = all(draws > 0 & draws < 1))
}

for (typical in c(.15, .5, .85)) {
  for (margin in list(copulaNaturalMarginLogitNormal(.6),
                      copulaNaturalMarginBeta(12))) {
    r <- probe(margin, typical)
    tag <- sprintf("%s @ typical=%.2f", margin$name, typical)
    ok(paste(tag, "quantile in (0,1)"), r$inside)
    ok(paste(tag, "cdf inverts quantile"), r$inverse < 1e-8,
      sprintf("%.2e", r$inverse))
    ok(paste(tag, "density = d/dx cdf"), r$density < 1e-5,
      sprintf("%.2e", r$density))
    ok(paste(tag, "anchor recovers typical"), r$anchor < 5e-3,
      sprintf("%.2e", r$anchor))
    ok(paste(tag, "draws in (0,1)"), r$range)
  }
}

## registry start functions recover a known truth and rank it best
typical <- .35
truth <- copulaNaturalMarginBeta(15)
x <- truth$random(20000L, rep(typical, 20000L), truth$parameters)
tv <- rep(typical, length(x))
startBeta <- saemix:::copulaNaturalMarginStart("beta", x, tv, "unit")
startLogit <- saemix:::copulaNaturalMarginStart("logitnormal", x, tv, "unit")
ok("beta start recovers precision", abs(startBeta$parameters - 15) < 1.5,
  sprintf("%.3f", startBeta$parameters))
ok("true family has higher start loglik",
  mean(startBeta$log_density(x, tv, startBeta$parameters)) >
  mean(startLogit$log_density(x, tv, startLogit$parameters)))

## screening grid builds on unit and mixed supports
grid <- saemix:::copulaNaturalFamilyGrid(NULL, c("unit", "unit"), 81L)
ok("unit grid is the family product", nrow(grid) == length(families)^2,
  paste(nrow(grid), "rows"))
mixed <- saemix:::copulaNaturalFamilyGrid(NULL,
  c("positive", "positive", "unit"), 81L)
ok("mixed support grid builds", nrow(mixed) ==
  length(copulaNaturalMarginFamilies("positive"))^2 * length(families),
  paste(nrow(mixed), "rows"))
ok("incompatible family rejected", inherits(try(
  saemix:::copulaNaturalFamilyGrid(list("gamma"), "unit", 81L),
  silent = TRUE), "try-error"))

if (nFail) stop(sprintf("unit-support margin checks failed: %d", nFail))
cat("unit-support natural margin checks passed\n")
