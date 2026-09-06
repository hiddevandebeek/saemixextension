## The bridge objective's gradient is analytic through the copula and the
## importance weights, with only the per-margin derivatives differenced. A wrong
## chain rule would still optimise, just to the wrong place, so it is checked
## against a difference quotient of the objective across families and points.

suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet=TRUE)})
suppressPackageStartupMessages(library(rvinecopulib))

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-54s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}

root <- "C:/package/saemix-copula/copula/experiments/combined-natural-frem-study"
source(file.path(root, "functions.R"))
simulation <- combined_simulate(1401001L, 150L)
sxdata <- combined_data(simulation$data)
std <- saemix(combined_model(), sxdata, combined_control(1401011L, 500L),
  population = combined_population(combined_initial_vine(),
    combined_standard_margins(simulation$crp), simulation$crp,
    "transformed-additive"))
s0 <- combined_fit_start(std); st0 <- copulaGet(std)
incumbent <- saemix(combined_model(s0$fixed, s0$residual, s0$omega), sxdata,
  combined_control(1401021L, 500L),
  population = combined_population(st0$vine,
    c(st0$margins[1:2], list(combined_flexible_covariate_margin(simulation$crp))),
    simulation$crp, "transformed-additive"))
state <- copulaGet(incumbent)
draws <- saemix:::copulaPosteriorEtaDraws(incumbent, 300L, 300L, 1401031L)
train <- saemix:::copulaNaturalBridgePrepare(
  saemix:::copulaNaturalPosteriorData(incumbent, draws), state)
supports <- c("positive", "positive")

worst <- 0
for (families in list(c("lognormal", "gamma"), c("lognormal", "weibull"),
                      c("gamma", "gamma"), c("lognormal", "lognormal"))) {
  st <- lapply(seq_along(families), function(j)
    saemix:::copulaNaturalMarginStart(families[j], train$natural[, j],
      train$typical[, j], supports[j]))
  layout <- saemix:::copulaMarginLayout(st)
  materialize <- function(v) saemix:::copulaMarginsWithParameters(st, layout,
    saemix:::copulaScoreFromInternal(v, layout$lower, layout$upper))
  objective <- function(v) -saemix:::copulaNaturalBridgeMetric(
    materialize(v), train, state)$deltaLogLik
  x0 <- saemix:::copulaScoreToInternal(layout$par, layout$lower, layout$upper)
  for (shift in list(0, .15, -.2)) {
    x <- x0 + shift
    analytic <- -saemix:::copulaNaturalBridgeGradient(materialize(x), train,
      state, materialize, x)
    numeric <- vapply(seq_along(x), function(k) {
      h <- 1e-4; p <- m <- x; p[k] <- p[k] + h; m[k] <- m[k] - h
      (objective(p) - objective(m)) / (2 * h) }, numeric(1))
    relative <- max(abs(analytic - numeric) / pmax(abs(numeric), 1e-6))
    worst <- max(worst, relative)
    ok(sprintf("%s at shift %+.2f", paste(families, collapse = "/"), shift),
      relative < 1e-4, sprintf("rel %.1e", relative))
  }
}
ok("worst relative error across all points", worst < 1e-4,
  sprintf("%.1e", worst))

## the shared state must not change the objective
families <- c("lognormal", "gamma")
st <- lapply(seq_along(families), function(j)
  saemix:::copulaNaturalMarginStart(families[j], train$natural[, j],
    train$typical[, j], supports[j]))
layout <- saemix:::copulaMarginLayout(st)
margins <- saemix:::copulaMarginsWithParameters(st, layout, layout$par)
direct <- saemix:::copulaNaturalBridgeCandidateLogPrior(margins, train, state)
shared <- saemix:::copulaNaturalBridgeState(margins, train)$logCandidate
ok("shared state reproduces the candidate log prior",
  identical(direct, shared))

if (nFail) stop(sprintf("bridge gradient checks failed: %d", nFail))
cat("bridge gradient checks passed\n")
