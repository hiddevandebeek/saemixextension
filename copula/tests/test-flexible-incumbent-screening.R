## Screening from a FLEXIBLE fit, not only a Normal/lognormal one.
##
## Natural-margin densities take (x, typical, par) where covariate margins take
## (x, par), so the Gaussian proposal and Gaussian baseline cannot be applied to
## a parameter-scale model. Both are now dispatched on the incumbent's own
## population scale. Without that, the posterior of a flexible fit cannot be
## sampled at all, and margin estimation is confined to one round off a
## Normal/lognormal model.
##
## The decisive check is self-consistency: an incumbent scored against its own
## margins must give a likelihood ratio of exactly one, i.e. a delta
## log-likelihood of zero. A wrong baseline would still run and quietly report
## a non-zero advantage for the incumbent's own model.

suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet=TRUE)})
suppressPackageStartupMessages(library(rvinecopulib))

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-56s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}

root <- "C:/package/saemix-copula/copula/experiments/combined-natural-frem-study"
source(file.path(root, "functions.R"))
simulation <- combined_simulate(1401001L, 180L)
sxdata <- combined_data(simulation$data)
cov <- combined_flexible_covariate_margin(simulation$crp)
standard <- saemix(combined_model(), sxdata, combined_control(1401011L, 500L),
  population = combined_population(combined_initial_vine(),
    combined_standard_margins(simulation$crp), simulation$crp,
    "transformed-additive"))
start <- combined_fit_start(standard); base <- copulaGet(standard)
scales <- vapply(base$margins[1:2], function(m)
  unname(m$parameters[["sd"]]), numeric(1))

## Scoring an incumbent against its own population law must give a ratio of
## exactly one. On the parameter scale the incumbent's margins are natural and
## can be passed straight back; on the transformed-additive scale they are
## centred margins with a different signature, so the natural equivalent of a
## centred Normal eta under a log map -- a median-anchored lognormal of the
## same spread -- is used instead.
selfConsistent <- function(fit, label, equivalent = NULL) {
  state <- copulaGet(fit)
  draws <- saemix:::copulaPosteriorEtaDraws(fit, 200L, 200L, 4242L)
  ok(paste(label, "posterior draws"), !inherits(draws, "try-error") &&
    all(dim(draws) == c(fit["data"]["N"], state$dEta, 200L)))
  data <- saemix:::copulaNaturalPosteriorData(fit, draws)
  bridge <- saemix:::copulaNaturalBridgePrepare(data, state)
  ok(paste(label, "baseline is finite"),
    all(is.finite(bridge$bridgeCache$baseline)))
  own <- if (is.null(equivalent)) state$margins[seq_len(state$dEta)] else
    equivalent
  self <- saemix:::copulaNaturalBridgeMetric(own, bridge, state)
  ok(paste(label, "scores its own law at zero"), abs(self$deltaLogLik) < 1e-6,
    sprintf("%.2e", self$deltaLogLik))
  invisible(state)
}

## a standard transformed-additive incumbent must still behave
selfConsistent(standard, "transformed-additive:",
  equivalent = lapply(scales, function(sd)
    copulaNaturalMarginLognormal(sd, anchor = "median")))

## and so must a fitted flexible one, on the parameter scale
for (spec in list(
    list(nm = "gamma margin:", m = copulaNaturalMarginGamma(3)),
    list(nm = "SNP margin:  ",
      m = copulaNaturalMarginSNP(2L, max(scales[2L], .1), c(.3, -.2),
        "positive")))) {
  flexible <- saemix(combined_model(start$fixed, start$residual, start$omega),
    sxdata, combined_control(1401061L, 500L),
    population = combined_population(base$vine,
      list(copulaNaturalMarginLognormal(scales[1L]), spec$m, cov),
      simulation$crp, "parameter"))
  ok(paste(spec$nm, "fit is on the parameter scale"),
    identical(copulaGet(flexible)$populationScale, "parameter"))
  selfConsistent(flexible, spec$nm)
  selection <- try(suppressWarnings(copulaSelectParameterMargins(flexible,
    supports = c("positive", "positive"),
    candidates = list("lognormal", c("lognormal", "gamma", "weibull")),
    posteriorDraws = 200L, max.iter = 200L, seed = 909L,
    optimizerMaxit = 100L, minimumEssFraction = 0,
    minimumPosteriorEssFraction = 0, maximumMcse = Inf)), silent = TRUE)
  ok(paste(spec$nm, "screen runs from this incumbent"),
    !inherits(selection, "try-error"),
    if (inherits(selection, "try-error")) "" else
      paste(selection$families, collapse = "/"))
}

## margin parameters must be freezable, which is what block ascent relies on
frozen <- copulaNaturalMarginSNP(2L, .3, c(.2, .1), "positive")
frozen$free[] <- FALSE
layout <- saemix:::copulaMarginLayout(list(frozen))
ok("a frozen margin contributes no free parameters",
  length(layout$par) == 0L, sprintf("%d free", length(layout$par)))

if (nFail) stop(sprintf("flexible-incumbent screening checks failed: %d", nFail))
cat("flexible-incumbent screening checks passed\n")
