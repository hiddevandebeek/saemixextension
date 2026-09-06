## Evidence about the marginals, for a person to read.
##
## Two routes to the same question, deliberately shown side by side, because
## each is unreliable in a way the other is not and their agreement is more
## informative than either alone.
##
## The pooled ranking fits each candidate to the pooled posterior draws of one
## coordinate. It is cheap and robust, and it is biased: the draws come from a
## fit that already assumed a margin, so they are pulled toward it, and the
## pull grows with shrinkage. Its "log-likelihood" is also over draws, not
## subjects -- tens of thousands of them -- so it must never be turned into a
## BIC or a p-value. It is a ranking, nothing more.
##
## The identity route uses the exact posterior likelihood-ratio identity,
##
##   g_m / g_0 = E_{pi_0}[ f_m / f_0 ] ,
##
## which holds for any incumbent however wrong, and penalises by BIC on the
## subject count. That quantity is a valid selection criterion. What it is not
## is reliable: the importance weights degenerate when candidates are far apart
## or when several coordinates change at once, and then its Monte Carlo error
## can exceed the differences being ranked.
##
## Measured on this implementation, neither dominates. Where the data are
## informative both agree and are right; by about twenty per cent shrinkage
## both stop discriminating. So the honest output is both, plus the
## diagnostics that say whether either can be believed, and an explicit note
## when they disagree.

copulaMarginEvidence <- function(object, supports, candidates = NULL,
                                 posteriorDraws = 800L, max.iter = NULL,
                                 seed = 935001L) {
  if (!inherits(object, "SaemixObject"))
    stop("object must be a fitted SaemixObject")
  state <- copulaGet(object)
  d <- state$dEta
  supports <- rep_len(as.character(supports), d)
  draws <- copulaPosteriorEtaDraws(object, posteriorDraws, max.iter,
    as.integer(seed))
  pool <- copulaNaturalPosteriorData(object, draws)
  subject <- rep(seq_len(pool$n), pool$samples)
  names <- as.character(object["model"]["name.modpar"])[
    as.integer(state$etaIndex)]
  if (!length(names) || anyNA(names)) names <- paste0("coord", seq_len(d))

  ## Shrinkage, because it is what decides whether either route can work.
  shrinkage <- vapply(seq_len(d), function(j) {
    inner <- if (identical(supports[j], "positive"))
      log(pool$natural[, j]) - log(pool$typical[, j]) else
      pool$natural[, j] - pool$typical[, j]
    population <- unname(state$margins[[j]]$parameters[1L])
    if (!is.finite(population) || population <= 0) return(NA_real_)
    1 - stats::sd(tapply(inner, subject, mean)) / population
  }, numeric(1))

  ## The pooled ranking, one coordinate at a time.
  rows <- list()
  pooledPick <- character(d)
  for (j in seq_len(d)) {
    families <- candidates %||% copulaNaturalMarginFamilies(supports[j])
    families <- intersect(families, copulaNaturalMarginFamilies(supports[j]))
    x <- as.numeric(pool$natural[, j]); typical <- mean(pool$typical[, j])
    ll <- vapply(families, function(f) {
      m <- try(copulaNaturalMarginStart(f, x, typical, supports[j]),
        silent = TRUE)
      if (inherits(m, "try-error")) return(NA_real_)
      value <- try(sum(m$log_density(x, typical, m$parameters)), silent = TRUE)
      if (inherits(value, "try-error") || !is.finite(value)) NA_real_ else value
    }, numeric(1))
    pooledPick[j] <- if (all(is.na(ll))) NA_character_ else
      families[which.max(ll)]
    rows[[j]] <- data.frame(coordinate = names[j], family = families,
      pooled_loglik = as.numeric(ll),
      gap_to_best = as.numeric(ll) - max(ll, na.rm = TRUE),
      row.names = NULL, stringsAsFactors = FALSE)
  }
  pooled <- do.call(rbind, rows)

  ## The identity route, which may legitimately refuse to answer.
  selection <- try(suppressWarnings(copulaSelectParameterMargins(object,
    supports, candidates = candidates, posteriorDraws = posteriorDraws,
    max.iter = max.iter, seed = as.integer(seed))), silent = TRUE)
  identityPick <- rep(NA_character_, d)
  identityTable <- NULL
  if (!inherits(selection, "try-error")) {
    identityPick <- as.character(selection$families)[seq_len(d)]
    identityTable <- selection$table
  }

  ## The gap to the runner-up, which is what says whether a choice was
  ## decisive. Two families can be near-identical on a coordinate -- a gamma of
  ## shape 16 is close to a lognormal of spread 0.25 -- and then picking either
  ## is harmless and the gap correctly says so. A large gap is a finding; a gap
  ## inside the Monte Carlo error is not.
  gaps <- vapply(seq_len(d), function(j) {
    sub <- pooled[pooled$coordinate == names[j], , drop = FALSE]
    ll <- sort(sub$pooled_loglik[is.finite(sub$pooled_loglik)],
      decreasing = TRUE)
    if (length(ll) < 2L) NA_real_ else ll[1L] - ll[2L]
  }, numeric(1))
  identityGap <- if (is.null(identityTable)) NA_real_ else {
    value <- sort(identityTable$validation_bic_advantage[
      is.finite(identityTable$validation_bic_advantage)], decreasing = TRUE)
    if (length(value) < 2L) NA_real_ else value[1L] - value[2L]
  }
  identityMcse <- if (is.null(identityTable)) NA_real_ else
    stats::median(identityTable$validation_mcse, na.rm = TRUE)

  structure(list(coordinates = names, supports = supports,
    pooledGap = stats::setNames(gaps, names),
    identityGap = identityGap, identityMcse = identityMcse,
    shrinkage = stats::setNames(shrinkage, names),
    pooled = pooled, pooledPick = stats::setNames(pooledPick, names),
    identityPick = stats::setNames(identityPick, names),
    identityTable = identityTable,
    agree = pooledPick == identityPick,
    posteriorDraws = posteriorDraws),
    class = c("saemixMarginEvidence", "list"))
}

print.saemixMarginEvidence <- function(x, ...) {
  cat("Marginal evidence, from", x$posteriorDraws, "posterior draws\n\n")
  cat(sprintf("%-10s %9s %11s %8s %13s   %s
", "coordinate", "shrinkage",
    "pooled", "gap", "identity", "agree"))
  for (j in seq_along(x$coordinates)) {
    agree <- x$agree[j]
    cat(sprintf("%-10s %8.0f%% %11s %8.0f %13s   %s
", x$coordinates[j],
      100 * x$shrinkage[j], x$pooledPick[j], x$pooledGap[j],
      x$identityPick[j] %||% "-",
      if (is.na(agree)) "?" else if (agree) "yes" else "**NO**"))
  }
  if (is.finite(x$identityGap))
    cat(sprintf("
Identity BIC gap to the runner-up %.2f, median Monte Carlo error %.2f.%s
",
      x$identityGap, x$identityMcse,
      if (x$identityGap < 3 * x$identityMcse)
        "
That gap is inside the noise: the choice is not decisive and should
not be reported as one." else ""))
  high <- which(x$shrinkage > .2)
  if (length(high))
    cat("\nShrinkage above 20% on ", paste(x$coordinates[high],
      collapse = ", "), ": at that level neither route reliably\n",
      "discriminates, and a selection here should not be acted on.\n",
      sep = "")
  if (any(!x$agree, na.rm = TRUE))
    cat("\nThe two routes disagree. That disagreement is the finding: the\n",
      "pooled ranking is biased toward the incumbent, the identity is\n",
      "unbiased but noisy, and neither is trustworthy where they differ.\n",
      sep = "")
  if (is.null(x$identityTable))
    cat("\nThe identity route did not return a table; only the pooled ranking\n",
      "is available, and it cannot support a BIC or a significance claim.\n",
      sep = "")
  cat("\nPooled log-likelihoods are over draws, not subjects. They rank\n",
    "candidates; they are not a likelihood for any test or penalty.\n",
    sep = "")
  invisible(x)
}
