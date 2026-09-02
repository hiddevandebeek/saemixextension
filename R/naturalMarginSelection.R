## Automatic natural-parameter family screening from full posterior draws of a
## standard transformed-additive incumbent. Families are support-based and do
## not depend on the saemix coordinate transform.

copulaNaturalFamilyGrid <- function(candidates, supports, maxModels) {
  d <- length(supports)
  if (is.null(candidates)) candidates <- lapply(supports,
    copulaNaturalMarginFamilies)
  if (is.character(candidates)) candidates <- rep(list(candidates), d)
  if (!is.list(candidates) || length(candidates) != d)
    stop("candidates must contain one family vector per parameter")
  registry <- copulaNaturalMarginRegistry()
  for (j in seq_len(d)) {
    candidates[[j]] <- unique(tolower(candidates[[j]]))
    if (!length(candidates[[j]]) || any(!candidates[[j]] %in% names(registry)))
      stop("candidate family is not registered")
    bad <- vapply(candidates[[j]], function(name)
      !identical(registry[[name]]$support, supports[j]), logical(1))
    if (any(bad)) stop("candidate family is incompatible with parameter support")
  }
  if (prod(lengths(candidates)) > maxModels)
    stop("natural-margin candidate grid exceeds maxModels")
  grid <- expand.grid(candidates, stringsAsFactors = FALSE)
  names(grid) <- paste0("parameter", seq_len(d)); grid
}

copulaNaturalPosteriorData <- function(object, draws) {
  state <- copulaGet(object); index <- as.integer(state$etaIndex)
  n <- dim(draws)[1L]; samples <- dim(draws)[3L]
  predictorSubject <- object["results"]["mean.phi"][, index, drop = FALSE]
  predictor <- predictorSubject[rep(seq_len(n), samples), , drop = FALSE]
  eta <- copulaFlattenEtaDraws(draws)
  phi <- predictor + eta
  transform <- as.integer(object["model"]["transform.par"][index])
  conditioning <- if ((state$dConditioning %||% 0L) > 0L)
    as.matrix(state$conditioning)[rep(seq_len(n), samples), , drop = FALSE] else
    matrix(numeric(), n * samples, 0L)
  list(n = n, samples = samples, eta = eta, predictor = predictor,
    typical = copulaWorkingToNatural(predictor, transform),
    natural = copulaWorkingToNatural(phi, transform), transform = transform,
    conditioning = conditioning)
}

copulaNaturalBridgePrepare <- function(data, incumbent) {
  hasConditioning <- (incumbent$dConditioning %||% 0L) > 0L
  baselineState <- if (hasConditioning)
    cbind(data$eta, data$conditioning) else data$eta
  baseline <- copulaGaussianFremLogPrior(baselineState, incumbent$vine,
    incumbent$margins, incumbent$dEta, "joint")
  R <- copulaGaussianRvineCor(incumbent$vine, incumbent$d)
  U <- chol(R)
  covariate <- if (hasConditioning) copulaGaussianFremEvaluateMargins(
    data$conditioning, incumbent$margins[incumbent$dEta +
      seq_len(incumbent$dConditioning)]) else list(
        z = matrix(numeric(), nrow(data$eta), 0L),
        logMargin = matrix(numeric(), nrow(data$eta), 0L),
        valid = rep(TRUE, nrow(data$eta)))
  if (any(!is.finite(baseline)) || !all(covariate$valid))
    stop("incumbent posterior pool has an invalid population density")
  data$bridgeCache <- list(baseline = baseline, covariate = covariate,
    U = U, logDet = 2 * sum(log(diag(U))),
    copulaQuadratic = solve(R) - diag(incumbent$d),
    logJacobian = rowSums(copulaWorkingLogJacobian(
      data$predictor + data$eta, data$transform)))
  data
}

copulaNaturalBridgeCandidateLogPrior <- function(margins, data, incumbent) {
  cache <- data$bridgeCache
  if (is.null(cache)) data <- copulaNaturalBridgePrepare(data, incumbent)
  cache <- data$bridgeCache
  evaluated <- copulaNaturalMarginsEvaluate(data$natural, data$typical, margins)
  valid <- evaluated$valid & cache$covariate$valid
  answer <- rep(-Inf, nrow(data$eta)); rows <- which(valid)
  if (!length(rows)) return(answer)
  z <- cbind(evaluated$z, cache$covariate$z)
  zr <- z[rows, , drop = FALSE]
  logCopula <- -.5 * cache$logDet - .5 * rowSums(
    (zr %*% cache$copulaQuadratic) * zr)
  answer[rows] <- logCopula +
    rowSums(evaluated$logMargin[rows, , drop = FALSE]) +
    rowSums(cache$covariate$logMargin[rows, , drop = FALSE]) +
    cache$logJacobian[rows]
  answer
}

copulaNaturalBridgeMetric <- function(margins, data, incumbent) {
  if (is.null(data$bridgeCache)) data <- copulaNaturalBridgePrepare(data, incumbent)
  candidate <- copulaNaturalBridgeCandidateLogPrior(margins, data, incumbent)
  baseline <- data$bridgeCache$baseline
  copulaPosteriorBridgeMetricFromLogWeight(matrix(candidate - baseline,
    nrow = data$n, ncol = data$samples))
}

copulaFitNaturalCandidate <- function(families, supports, train, validation,
                                      incumbent, maxit) {
  start <- lapply(seq_along(families), function(j)
    copulaNaturalMarginStart(families[j], train$natural[, j],
      train$typical[, j], supports[j]))
  layout <- copulaMarginLayout(start)
  internal <- copulaScoreToInternal(layout$par, layout$lower, layout$upper)
  materialize <- function(value) copulaMarginsWithParameters(start, layout,
    copulaScoreFromInternal(value, layout$lower, layout$upper))
  objective <- function(value) {
    margins <- try(materialize(value), silent = TRUE)
    if (inherits(margins, "try-error")) return(.Machine$double.xmax / 100)
    metric <- try(copulaNaturalBridgeMetric(margins, train, incumbent),
      silent = TRUE)
    if (inherits(metric, "try-error") || !is.finite(metric$deltaLogLik))
      return(.Machine$double.xmax / 100)
    -metric$deltaLogLik
  }
  fit <- stats::optim(internal, objective, method = "BFGS",
    control = list(maxit = as.integer(maxit), reltol = 1e-7))
  margins <- materialize(fit$par)
  training <- copulaNaturalBridgeMetric(margins, train, incumbent)
  validationMetric <- copulaNaturalBridgeMetric(margins, validation, incumbent)
  incumbentParameters <- sum(vapply(incumbent$margins[seq_len(incumbent$dEta)],
    function(margin) sum(margin$free), integer(1)))
  candidateParameters <- sum(vapply(margins,
    function(margin) sum(margin$free), integer(1)))
  deltaParameters <- candidateParameters - incumbentParameters
  list(families = families, margins = margins, convergence = fit$convergence,
    training = training, validation = validationMetric,
    deltaParameters = deltaParameters,
    validationBicAdvantage = 2 * validationMetric$deltaLogLik -
      deltaParameters * log(train$n))
}

copulaSelectParameterMargins <- function(object, supports, candidates = NULL,
    posteriorDraws = 500L, max.iter = NULL, seed = 935001L,
    maxModels = 81L, optimizerMaxit = 100L,
    minimumEssFraction = .01, minimumPosteriorEssFraction = .05,
    maximumMcse = .5) {
  if (!inherits(object, "SaemixObject")) stop("object must be a fitted SaemixObject")
  state <- copulaGet(object)
  if (!identical(state$populationScale, "transformed-additive") ||
      !copulaIsFullGaussianVine(state$vine, state$d))
    stop("natural-margin screening requires a standard transformed-additive incumbent")
  if ((state$dConditioning %||% 0L) > 0L &&
      (anyNA(state$conditioning) || any(vapply(
        state$margins[state$dEta + seq_len(state$dConditioning)],
        function(m) !identical(m$type, "continuous"), logical(1)))))
    stop("natural-margin screening with covariates requires complete continuous conditioning")
  supports <- rep_len(as.character(supports), state$dEta)
  if (any(!supports %in% c("real", "positive", "unit")))
    stop("supports must be real, positive, or unit")
  grid <- copulaNaturalFamilyGrid(candidates, supports, as.integer(maxModels))
  trainingDraws <- copulaPosteriorEtaDraws(object, posteriorDraws, max.iter,
    as.integer(seed))
  validationDraws <- copulaPosteriorEtaDraws(object, posteriorDraws, max.iter,
    as.integer(seed) + 1L)
  training <- copulaNaturalPosteriorData(object, trainingDraws)
  validation <- copulaNaturalPosteriorData(object, validationDraws)
  training <- copulaNaturalBridgePrepare(training, state)
  validation <- copulaNaturalBridgePrepare(validation, state)
  fitted <- lapply(seq_len(nrow(grid)), function(row)
    copulaFitNaturalCandidate(as.character(grid[row, ]), supports, training,
      validation, state, optimizerMaxit))
  table <- do.call(rbind, lapply(seq_along(fitted), function(i) {
    result <- fitted[[i]]
    data.frame(model = i, families = paste(result$families, collapse = "/"),
      convergence = result$convergence,
      delta_parameters = result$deltaParameters,
      training_delta_loglik = result$training$deltaLogLik,
      validation_delta_loglik = result$validation$deltaLogLik,
      validation_bic_advantage = result$validationBicAdvantage,
      validation_mcse = result$validation$mcse,
      minimum_ess_fraction = result$validation$minimumEssFraction,
      median_ess_fraction = result$validation$medianEssFraction,
      eligible = result$convergence == 0L &&
        is.finite(result$validationBicAdvantage) &&
        result$validation$minimumEssFraction >= minimumEssFraction &&
        result$validation$mcse <= maximumMcse)
  }))
  eligible <- which(table$eligible)
  if (!length(eligible)) stop("no natural-margin candidate passed overlap diagnostics")
  selected <- eligible[which.max(table$validation_bic_advantage[eligible])]
  alternatives <- setdiff(eligible, selected)
  gap <- if (length(alternatives)) table$validation_bic_advantage[selected] -
    max(table$validation_bic_advantage[alternatives]) else Inf
  runner <- if (length(alternatives)) alternatives[
    which.max(table$validation_bic_advantage[alternatives])] else NA_integer_
  selectionMcse <- if (is.na(runner)) 0 else 2 * sqrt(
    table$validation_mcse[selected]^2 + table$validation_mcse[runner]^2)
  trainDiag <- attr(trainingDraws, "diagnostics")
  validationDiag <- attr(validationDraws, "diagnostics")
  mixing <- trainDiag$minimumMcmcEssFraction >= minimumPosteriorEssFraction &&
    validationDiag$minimumMcmcEssFraction >= minimumPosteriorEssFraction
  resolved <- all(table$eligible) && mixing && gap > 2 * selectionMcse
  margins <- fitted[[selected]]$margins
  finalMargins <- if ((state$dConditioning %||% 0L) > 0L)
    c(margins, state$margins[state$dEta + seq_len(state$dConditioning)]) else
    margins
  population <- copulaPopulation(state$vine, margins = finalMargins,
    scale = "parameter",
    conditioning = if ((state$dConditioning %||% 0L) > 0L)
      list(values = state$conditioning,
        variableName = state$conditioningName) else NULL,
    populationAlgorithm = "score-sa", scoreScale = "auto",
    scoreBurn = state$scoreBurn %||% 50L,
    scoreGainScale = state$scoreGainScale %||% .2,
    scoreGainPower = state$scoreGainPower %||% .8,
    scoreFiniteDifference = state$scoreFiniteDifference %||% 1e-4,
    scoreProjection = state$scoreProjection %||% 24)
  answer <- list(selected = selected, table = table,
    families = fitted[[selected]]$families, margins = margins,
    population = population, posteriorDrawsPerPool = posteriorDraws,
    diagnostics = list(selectionResolved = resolved,
      posteriorMixingAdequate = mixing, selectionGap = gap,
      selectionMcse = selectionMcse, allCandidatesEligible = all(table$eligible),
      fixedDuringFinalFit = TRUE, selectionOutsideScoreSa = TRUE))
  class(answer) <- c("saemixNaturalMarginSelection", "list")
  if (!resolved) warning("natural-margin ranking is unresolved; increase posteriorDraws or max.iter",
    call. = FALSE)
  answer
}

print.saemixNaturalMarginSelection <- function(x, ...) {
  cat("Natural parameter-margin selection\n  selected:",
    paste(x$families, collapse = " / "), "\n")
  print(x$table[order(-x$table$validation_bic_advantage), ], row.names = FALSE)
  invisible(x)
}
