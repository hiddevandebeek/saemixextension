## Eta-family selection around, not inside, the fixed-model score estimator.

copulaEtaCandidateMargin <- function(family, values) {
  scale <- sqrt(mean(values^2))
  if (!is.finite(scale) || scale <= 0) stop("posterior eta scale is invalid")
  switch(family,
    normal = copulaMarginNormal(scale),
    student = copulaMarginStudent(scale, 6),
    laplace = copulaMarginLaplace(scale),
    stop("unsupported automatic eta family: ", family))
}

copulaPosteriorMcmcEss <- function(draws) {
  samples <- dim(draws)[3L]
  value <- numeric(dim(draws)[1L] * dim(draws)[2L])
  cursor <- 0L
  for (subject in seq_len(dim(draws)[1L]))
    for (coordinate in seq_len(dim(draws)[2L])) {
      cursor <- cursor + 1L
      x <- draws[subject, coordinate, ]
      if (!is.finite(stats::var(x)) || stats::var(x) <= 0) {
        value[cursor] <- 0
        next
      }
      correlation <- as.numeric(stats::acf(x, plot = FALSE,
        lag.max = min(50L, samples - 1L), demean = TRUE)$acf)[-1L]
      stop <- which(!is.finite(correlation) | correlation <= 0)
      if (length(stop)) correlation <- correlation[seq_len(stop[1L] - 1L)]
      value[cursor] <- samples / (1 + 2 * sum(correlation))
    }
  pmin(samples, pmax(1, value))
}

copulaPosteriorEtaDraws <- function(object, nsamp, max.iter, seed,
                                    thin = 1L) {
  state <- copulaGet(object)
  index <- as.integer(state$etaIndex)
  n <- object["data"]["N"]
  dEta <- length(index)
  centre <- object["results"]["mean.phi"][, index, drop = FALSE]
  current <- object["results"]["phi"][, index, drop = FALSE] - centre
  if (any(!is.finite(current))) current[,] <- 0
  conditioning <- if (state$dConditioning > 0L)
    as.matrix(state$conditioning) else NULL
  if (is.null(conditioning)) {
    priorNegative <- function(eta) -copulaGaussianFremLogPrior(eta,
      state$vine, state$margins, state$dEta, "joint")
    priorRandom <- function() copulaMarginsQuantile(
      rvinecopulib::rvinecop(n, state$vine), state$margins)
  } else {
    kernel <- copulaGaussianFremConditionalKernel(conditioning,
      state$vine, state$margins, state$dEta)
    priorNegative <- kernel$negative
    priorRandom <- kernel$random
  }
  response <- function(eta) {
    phi <- object["results"]["mean.phi"]
    phi[, index] <- centre + eta
    as.numeric(copulaResponseLogLikBatch(object, phi, 1L))
  }
  proposal <- as.matrix(state$proposalOmega)
  if (any(dim(proposal) != c(dEta, dEta)) ||
      inherits(try(chol(proposal), silent = TRUE), "try-error"))
    proposal <- diag(copulaMarginScales(
      state$margins[seq_len(dEta)])^2, dEta)
  proposalChol <- chol(proposal) * .45
  burn <- if (is.null(max.iter)) max(100L,
    min(500L, as.integer(sum(object["options"]$nbiter.saemix) / 4L))) else
    max(50L, as.integer(max.iter))
  thin <- max(1L, as.integer(thin))
  draws <- array(NA_real_, c(n, dEta, nsamp))
  independentAccepted <- randomWalkAccepted <- integer(n)
  withSeed(seed, {
    responseCurrent <- response(current)
    priorCurrent <- priorNegative(current)
    retained <- 0L
    total <- burn + nsamp * thin
    for (iteration in seq_len(total)) {
      independent <- priorRandom()
      responseCandidate <- response(independent)
      priorIndependent <- priorNegative(independent)
      accept <- is.finite(responseCandidate) &
        log(stats::runif(n)) < pmin(0, responseCandidate - responseCurrent)
      if (any(accept)) {
        independentAccepted <- independentAccepted + as.integer(accept)
        current[accept, ] <- independent[accept, , drop = FALSE]
        responseCurrent[accept] <- responseCandidate[accept]
        priorCurrent[accept] <- priorIndependent[accept]
      }

      randomWalk <- current + matrix(stats::rnorm(n * dEta), ncol = dEta) %*%
        proposalChol
      responseCandidate <- response(randomWalk)
      priorCandidate <- priorNegative(randomWalk)
      logRatio <- responseCandidate - priorCandidate -
        responseCurrent + priorCurrent
      accept <- is.finite(logRatio) &
        log(stats::runif(n)) < pmin(0, logRatio)
      if (any(accept)) {
        randomWalkAccepted <- randomWalkAccepted + as.integer(accept)
        current[accept, ] <- randomWalk[accept, , drop = FALSE]
        responseCurrent[accept] <- responseCandidate[accept]
        priorCurrent[accept] <- priorCandidate[accept]
      }
      if (iteration > burn && (iteration - burn) %% thin == 0L) {
        retained <- retained + 1L
        draws[, , retained] <- current
      }
    }
  })
  mcmcEss <- copulaPosteriorMcmcEss(draws)
  attr(draws, "diagnostics") <- list(
    burn = burn, thin = thin,
    independentAcceptance = independentAccepted / total,
    randomWalkAcceptance = randomWalkAccepted / total,
    minimumMcmcEssFraction = min(mcmcEss) / nsamp,
    medianMcmcEssFraction = stats::median(mcmcEss) / nsamp)
  draws
}

copulaFlattenEtaDraws <- function(draws, conditioning = NULL) {
  n <- dim(draws)[1L]
  dEta <- dim(draws)[2L]
  samples <- dim(draws)[3L]
  eta <- do.call(rbind, lapply(seq_len(samples), function(sample)
    matrix(draws[, , sample], nrow = n, ncol = dEta)))
  if (is.null(conditioning)) return(eta)
  conditioning <- as.matrix(conditioning)
  if (nrow(conditioning) != n)
    stop("posterior draws and conditioning rows do not align")
  cbind(eta, conditioning[rep(seq_len(n), samples), , drop = FALSE])
}

copulaPosteriorBridgeMetricFromLogWeight <- function(logWeight) {
  logWeight <- as.matrix(logWeight)
  nSubject <- nrow(logWeight)
  nSample <- ncol(logWeight)
  maximum <- apply(logWeight, 1L, max)
  valid <- is.finite(maximum)
  ratio <- rep(-Inf, nSubject)
  ess <- rep(0, nSubject)
  relativeVariance <- rep(Inf, nSubject)
  if (any(valid)) {
    shifted <- exp(logWeight[valid, , drop = FALSE] - maximum[valid])
    average <- rowMeans(shifted)
    ratio[valid] <- maximum[valid] + log(average)
    ess[valid] <- rowSums(shifted)^2 / rowSums(shifted^2)
    batchLength <- max(2L, floor(sqrt(nSample)))
    nBatch <- floor(nSample / batchLength)
    if (nBatch >= 2L) {
      used <- seq_len(nBatch * batchLength)
      batchMean <- vapply(seq_len(nBatch), function(batch) {
        columns <- used[(batch - 1L) * batchLength + seq_len(batchLength)]
        rowMeans(shifted[, columns, drop = FALSE])
      }, numeric(sum(valid)))
      if (nBatch == 1L) batchMean <- matrix(batchMean, ncol = 1L)
      relativeVariance[valid] <- apply(batchMean, 1L, stats::var) /
        (nBatch * average^2)
    } else {
      relativeVariance[valid] <- apply(shifted, 1L, stats::var) /
        (nSample * average^2)
    }
  }
  list(deltaLogLik = sum(ratio), mcse = sqrt(sum(relativeVariance)),
    minimumEssFraction = min(ess) / nSample,
    medianEssFraction = stats::median(ess) / nSample,
    perSubject = ratio)
}

copulaPosteriorBridgeCache <- function(state, values, nSubject, nSample) {
  dEta <- state$dEta
  conditioningMargins <- if (state$dConditioning > 0L)
    state$margins[dEta + seq_len(state$dConditioning)] else list()
  fast <- all(vapply(conditioningMargins, function(margin)
    identical(margin$type, "continuous"), logical(1)))
  if (!fast) {
    baseline <- copulaGaussianFremLogPrior(values, state$vine,
      state$margins, dEta, "joint")
    return(list(fast = FALSE, values = values, baseline = baseline,
      nSubject = nSubject, nSample = nSample))
  }
  eta <- values[, seq_len(dEta), drop = FALSE]
  baselineEta <- copulaGaussianFremEvaluateMargins(eta,
    state$margins[seq_len(dEta)])
  if (!all(baselineEta$valid)) stop("incumbent posterior eta draw is invalid")
  conditioning <- if (state$dConditioning > 0L)
    values[, dEta + seq_len(state$dConditioning), drop = FALSE] else
    matrix(numeric(), nrow(values), 0L)
  pattern <- if (ncol(conditioning))
    apply(!is.na(conditioning), 1L, paste0, collapse = "") else
    rep("", nrow(values))
  correlation <- copulaGaussianRvineCor(state$vine, state$d)
  patterns <- lapply(unique(pattern), function(key) {
    rows <- which(pattern == key)
    observed <- if (ncol(conditioning))
      which(!is.na(conditioning[rows[1L], ])) else integer()
    index <- c(seq_len(dEta), dEta + observed)
    K <- solve(correlation[index, index, drop = FALSE]) - diag(length(index))
    zCov <- if (length(observed))
      copulaGaussianFremEvaluateMargins(
        conditioning[rows, observed, drop = FALSE],
        state$margins[dEta + observed])$z else
      matrix(numeric(), length(rows), 0L)
    list(rows = rows, etaK = K[seq_len(dEta), seq_len(dEta), drop = FALSE],
      crossK = if (length(observed)) K[seq_len(dEta),
        dEta + seq_along(observed), drop = FALSE] else
        matrix(numeric(), dEta, 0L), zCov = zCov)
  })
  contribution <- function(evaluated, item) {
    rows <- item$rows
    z <- evaluated$z[rows, , drop = FALSE]
    quadratic <- rowSums((z %*% item$etaK) * z)
    cross <- if (ncol(item$crossK))
      rowSums((z %*% item$crossK) * item$zCov) else 0
    -.5 * quadratic - cross +
      rowSums(evaluated$logMargin[rows, , drop = FALSE])
  }
  baselineContribution <- numeric(nrow(values))
  for (item in patterns)
    baselineContribution[item$rows] <- contribution(baselineEta, item)
  list(fast = TRUE, eta = eta, patterns = patterns,
    baselineContribution = baselineContribution, contribution = contribution,
    nSubject = nSubject, nSample = nSample)
}

copulaPosteriorBridgeCachedMetric <- function(etaMargins, cache, state) {
  if (!cache$fast) {
    margins <- c(etaMargins, state$margins[-seq_len(state$dEta)])
    candidate <- copulaGaussianFremLogPrior(cache$values, state$vine,
      margins, state$dEta, "joint")
    logWeight <- candidate - cache$baseline
  } else {
    evaluated <- copulaGaussianFremEvaluateMargins(cache$eta, etaMargins)
    contribution <- rep(-Inf, nrow(cache$eta))
    if (any(evaluated$valid)) for (item in cache$patterns) {
      rows <- item$rows
      validRows <- rows[evaluated$valid[rows]]
      if (length(validRows) == length(rows))
        contribution[rows] <- cache$contribution(evaluated, item)
    }
    logWeight <- contribution - cache$baselineContribution
  }
  copulaPosteriorBridgeMetricFromLogWeight(matrix(logWeight,
    nrow = cache$nSubject, ncol = cache$nSample))
}

copulaEtaFamilyGrid <- function(candidates, dEta, maxModels) {
  allowed <- c("normal", "student", "laplace")
  if (is.character(candidates)) candidates <- rep(list(candidates), dEta)
  if (!is.list(candidates) || length(candidates) != dEta)
    stop("candidates must be a character vector or one character vector per eta")
  candidates <- lapply(candidates, function(value) {
    value <- unique(tolower(as.character(value)))
    if (!length(value) || any(!value %in% allowed))
      stop("eta candidates must be Normal, Student, or Laplace")
    value
  })
  count <- prod(lengths(candidates))
  if (count > maxModels)
    stop("eta candidate grid has ", count, " models; reduce candidates or raise maxModels")
  grid <- expand.grid(candidates, stringsAsFactors = FALSE)
  names(grid) <- paste0("eta", seq_len(dEta))
  grid
}

copulaFitPosteriorEtaCandidate <- function(families, state, trainDraws,
                                           validationDraws, maxit) {
  dEta <- state$dEta
  n <- dim(trainDraws)[1L]
  trainSample <- dim(trainDraws)[3L]
  validationSample <- dim(validationDraws)[3L]
  start <- lapply(seq_len(dEta), function(j)
    copulaEtaCandidateMargin(families[j], as.numeric(trainDraws[, j, ])))
  layout <- copulaMarginLayout(start)
  internal <- copulaScoreToInternal(layout$par, layout$lower, layout$upper)
  conditioning <- state$conditioning
  trainValues <- copulaFlattenEtaDraws(trainDraws, conditioning)
  validationValues <- copulaFlattenEtaDraws(validationDraws, conditioning)
  trainCache <- copulaPosteriorBridgeCache(state, trainValues, n, trainSample)
  validationCache <- copulaPosteriorBridgeCache(state, validationValues, n,
    validationSample)
  materialize <- function(value) {
    native <- copulaScoreFromInternal(value, layout$lower, layout$upper)
    copulaMarginsWithParameters(start, layout, native)
  }
  objective <- function(value) {
    margins <- try(materialize(value), silent = TRUE)
    if (inherits(margins, "try-error")) return(.Machine$double.xmax / 100)
    metric <- try(copulaPosteriorBridgeCachedMetric(margins, trainCache,
      state), silent = TRUE)
    if (inherits(metric, "try-error") || !is.finite(metric$deltaLogLik))
      return(.Machine$double.xmax / 100)
    -metric$deltaLogLik
  }
  fit <- stats::optim(internal, objective, method = "BFGS",
    control = list(maxit = as.integer(maxit), reltol = 1e-7))
  etaMargins <- materialize(fit$par)
  margins <- c(etaMargins, state$margins[-seq_len(dEta)])
  training <- copulaPosteriorBridgeCachedMetric(etaMargins, trainCache, state)
  validation <- copulaPosteriorBridgeCachedMetric(etaMargins,
    validationCache, state)
  incumbentParameters <- sum(vapply(state$margins[seq_len(dEta)],
    function(margin) sum(margin$free), integer(1)))
  candidateParameters <- sum(vapply(margins[seq_len(dEta)],
    function(margin) sum(margin$free), integer(1)))
  deltaParameters <- candidateParameters - incumbentParameters
  list(families = families, margins = margins, convergence = fit$convergence,
    training = training, validation = validation,
    deltaParameters = deltaParameters,
    validationBicAdvantage = 2 * validation$deltaLogLik -
      deltaParameters * log(n))
}

copulaPopulationFromEtaSelection <- function(state, margins) {
  conditioning <- if (state$dConditioning > 0L)
    list(values = state$conditioning,
      variableName = state$conditioningName) else NULL
  copulaPopulation(state$vine, margins = margins,
    scale = "transformed-additive", conditioning = conditioning,
    scoreScale = state$scoreScale %||% 1,
    scoreFiniteDifference = state$scoreFiniteDifference %||% 1e-4,
    scoreProjection = state$scoreProjection %||% 24,
    scoreGainScale = state$scoreGainScale %||% .2,
    scoreGainPower = state$scoreGainPower %||% .8,
    scoreGainOffset = state$scoreGainOffset %||% 30,
    scoreBurn = state$scoreBurn %||% 50L)
}

copulaSelectEtaMargins <- function(object,
    candidates = c("normal", "student", "laplace"),
    posteriorDraws = 500L, max.iter = NULL, seed = 930001L,
    maxModels = 81L, optimizerMaxit = 100L,
    minimumEssFraction = .01, minimumPosteriorEssFraction = .05,
    maximumMcse = .5) {
  if (!inherits(object, "SaemixObject"))
    stop("object must be a fitted SaemixObject")
  state <- copulaGet(object)
  if (!copulaIsFullGaussianVine(state$vine, state$d) ||
      !identical(state$populationAlgorithm, "score-sa"))
    stop("eta selection requires a fitted fixed Gaussian-copula FREM model")
  posteriorDraws <- as.integer(posteriorDraws)
  if (length(posteriorDraws) != 1L || is.na(posteriorDraws) ||
      posteriorDraws < 50L)
    stop("posteriorDraws must be at least 50 per independent pool")
  grid <- copulaEtaFamilyGrid(candidates, state$dEta,
    as.integer(maxModels))
  training <- copulaPosteriorEtaDraws(object, posteriorDraws, max.iter,
    as.integer(seed))
  validation <- copulaPosteriorEtaDraws(object, posteriorDraws, max.iter,
    as.integer(seed) + 1L)
  trainingMcmc <- attr(training, "diagnostics")
  validationMcmc <- attr(validation, "diagnostics")
  fitted <- lapply(seq_len(nrow(grid)), function(row)
    copulaFitPosteriorEtaCandidate(as.character(grid[row, ]), state,
      training, validation, optimizerMaxit))
  table <- do.call(rbind, lapply(seq_along(fitted), function(i) {
    result <- fitted[[i]]
    data.frame(model = i,
      families = paste(result$families, collapse = "/"),
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
  if (!length(eligible))
    stop("no eta-family candidate passed the posterior-overlap diagnostics")
  selected <- eligible[which.max(table$validation_bic_advantage[eligible])]
  alternatives <- setdiff(eligible, selected)
  selectionGap <- if (length(alternatives))
    table$validation_bic_advantage[selected] -
      max(table$validation_bic_advantage[alternatives]) else Inf
  runnerUp <- if (length(alternatives)) alternatives[
    which.max(table$validation_bic_advantage[alternatives])] else NA_integer_
  selectionMcse <- if (is.na(runnerUp)) 0 else 2 * sqrt(
    table$validation_mcse[selected]^2 + table$validation_mcse[runnerUp]^2)
  allCandidatesEligible <- all(table$eligible)
  posteriorMixingAdequate <-
    trainingMcmc$minimumMcmcEssFraction >= minimumPosteriorEssFraction &&
    validationMcmc$minimumMcmcEssFraction >= minimumPosteriorEssFraction
  selectionResolved <- allCandidatesEligible && posteriorMixingAdequate &&
    selectionGap > 2 * selectionMcse
  populations <- lapply(fitted, function(result)
    copulaPopulationFromEtaSelection(state, result$margins))
  answer <- list(selected = selected, table = table,
    families = fitted[[selected]]$families,
    margins = fitted[[selected]]$margins,
    population = populations[[selected]], populations = populations,
    posteriorDrawsPerPool = posteriorDraws,
    trainingSeed = as.integer(seed), validationSeed = as.integer(seed) + 1L,
    diagnostics = list(minimumEssFraction = minimumEssFraction,
      maximumMcse = maximumMcse,
      minimumPosteriorEssFraction = minimumPosteriorEssFraction,
      trainingMcmc = trainingMcmc, validationMcmc = validationMcmc,
      posteriorMixingAdequate = posteriorMixingAdequate,
      selectionGap = selectionGap, selectionMcse = selectionMcse,
      selectionResolved = selectionResolved,
      allCandidatesEligible = allCandidatesEligible,
      fixedDuringFinalFit = TRUE,
      selectionOutsideScoreSa = TRUE))
  class(answer) <- c("saemixEtaMarginSelection", "list")
  if (!selectionResolved)
    warning(paste0("eta-family ranking is unresolved because of Monte Carlo ",
      "error, candidate overlap, or posterior mixing; increase posteriorDraws ",
      "or max.iter"),
      call. = FALSE)
  answer
}

print.saemixEtaMarginSelection <- function(x, ...) {
  cat("Posterior eta-margin selection\n")
  cat("  selected:", paste(x$families, collapse = " / "), "\n")
  cat("  posterior draws per independent pool:",
    x$posteriorDrawsPerPool, "\n")
  print(x$table[order(-x$table$validation_bic_advantage), ], row.names = FALSE)
  invisible(x)
}
