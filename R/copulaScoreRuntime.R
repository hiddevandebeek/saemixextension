## Runtime integration of the pure score step with saemix.

copulaScoreBatchUpdate <- function(eta, nchains, phi = NULL, design = NULL,
                                   locationMap = NULL, beta = NULL,
                                   betaFree = NULL, subject = NULL) {
  eta <- as.matrix(eta)
  if (!nrow(eta) || any(!is.finite(eta)))
    stop("score batch requires finite eta rows")
  nchains <- as.integer(nchains)
  if (length(nchains) != 1L || is.na(nchains) || nchains < 1L)
    stop("score batch requires a positive chain count")
  if (is.null(subject) || length(subject) != nrow(eta) ||
      any(!is.finite(subject)) || any(subject != as.integer(subject)) ||
      any(subject < 1L))
    stop("score batch requires one positive integer subject label per row")
  subject <- as.integer(subject)
  if (any(tabulate(match(subject, unique(subject))) != nchains))
    stop("score batch requires every subject exactly once per chain")
  .cop$curEta <- eta

  conditioning <- NULL
  if ((.cop$dConditioning %||% 0L) > 0L) {
    if (nrow(.cop$conditioning) < max(subject))
      stop("score batch conditioning matrix has fewer rows than subjects")
    conditioning <- .cop$conditioning[subject, , drop = FALSE]
  }
  .cop$curConditioning <- conditioning
  conditioningMissing <- if (is.null(conditioning)) NULL else is.na(conditioning)

  designParts <- c(phi = !is.null(phi), design = !is.null(design),
    locationMap = !is.null(locationMap), beta = !is.null(beta))
  if (any(designParts) && !all(designParts))
    stop("score batch requires phi, design, locationMap, and beta together")
  hasDesign <- all(designParts)
  if (hasDesign) {
    phi <- as.matrix(phi)
    design <- as.matrix(design)
    locationMap <- as.matrix(locationMap)
    if (nrow(phi) != nrow(eta) || nrow(design) != nrow(eta) ||
        ncol(design) != nrow(locationMap) ||
        ncol(phi) != ncol(locationMap) || length(beta) != ncol(design))
      stop("score batch has incompatible phi/design/location dimensions")
    .cop$curPhi <- phi
    .cop$curX <- design
    .cop$locMap <- locationMap
    .cop$betaCurrent <- as.numeric(beta)
    .cop$betaFree <- if (is.null(betaFree)) seq_along(beta) else
      as.integer(betaFree)
    if (!isTRUE(.cop$locationRankChecked)) {
      active <- intersect(.cop$betaFree,
        which(rowSums(abs(locationMap)) > 0))
      if (length(active)) {
        effective <- vapply(active, function(j)
          as.vector(outer(design[, j], locationMap[j, ])),
          numeric(nrow(design) * ncol(locationMap)))
        if (qr(effective)$rank < length(active))
          stop("score batch beta-to-eta location design is rank deficient")
      }
      .cop$locationRankChecked <- TRUE
    }
  } else {
    .cop$curPhi <- NULL
    .cop$curX <- NULL
    .cop$locMap <- NULL
    .cop$betaCurrent <- NULL
    .cop$betaFree <- NULL
  }

  hasCategorical <- !is.null(conditioning) && any(vapply(
    .cop$margins[(.cop$dEta + 1L):.cop$d], function(m)
      identical(m$type, "discrete"), logical(1)))
  .cop$curCategoricalUniform <- NULL
  if (hasCategorical) {
    augmented <- copulaGaussianFremAugmentMixedConditioning(
      eta, conditioning, .cop$vine, .cop$margins,
      .cop$dEta %||% ncol(eta))
    conditioning <- augmented$conditioning
    .cop$curCategoricalUniform <- augmented$categoricalUniform
  } else if (copulaIsFullGaussianVine(.cop$vine, .cop$d) &&
      isTRUE(.cop$augmentMissingGaussian) &&
      !is.null(conditioning) && anyNA(conditioning))
    conditioning <- copulaGaussianFremImputeMissingConditioning(
      eta, conditioning, .cop$vine, .cop$margins,
      .cop$dEta %||% ncol(eta))
  .cop$curConditioningComplete <- conditioning
  referenceUniform <- matrix(NA_real_, nrow(eta), .cop$d)
  movingEta <- which(vapply(.cop$margins[seq_len(.cop$dEta)],
    copulaMarginHasMovingSupport, logical(1)))
  for (j in movingEta) {
    u <- .cop$margins[[j]]$cdf(eta[, j], .cop$margins[[j]]$parameters)
    if (any(!is.finite(u)) || any(u <= 0 | u >= 1))
      stop("moving-support eta could not be mapped to the fixed reference interior")
    referenceUniform[, j] <- u
  }
  if (!is.null(conditioning) && any(conditioningMissing)) {
    for (local in seq_len(ncol(conditioning))) {
      j <- .cop$dEta + local
      rows <- conditioningMissing[, local] &
        copulaMarginHasMovingSupport(.cop$margins[[j]])
      if (any(rows)) {
        u <- .cop$margins[[j]]$cdf(conditioning[rows, local],
          .cop$margins[[j]]$parameters)
        if (any(!is.finite(u)) || any(u <= 0 | u >= 1))
          stop("missing moving-support covariate left the fixed reference interior")
        referenceUniform[rows, j] <- u
      }
    }
  }
  .cop$curReferenceUniform <- referenceUniform
  invisible(NULL)
}

copulaScoreResidualIndices <- function(errorModel) {
  unlist(lapply(seq_along(errorModel), function(j) {
    base <- 2L * (j - 1L)
    switch(errorModel[j], constant = base + 1L,
      exponential = base + 1L, proportional = base + 2L,
      combined = base + c(1L, 2L), base + c(1L, 2L))
  }), use.names = FALSE)
}

copulaScoreResponseBlock <- function(phi, randomIndex, transform, id, x, y,
                                     errorModel, exponentialType,
                                     structuralModel, predictions, residual,
                                     residualSum, nchains, nobs,
                                     hasFixedOnly = FALSE) {
  if (isTRUE(hasFixedOnly))
    stop("score-sa requires every estimated structural parameter in the population location block")
  evaluate <- local({
    template <- phi
    function(phiRandom, candidateResidual) {
      candidatePhi <- template
      candidatePhi[, randomIndex] <- phiRandom
      candidatePrediction <- structuralModel(
        transphi(candidatePhi, transform), id, x)
      for (type in exponentialType)
        candidatePrediction[x$ytype == type] <- log(cutoff(
          candidatePrediction[x$ytype == type]))
      candidateSd <- error(candidatePrediction, candidateResidual, x$ytype)
      if (any(!is.finite(candidateSd)) || any(candidateSd <= 0))
        return(rep(-Inf, length(candidatePrediction)))
      -.5 * ((y - candidatePrediction) / candidateSd)^2 -
        log(candidateSd) - .5 * log(2 * pi)
    }
  })
  gradient <- local({
    evaluateLog <- evaluate
    observationId <- as.integer(id)
    function(phiRandom, candidateResidual, coordinates, step) {
      phiRandom <- as.matrix(phiRandom)
      coordinates <- as.integer(coordinates)
      if (!length(coordinates))
        return(matrix(0, nrow(phiRandom), ncol(phiRandom)))
      if (any(observationId < 1L | observationId > nrow(phiRandom)))
        stop("response-gradient subject index is incompatible with the latent batch")
      answer <- matrix(0, nrow(phiRandom), ncol(phiRandom))
      for (coordinate in coordinates) {
        plus <- minus <- phiRandom
        plus[, coordinate] <- plus[, coordinate] + step
        minus[, coordinate] <- minus[, coordinate] - step
        derivative <- (evaluateLog(plus, candidateResidual) -
          evaluateLog(minus, candidateResidual)) / (2 * step)
        if (any(!is.finite(derivative)))
          stop("non-finite fixed-reference response derivative")
        summed <- rowsum(matrix(derivative, ncol = 1L), observationId,
          reorder = FALSE)
        answer[as.integer(rownames(summed)), coordinate] <- summed[, 1L]
      }
      answer
    }
  })
  list(y = y, f = predictions, etype = x$ytype, pres = residual,
    free = copulaScoreResidualIndices(errorModel), evaluate = evaluate,
    gradient = gradient,
    batchResidualMle = if (length(errorModel) == 1L &&
      errorModel %in% c("constant", "exponential", "proportional"))
      sqrt(residualSum / (nchains * nobs)) else NA_real_)
}

## Apply one score update to the current controlled-MCMC batch.
copulaScoreMstep <- function(kiter, final = FALSE, response = NULL) {
  if (!identical(.cop$mode, "joint") || !isTRUE(.cop$modelFrozen))
    stop("score-sa population update requires a frozen joint model")
  hasDesign <- !is.null(.cop$curPhi) && !is.null(.cop$curX) &&
    !is.null(.cop$locMap) && !is.null(.cop$betaCurrent)
  values <- if (hasDesign) .cop$curPhi else .cop$curEta
  design <- if (hasDesign) .cop$curX else NULL
  locationMap <- if (hasDesign) .cop$locMap else NULL
  betaStart <- if (hasDesign)
    copulaAcceptedBeta(.cop$betaCurrent, .cop$betaFitted,
      .cop$locMap, .cop$betaFree) else NULL
  conditioning <- .cop$curConditioningComplete %||% .cop$curConditioning
  if ((.cop$dConditioning %||% 0L) > 0L) {
    if (is.null(conditioning) || nrow(conditioning) != nrow(values))
      stop("score-sa requires aligned current conditioning rows")
    values <- cbind(values, conditioning)
    if (hasDesign)
      locationMap <- cbind(locationMap,
        matrix(0, nrow = nrow(locationMap), ncol = .cop$dConditioning))
  }
  scoreBurn <- .cop$scoreBurn %||% 50L
  gain <- .cop$scoreGainScale *
    (kiter + .cop$scoreGainOffset)^(-.cop$scoreGainPower)
  answer <- copulaGaussianFremPopulationScoreStep(
    values, rep(1 / nrow(values), nrow(values)), .cop$margins, .cop$vine,
    .cop$d, .cop$dEta, gain,
    X = design, locMap = locationMap, beta0 = betaStart,
    betaFree = if (hasDesign) .cop$betaFree else NULL,
    state = .cop$scoreState, scoreScale = .cop$scoreScale,
    finiteDifference = .cop$scoreFiniteDifference,
    projection = .cop$scoreProjection,
    adaptMetric = kiter <= scoreBurn, average = kiter > scoreBurn,
    useAverage = isTRUE(final), response = response,
    categoricalUniform = .cop$curCategoricalUniform,
    referenceUniform = .cop$curReferenceUniform)

  .cop$scoreState <- answer$state
  .cop$margins <- answer$margins
  .cop$vine <- answer$vine
  .cop$sdPrev <- answer$sd
  .cop$sd <- answer$sd
  .cop$delta <- answer$delta[seq_len(.cop$dEta)]
  .cop$betaJoint <- answer$beta
  .cop$betaFitted <- answer$beta
  .cop$residualJoint <- answer$residual
  copulaAssertFrozen(.cop$vine)

  theory <- answer$scoreTheory
  theory$proposalScaleFrozen <-
    !is.null(.cop$rwProposalFrozen) && kiter > scoreBurn
  theory$blockScheduleFrozen <-
    !is.null(.cop$rwBlockSizeFrozen) && kiter > scoreBurn
  theory$individualKernel <- if (isTRUE(.cop$multiCategoricalExactKernel))
    "exact conditional-prior independence; rectangle-dependent random walks disabled" else
    "exact conditional-prior independence plus invariant random-walk kernels"
  theory$rectangleProbabilityUsedInScore <- FALSE
  theory$gainPower <- .cop$scoreGainPower
  theory$fortH6CompatibleForLipschitzP2 <-
    is.finite(.cop$scoreGainPower) && .cop$scoreGainPower > .75 &&
    .cop$scoreGainPower <= 1
  theory$runtimeConditionsObserved <-
    isTRUE(theory$runtimeConditionsObserved) &&
    isTRUE(theory$proposalScaleFrozen) &&
    isTRUE(theory$blockScheduleFrozen) &&
    isTRUE(theory$fortH6CompatibleForLipschitzP2)
  theory$convergenceClaim <- paste0(
    "Conditional on the stated smoothness, controlled-Markov, moment, ",
    "stability, critical-value, and summable-random-error assumptions; ",
    "runtime diagnostics are necessary evidence, not a proof of those ",
    "assumptions.")

  .cop$lastJoint <- list(kiter = kiter, conv = 0L,
    value = answer$value, final = isTRUE(final), backend = answer$backend,
    scoreMax = answer$scoreMax, score = answer$score,
    scoreAverageMax = answer$scoreAverageMax,
    scoreAverage = answer$scoreAverage,
    gain = answer$gain, scoreScale = answer$scoreScale,
    finiteDifference = answer$finiteDifference,
    scoreMethod = answer$scoreMethod,
    numericalScoreComponents = answer$numericalScoreComponents,
    residual = answer$residual,
    batchResidualMle = response$batchResidualMle %||% NA_real_,
    preconditionerMin = answer$preconditionerMin,
    preconditionerMax = answer$preconditionerMax,
    metricGain = answer$metricGain,
    averaging = answer$averaging, usedAverage = answer$usedAverage,
    projectionCount = answer$state$projectionCount,
    backtrackCount = answer$state$backtrackCount,
    noMoveCount = answer$state$noMoveCount,
    postFreezeProjectionCount = answer$state$postFreezeProjectionCount,
    postFreezeBacktrackCount = answer$state$postFreezeBacktrackCount,
    postFreezeNoMoveCount = answer$state$postFreezeNoMoveCount,
    rwProposalFreezeIteration = .cop$rwProposalFreezeIteration %||% NA_integer_,
    rwBlockSizeFrozen = .cop$rwBlockSizeFrozen %||% NA_integer_,
    scoreTheory = theory)
  .cop$trace[[length(.cop$trace) + 1L]] <- list(
    kiter = as.integer(kiter), gamma = as.numeric(gain),
    beta = answer$beta,
    margin = unlist(lapply(answer$margins, `[[`, "parameters"),
      use.names = TRUE),
    copula = unlist(lapply(copulaPadFlat(answer$vine, .cop$d),
      `[[`, "parameters"), use.names = FALSE),
    value = answer$value, scoreMax = answer$scoreMax,
    residual = answer$residual,
    batchResidualMle = response$batchResidualMle %||% NA_real_,
    backend = answer$backend, final = isTRUE(final))
  invisible(NULL)
}
