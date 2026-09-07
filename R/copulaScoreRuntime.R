## Runtime integration of the pure score step with saemix.

copulaScoreBatchUpdate <- function(eta, nchains, phi = NULL, design = NULL,
                                   locationMap = NULL, beta = NULL,
                                   betaFree = NULL, transform = NULL,
                                   subject = NULL) {
  eta <- as.matrix(eta)
  if (!nrow(eta) || any(!is.finite(eta)))
    stop("score batch requires finite eta rows")
  nchains <- as.integer(nchains)
  if (length(nchains) != 1L || is.na(nchains) || nchains < 1L)
    stop("score batch requires a positive chain count")
  ## The subject layout is fixed for a run but this is called once per batch
  ## per iteration, so the checks are cached against the vector itself. They
  ## are also written to avoid unique() and match(), which were together the
  ## single largest self-time in a profile of the fit: subject labels are
  ## positive integers, so tabulate() counts them directly.
  if (is.null(subject) || length(subject) != nrow(eta))
    stop("score batch requires one positive integer subject label per row")
  if (!identical(.cop$subjectChecked, subject)) {
    if (any(!is.finite(subject)) || any(subject != as.integer(subject)) ||
        any(subject < 1L))
      stop("score batch requires one positive integer subject label per row")
    counts <- tabulate(as.integer(subject))
    if (any(counts[counts > 0L] != nchains))
      stop("score batch requires every subject exactly once per chain")
    .cop$subjectChecked <- subject
  }
  subject <- as.integer(subject)
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
    .cop$transform <- if (is.null(transform)) rep(0L, ncol(phi)) else
      rep_len(as.integer(transform), ncol(phi))
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
    .cop$transform <- NULL
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
  ## The parameter-scale score holds the natural parameter psi fixed, so the
  ## population score is the ordinary complete-data score d/dtheta log p(psi;
  ## theta): the response term does not enter because psi does not move with
  ## theta. This is what Fisher's identity needs for every family in the
  ## registry, whose supports do not depend on their parameters.
  if (identical(.cop$populationScale, "parameter") && !hasDesign)
    stop("parameter-scale score fitting requires a complete design-aware eta block")
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
  ## Which latent row each observation belongs to. The aggregate score never
  ## needed it -- it sums over observations and divides -- but a per-subject
  ## score does, because the residual block has to be attributed back to the
  ## row whose random effects produced the prediction.
  list(y = y, f = predictions, etype = x$ytype, pres = residual,
    row = as.integer(id),
    free = copulaScoreResidualIndices(errorModel),
    batchResidualMle = if (length(errorModel) == 1L &&
      errorModel %in% c("constant", "exponential", "proportional"))
      sqrt(residualSum / (nchains * nobs)) else NA_real_)
}

## Apply one score update to the current controlled-MCMC batch.
copulaScoreMstep <- function(kiter, final = FALSE, response = NULL,
                             total = NA_integer_, explore = NA_integer_) {
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
  ## Gain schedule. Single phase by default: it decays from the first
  ## iteration, which is the safe choice for a gradient step whose scale is not
  ## known in advance.
  ##
  ## The two-phase schedule of Kuhn and Lavielle (2004) holds the gain constant
  ## while exploring and decays afterwards, the decaying phase being the one
  ## the convergence conditions govern. Their constant is one, but that is a
  ## gain on SAEM's averaging of sufficient statistics, where one means "no
  ## memory"; it is not a step size. Taking it literally here -- a full step
  ## against an information matrix estimated from a handful of iterations --
  ## broke six fits in eight and left the survivors no better. What is kept is
  ## the shape rather than the constant: the initial gain, held rather than
  ## decayed, for the exploration phase only.
  ## The learning step of Baey et al. (2023), section 3.4.1, reproduced as
  ## published:
  ##
  ##   pre-heating   gamma_k = gamma_0^(1 - k/K_pre),  rising from gamma_0 to 1
  ##   heating       gamma_k = 1
  ##   decreasing    gamma_k = (k - K_heat)^-alpha
  ##
  ## with their proposed gamma_0 = 1e-4, K_pre = 1000 and alpha = 2/3. The
  ## constant gain of one is meaningful because the step is I_k^{-1} v_k, a
  ## Newton step, whose natural length is one.
  ##
  ## The end of the heating phase is theirs too, and adaptive: "averaging the
  ## norms of the gradients calculated with a third order filter of constant
  ## 1/1000 and to stop the heating phase when the norm of the averaged
  ## gradient does not decrease anymore".
  ## How long to ramp for is the one place where their design does not
  ## transfer, and the adaptation is ours rather than theirs.
  ##
  ## Their algorithm has no iteration budget. `estim()` loops over
  ## `itertools.count()`: it runs until the heating phase ends, then until
  ## `step_mean @ grad_mean * 1000 <= 1e-6`, and only then takes a fixed
  ## smoothing window of five or ten thousand iterations. The run length is an
  ## output. Their simulation study uses `smart_start = 2000` with
  ## `N_smooth = 5000`, on top of a thousand MCMC-only steps, so a run is
  ## eight to thirteen thousand iterations and the ramp is about a fifth of it.
  ##
  ## saemix has a fixed `nbiter.saemix`, so a fraction has to be chosen, and
  ## theirs is the sensible one to copy: a fifth, capped at the 2000 they
  ## actually use. Leaving it at a flat 1000, as the paper's text proposes,
  ## spends two thirds of a 1500-iteration run on a ramp whose gain is
  ## negligible until its last fifth -- measured, the difference between 53 and
  ## 96 per cent coverage, and between 1.05 and 0.27 standard deviations of
  ## bias.
  ## A re-initialisation (see copulaScoreSa.R) restarts the schedule: the
  ## iteration count and the run length are taken relative to the start of
  ## the current epoch, exactly as the step-size index is reset in Algorithm 2
  ## of Fort et al. (2016). The terminal pass keeps its absolute position.
  epoch <- .cop$scoreEpochStart %||% 0L
  kEff <- kiter - epoch
  totalEff <- if (is.finite(total)) total - epoch else total
  preheat <- max(1L, as.integer(.cop$scorePreheat %||%
    (if (is.finite(totalEff)) min(2000L, max(50L, as.integer(totalEff %/% 5L)))
     else 1000L)))
  start <- .cop$scoreGainStart %||% 1e-4
  alpha <- .cop$scoreGainPower %||% 0.8
  heatEnd <- .cop$scoreHeatEnd %||% NA_integer_
  gain <- if (kEff <= preheat) start^(1 - kEff / preheat) else
    if (!is.finite(heatEnd)) 1 else (kiter - heatEnd)^(-alpha)
  preheating <- kEff <= preheat
  ## Terminal averaging pass. The parameter step is frozen and the per-subject
  ## mean scores are averaged with gain 1/k from the start of the pass, so the
  ## information reported at the end is the outer product of scores averaged
  ## over the whole pass at one fixed parameter value.
  terminalLength <- .cop$scoreTerminal %||% 0L
  terminalStart <- if (is.finite(total) && terminalLength > 0L)
    total - terminalLength + 1L else Inf
  terminal <- kiter >= terminalStart
  deltaGain <- if (terminal) 1 / (kiter - terminalStart + 1) else gain
  ## Their reference implementation averages theta over the smoothing phase,
  ## which the algorithm box in the paper does not show: `estim()` returns the
  ## mean of theta over its final ten thousand iterations, alongside the mean
  ## of the information over the same window. Averaging therefore begins when
  ## their convergence criterion is met, not at a fixed iteration.
  averageFrom <- if (isTRUE(.cop$scoreSmoothing)) 0L else Inf
  if (terminal && !isTRUE(.cop$scoreTerminalStarted)) {
    ## Move to the Polyak average before freezing: the terminal pass has to
    ## average scores at the estimate that is reported, not at the last
    ## iterate.
    st <- .cop$scoreState
    if (!is.null(st) && isTRUE(st$averagingStarted) && !is.null(st$average)) {
      layoutNow <- copulaGaussianFremScoreLayout(.cop$margins, .cop$vine,
        .cop$d, .cop$dEta, design, locationMap, betaStart,
        if (hasDesign) .cop$betaFree else NULL, TRUE, response)
      moved <- copulaGaussianFremScoreMaterialize(st$average, layoutNow)
      .cop$margins <- moved$margins; .cop$vine <- moved$vine
      if (hasDesign) { .cop$betaCurrent <- moved$beta; betaStart <- moved$beta }
      if (!is.null(moved$residual)) {
        .cop$residualJoint <- moved$residual
        if (!is.null(response)) response$pres <- moved$residual
      }
      st$internal <- st$average
      copulaAssertFrozen(.cop$vine)
    }
    ## Restart the per-subject averages from scratch at the frozen value,
    ## whether or not a Polyak average was available to move to.
    if (!is.null(st)) {
      st$deltaSubject <- NULL; st$deltaSubjectA <- NULL; st$deltaSubjectB <- NULL
      st$fisherCount <- 0L
      .cop$scoreState <- st
    }
    .cop$scoreTerminalStarted <- TRUE
  }
  answer <- copulaGaussianFremPopulationScoreStep(
    values, rep(1 / nrow(values), nrow(values)), .cop$margins, .cop$vine,
    .cop$d, .cop$dEta, gain,
    X = design, locMap = locationMap, beta0 = betaStart,
    betaFree = if (hasDesign) .cop$betaFree else NULL,
    state = .cop$scoreState, scoreScale = .cop$scoreScale,
    finiteDifference = .cop$scoreFiniteDifference,
    projection = .cop$scoreProjection,
    adaptMetric = kiter <= scoreBurn,
    average = kiter > averageFrom && !terminal,
    preheating = preheating,
    smoothing = isTRUE(.cop$scoreSmoothing) && !terminal,
    useAverage = isTRUE(final), response = response,
    freeze = terminal, deltaGain = deltaGain,
    categoricalUniform = .cop$curCategoricalUniform,
    populationScale = .cop$populationScale,
    transform = .cop$transform,
    subject = .cop$subjectChecked)

  ## The end of heating, and the start of smoothing, exactly as their
  ## reference implementation computes them. Three details differ from the
  ## prose in the paper and all three matter.
  ##
  ## The filter runs on the *preconditioned step* `I^-1 g`, not on the
  ## gradient norm. It is initialised at zero and bias corrected, Adam style,
  ## by dividing by an accumulating `mone`; initialising it at the current
  ## value instead gives it nothing to descend from and it stops immediately.
  ## Its time constant is 100, not 1000. Their code:
  ##
  ##   factor = -expm1(-1/tc);  mone += factor*(1-mone)
  ##   m1 += factor*(val - m1); m2 += factor*(m1/mone - m2)
  ##   m3 += factor*(m2/mone - m3);  unbiased = m3/mone
  ##
  ## and heating ends when the squared norm of the unbiased third moment stops
  ## falling.
  ##
  ## The comparison is only meaningful once the filter has forgotten its
  ## initial state. A cascade of three first-order filters with time constant
  ## `tc` has an impulse response peaking at 2 tc and settling by about 3 tc;
  ## before that the "unbiased" third moment is a weighted average of a handful
  ## of steps, and two successive values differ by sampling noise, not by a
  ## trend. Without a floor the criterion fired within three iterations of
  ## heating in 115 of the 500 Gaussian-arm fits of the joint study and in 18
  ## of the 820 primary-study fits (measured from their traces on 2026-09-05),
  ## and a large-dataset run started twice the generating clearance stalled at
  ## 4.2 against 3.5 because the decreasing gain arrived while the score was
  ## still far from zero. The floor `3 tc` costs nothing when heating would
  ## have lasted longer, which was the case in every remaining fit (the
  ## shortest of them ran 380 iterations). Convergence theory (Fort et al.
  ## 2016) constrains only the decreasing phase, so a longer heating phase is
  ## always admissible; a shorter one is what stalls.
  ##
  ## Heating is also ended, at the latest, when the averaging floor below
  ## begins, so that the Polyak average is always taken over decreasing-gain
  ## iterates rather than over a constant-gain random walk.
  restarted <- FALSE
  if (isTRUE(answer$state$restartRequested)) {
    answer$state$restartRequested <- FALSE
    restarted <- TRUE
    .cop$scoreEpochStart <- as.integer(kiter)
    .cop$scoreRestarts <- (.cop$scoreRestarts %||% 0L) + 1L
    .cop$chainResetRequested <- TRUE
    .cop$scoreHeatEnd <- NA_integer_; .cop$scoreFilter <- NULL
    .cop$scoreSmoothing <- FALSE; .cop$scoreSmoothingStart <- NULL
    .cop$scoreSmoothingSource <- NULL
    .cop$scoreStepMean <- NULL; .cop$scoreGradMean <- NULL
  }
  if (!restarted && kEff > preheat &&
      !is.finite(.cop$scoreHeatEnd %||% NA_integer_)) {
    step <- as.numeric(answer$preconditionedStep)
    tc <- .cop$scoreFilterTime %||% 100
    heatMin <- as.integer(ceiling(3 * tc))
    f <- -expm1(-1 / tc)
    st <- .cop$scoreFilter %||% list(mone = 0, m1 = 0 * step, m2 = 0 * step,
      m3 = 0 * step, previous = NULL)
    st$mone <- st$mone + f * (1 - st$mone)
    st$m1 <- st$m1 + f * (step - st$m1)
    st$m2 <- st$m2 + f * (st$m1 / st$mone - st$m2)
    st$m3 <- st$m3 + f * (st$m2 / st$mone - st$m3)
    unbiased <- st$m3 / st$mone
    if (!is.null(st$previous) && kEff - preheat > heatMin &&
        sum(unbiased^2) > sum(st$previous^2))
      .cop$scoreHeatEnd <- as.integer(kiter)
    st$previous <- unbiased
    .cop$scoreFilter <- st
  }
  if (!restarted && !is.finite(.cop$scoreHeatEnd %||% NA_integer_) &&
      is.finite(totalEff) && totalEff > 0 &&
      kEff >= as.integer(ceiling(.6 * totalEff)))
    .cop$scoreHeatEnd <- as.integer(kiter)
  ## Their convergence criterion, which gates when the smoothing phase begins:
  ## `estim()` drops iterations while `end_heating is None` or
  ## `step_mean @ grad_mean * 1000 > 1e-6`, with both quantities formed as
  ## `x = factor * (value - x)` -- written that way in their code, not as a
  ## running mean -- and then averages theta and the information over the rest.
  if (is.finite(.cop$scoreHeatEnd %||% NA_integer_)) {
    sm <- .cop$scoreStepMean %||% rep(0, length(answer$score))
    gm <- .cop$scoreGradMean %||% rep(0, length(answer$score))
    .cop$scoreStepMean <- gain * (as.numeric(answer$preconditionedStep) - sm)
    .cop$scoreGradMean <- gain * (as.numeric(answer$score) - gm)
    if (!isTRUE(.cop$scoreSmoothing) &&
        sum(.cop$scoreStepMean * .cop$scoreGradMean) * 1000 <= 1e-6) {
      .cop$scoreSmoothing <- TRUE
      .cop$scoreSmoothingStart <- as.integer(kiter)
      .cop$scoreSmoothingSource <- "criterion"
    }
  }
  ## A floor under the averaging phase.
  ##
  ## The criterion above is theirs and is the right test, but it is a test, and
  ## on these models it was observed never to pass: across every replicate of
  ## several studies the information was averaged over exactly zero iterations,
  ## so what got reported was the last single Istar. That is the case this
  ## file's own comment warns about -- the last matrix "leaves all of Delta's
  ## sampling noise in it, which inflates the information and shrinks every
  ## interval" -- and it is why reported standard errors came out below the
  ## complete-data bound, which is impossible.
  ##
  ## Their reference implementation does not rely on the criterion alone
  ## either: `estim()` waits for it and then always takes a further fixed block
  ## of iterations to average over. Averaging over the last two fifths of the
  ## run reproduces that, and it can only help -- if the criterion fires first
  ## it still wins, and if it never fires there is an average instead of a
  ## single draw. It affects the reported information only; the step continues
  ## to use Istar, so the iterates are unchanged.
  if (!restarted && !isTRUE(.cop$scoreSmoothing) && is.finite(totalEff) &&
      totalEff > 0 && kEff >= as.integer(ceiling(.6 * totalEff))) {
    .cop$scoreSmoothing <- TRUE
    .cop$scoreSmoothingStart <- as.integer(kiter)
    .cop$scoreSmoothingSource <- "floor"
  }
  .cop$scoreState <- answer$state
  .cop$fisherCovariance <- answer$fisherCovariance
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
  theory$scoreScaleSource <- .cop$scoreScaleSource %||% "numeric-override"
  theory$scoreMetric <- paste0("Delattre-Kuhn empirical score-",
    "information inverse; adapted only through iteration ", scoreBurn)
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

  theory$restartCount <- .cop$scoreRestarts %||% 0L
  .cop$lastJoint <- list(kiter = kiter, conv = 0L,
    restartCount = .cop$scoreRestarts %||% 0L,
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
    metricUpdateCount = answer$metricUpdateCount,
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
