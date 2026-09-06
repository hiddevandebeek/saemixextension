## Automatic natural-parameter family screening from full posterior draws of a
## standard transformed-additive incumbent. Families are support-based and do
## not depend on the saemix coordinate transform.

## Evaluate the candidate grid. Candidate fitting is a deterministic BFGS
## optimisation over pre-generated posterior pools with no random-number use,
## so the result does not depend on evaluation order or on how the work is
## distributed. When mirai daemons are already running the grid is dispatched
## across them; otherwise evaluation is serial and identical.

copulaMapCandidates <- function(n, fun) {
  n <- as.integer(n)
  if (n <= 1L || !copulaMiraiDaemons()) return(lapply(seq_len(n), fun))
  results <- mirai::mirai_map(seq_len(n), fun)[.progress = FALSE]
  if (length(results) != n || any(vapply(results, function(x)
      inherits(x, "miraiError") || inherits(x, "errorValue"), logical(1))))
    return(lapply(seq_len(n), fun))
  results
}

copulaMiraiDaemons <- function() {
  if (!requireNamespace("mirai", quietly = TRUE)) return(FALSE)
  status <- try(mirai::status(), silent = TRUE)
  if (inherits(status, "try-error")) return(FALSE)
  connections <- suppressWarnings(as.integer(status$connections)[1L])
  isTRUE(!is.na(connections) && connections > 0L)
}

## Per-coordinate candidate families, validated against the declared support.
copulaNaturalFamilyOptions <- function(candidates, supports) {
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
  candidates
}

## The product grid enumerates K^d combinations, which becomes unusable beyond
## a couple of eta coordinates. `copulaSelectParameterMargins(search =
## "coordinate")` visits K*d per sweep instead; see there.
copulaNaturalFamilyGrid <- function(candidates, supports, maxModels) {
  candidates <- copulaNaturalFamilyOptions(candidates, supports)
  if (prod(lengths(candidates)) > maxModels)
    stop("natural-margin candidate grid exceeds maxModels")
  grid <- expand.grid(candidates, stringsAsFactors = FALSE)
  names(grid) <- paste0("parameter", seq_len(length(supports))); grid
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
  ## The baseline is the density of the drawn eta under the incumbent itself,
  ## so it must match the incumbent's own population scale. A standard
  ## transformed-additive incumbent uses the Gaussian FREM prior; a fitted
  ## flexible model already carries natural-scale margins and needs the natural
  ## prior. Dispatching here is what allows the margin to be re-estimated from
  ## the posterior of a flexible fit, not only of a Normal/lognormal one.
  baseline <- if (identical(incumbent$populationScale, "parameter"))
    copulaNaturalFremLogPrior(data$eta,
      if (hasConditioning) data$conditioning else
        matrix(numeric(), nrow(data$eta), 0L),
      incumbent$vine, incumbent$margins, incumbent$dEta,
      data$predictor, data$transform, "joint") else
    copulaGaussianFremLogPrior(baselineState, incumbent$vine,
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
  data$columnCache <- new.env(parent = emptyenv())
  data
}

copulaNaturalBridgeCandidateLogPrior <- function(margins, data, incumbent) {
  cache <- data$bridgeCache
  if (is.null(cache)) data <- copulaNaturalBridgePrepare(data, incumbent)
  copulaNaturalBridgeState(margins, data)$logCandidate
}

## Shared state behind both the bridge objective and its gradient: the margin
## evaluation, the valid rows, the copula product Qz on those rows, and the
## candidate log prior. The optimiser asks for value and gradient at the same
## point, so computing this once and reusing it removes a full margin sweep
## from every gradient call.
## One column of the margin evaluation, reused when that coordinate's margin
## has not moved. Both callers vary one coordinate at a time -- the ladder
## raises a single degree, the bridge gradient perturbs a single parameter --
## so every other column repeats work on a pool of subjects times draws. One
## entry per column is enough to catch that, and the key is compared exactly,
## so a miss only costs a recomputation.
copulaNaturalBridgeColumns <- function(margins, data) {
  store <- data$columnCache
  if (is.null(store))
    return(copulaNaturalMarginsEvaluate(data$natural, data$typical, margins))
  rows <- nrow(data$natural)
  z <- logMargin <- matrix(NA_real_, rows, length(margins))
  valid <- rep(TRUE, rows)
  for (j in seq_along(margins)) {
    margin <- margins[[j]]
    key <- c(margin$name, margin$support,
      format(margin$parameters, digits = 17))
    slot <- as.character(j)
    entry <- store[[slot]]
    if (is.null(entry) || length(entry$key) != length(key) ||
        any(entry$key != key)) {
      entry <- copulaMarginColumn(margin, data$natural[, j], data$typical[, j])
      entry$key <- key
      assign(slot, entry, envir = store)
    }
    valid <- valid & entry$valid
    z[entry$valid, j] <- entry$z[entry$valid]
    logMargin[entry$valid, j] <- entry$logMargin[entry$valid]
  }
  list(z = z, logMargin = logMargin, valid = valid)
}

copulaNaturalBridgeState <- function(margins, data) {
  cache <- data$bridgeCache
  evaluated <- copulaNaturalBridgeColumns(margins, data)
  valid <- evaluated$valid & cache$covariate$valid
  rows <- which(valid)
  logCandidate <- rep(-Inf, nrow(data$eta))
  if (!length(rows))
    return(list(evaluated = evaluated, rows = rows, Qzr = NULL,
      logCandidate = logCandidate))
  z <- cbind(evaluated$z, cache$covariate$z)
  zr <- z[rows, , drop = FALSE]
  Qzr <- zr %*% cache$copulaQuadratic
  logCandidate[rows] <- -.5 * cache$logDet - .5 * rowSums(Qzr * zr) +
    rowSums(evaluated$logMargin[rows, , drop = FALSE]) +
    rowSums(cache$covariate$logMargin[rows, , drop = FALSE]) +
    cache$logJacobian[rows]
  list(evaluated = evaluated, rows = rows, Qzr = Qzr,
    logCandidate = logCandidate)
}

## Gradient of the bridge objective in internal coordinates.
##
## deltaLogLik = sum_i log mean_r exp(w_ir) with w = candidate - baseline, and
## the candidate log prior is
##   -0.5 logdet - 0.5 z' Q z + sum_j log f_j + (terms free of the candidate).
## Differentiating gives
##   d deltaLogLik / d theta = sum_ir p_ir { d log f_j - (Qz)_j d xi_j }
## with p_ir the per-subject softmax of w. The chain rule through the copula is
## therefore exact; only d log f_j and d xi_j are differenced, and a single
## internal coordinate perturbs exactly one margin. That replaces 2p full
## objective evaluations -- each of which re-evaluates every margin, the
## quadratic form and the log-sum-exp -- with 2p single-margin evaluations.

copulaNaturalBridgeGradient <- function(margins, data, incumbent, materialize,
                                        internal, step = 1e-5, state = NULL) {
  cache <- data$bridgeCache
  dEta <- length(margins)
  if (is.null(state)) state <- copulaNaturalBridgeState(margins, data)
  rows <- state$rows
  if (!length(rows)) return(NULL)
  Qzr <- state$Qzr
  w <- matrix(state$logCandidate - cache$baseline,
    nrow = data$n, ncol = data$samples)
  maximum <- apply(w, 1L, max)
  if (any(!is.finite(maximum))) return(NULL)
  shifted <- exp(w - maximum)
  weight <- as.numeric(shifted / rowSums(shifted))[rows]

  gradient <- numeric(length(internal))
  for (m in seq_along(internal)) {
    plus <- minus <- internal
    plus[m] <- plus[m] + step; minus[m] <- minus[m] - step
    marginsPlus <- try(materialize(plus), silent = TRUE)
    marginsMinus <- try(materialize(minus), silent = TRUE)
    if (inherits(marginsPlus, "try-error") ||
        inherits(marginsMinus, "try-error")) return(NULL)
    changed <- which(vapply(seq_len(dEta), function(j)
      !identical(marginsPlus[[j]]$parameters, marginsMinus[[j]]$parameters),
      logical(1)))
    for (j in changed) {
      natural <- data$natural[, j]; typical <- data$typical[, j]
      up <- copulaMarginColumn(marginsPlus[[j]], natural, typical)
      down <- copulaMarginColumn(marginsMinus[[j]], natural, typical)
      if (!all(up$valid[rows]) || !all(down$valid[rows])) return(NULL)
      dScore <- (up$z[rows] - down$z[rows]) / (2 * step)
      dDensity <- (up$logMargin[rows] - down$logMargin[rows]) / (2 * step)
      gradient[m] <- gradient[m] +
        sum(weight * (dDensity - Qzr[, j] * dScore))
    }
  }
  if (any(!is.finite(gradient))) return(NULL)
  gradient
}

copulaNaturalBridgeMetric <- function(margins, data, incumbent) {
  if (is.null(data$bridgeCache)) data <- copulaNaturalBridgePrepare(data, incumbent)
  candidate <- copulaNaturalBridgeCandidateLogPrior(margins, data, incumbent)
  baseline <- data$bridgeCache$baseline
  copulaPosteriorBridgeMetricFromLogWeight(matrix(candidate - baseline,
    nrow = data$n, ncol = data$samples))
}

## Candidate margin parameters are fitted by an unconstrained quasi-Newton
## search over the internal coordinates, both methods supplied by stats::optim
## so the choice adds no dependency. L-BFGS-B is the default: over a five-grid
## benchmark it reached the same optimum as BFGS to within 1e-8 on every
## candidate while using 170 objective evaluations against 196. Tightening its
## convergence factor is counterproductive -- factr 1e5 needed 235 evaluations
## and 1e1 needed 1175, for no gain in the attained optimum -- so the stock
## tolerance is kept. Gradient-based nloptr methods were slower here, because
## a finite-difference gradient at R level costs more than the internal one.

copulaCandidateOptimise <- function(start, objective, maxit,
                                    optimiser = c("L-BFGS-B", "BFGS"),
                                    gradient = NULL) {
  optimiser <- match.arg(optimiser)
  control <- if (identical(optimiser, "BFGS"))
    list(maxit = as.integer(maxit), reltol = 1e-7) else
    list(maxit = as.integer(maxit))
  stats::optim(start, objective, gr = gradient, method = optimiser,
    control = control)
}

copulaFitNaturalCandidate <- function(families, supports, train, validation,
                                      incumbent, maxit, optimiser = "L-BFGS-B") {
  start <- lapply(seq_along(families), function(j)
    copulaNaturalMarginStart(families[j], train$natural[, j],
      train$typical[, j], supports[j]))
  layout <- copulaMarginLayout(start)
  internal <- copulaScoreToInternal(layout$par, layout$lower, layout$upper)
  materialize <- function(value) copulaMarginsWithParameters(start, layout,
    copulaScoreFromInternal(value, layout$lower, layout$upper))
  ## The optimiser evaluates value and gradient at the same point, so the
  ## shared bridge state is kept from the objective and reused by the gradient.
  memo <- new.env(parent = emptyenv()); memo$value <- NULL
  bridgeState <- function(value, margins) {
    if (!is.null(memo$value) && length(memo$value) == length(value) &&
        all(memo$value == value)) return(memo$state)
    state <- copulaNaturalBridgeState(margins, train)
    memo$value <- value; memo$state <- state
    state
  }
  objective <- function(value) {
    margins <- try(materialize(value), silent = TRUE)
    if (inherits(margins, "try-error")) return(.Machine$double.xmax / 100)
    state <- try(bridgeState(value, margins), silent = TRUE)
    if (inherits(state, "try-error")) return(.Machine$double.xmax / 100)
    metric <- try(copulaPosteriorBridgeMetricFromLogWeight(matrix(
      state$logCandidate - train$bridgeCache$baseline,
      nrow = train$n, ncol = train$samples)), silent = TRUE)
    if (inherits(metric, "try-error") || !is.finite(metric$deltaLogLik))
      return(.Machine$double.xmax / 100)
    -metric$deltaLogLik
  }
  ## Supplied gradient, with the difference quotient of the objective as a
  ## fallback: the analytic form needs every posterior draw to stay inside the
  ## candidate support, which a far-off trial point can violate.
  numericalGradient <- function(value, h = 1e-5) vapply(seq_along(value),
    function(k) {
      plus <- minus <- value; plus[k] <- plus[k] + h; minus[k] <- minus[k] - h
      (objective(plus) - objective(minus)) / (2 * h)
    }, numeric(1))
  gradient <- function(value) {
    margins <- try(materialize(value), silent = TRUE)
    if (inherits(margins, "try-error")) return(numericalGradient(value))
    state <- try(bridgeState(value, margins), silent = TRUE)
    if (inherits(state, "try-error")) return(numericalGradient(value))
    exact <- try(copulaNaturalBridgeGradient(margins, train, incumbent,
      materialize, value, state = state), silent = TRUE)
    if (inherits(exact, "try-error") || is.null(exact))
      return(numericalGradient(value))
    -exact
  }
  fit <- copulaCandidateOptimise(internal, objective, maxit, optimiser,
    gradient)
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
    maxModels = 81L, optimizerMaxit = 100L, optimiser = "L-BFGS-B",
    search = c("coordinate", "product"), maxSweeps = 3L,
    minimumEssFraction = .01, minimumPosteriorEssFraction = .05,
    maximumMcse = .5, maximumParetoK = Inf) {
  search <- match.arg(search)
  if (!inherits(object, "SaemixObject")) stop("object must be a fitted SaemixObject")
  state <- copulaGet(object)
  if (!state$populationScale %in% c("transformed-additive", "parameter") ||
      !copulaIsFullGaussianVine(state$vine, state$d))
    stop("natural-margin screening requires a Gaussian-copula incumbent")
  if ((state$dConditioning %||% 0L) > 0L &&
      (anyNA(state$conditioning) || any(vapply(
        state$margins[state$dEta + seq_len(state$dConditioning)],
        function(m) !identical(m$type, "continuous"), logical(1)))))
    stop("natural-margin screening with covariates requires complete continuous conditioning")
  supports <- rep_len(as.character(supports), state$dEta)
  if (any(!supports %in% c("real", "positive", "unit")))
    stop("supports must be real, positive, or unit")
  options <- copulaNaturalFamilyOptions(candidates, supports)
  trainingDraws <- copulaPosteriorEtaDraws(object, posteriorDraws, max.iter,
    as.integer(seed))
  validationDraws <- copulaPosteriorEtaDraws(object, posteriorDraws, max.iter,
    as.integer(seed) + 1L)
  training <- copulaNaturalPosteriorData(object, trainingDraws)
  validation <- copulaNaturalPosteriorData(object, validationDraws)
  training <- copulaNaturalBridgePrepare(training, state)
  validation <- copulaNaturalBridgePrepare(validation, state)
  evaluateCombination <- function(families)
    copulaFitNaturalCandidate(families, supports, training, validation, state,
      optimizerMaxit, optimiser)

  if (identical(search, "product")) {
    ## Every combination of families: K^d evaluations.
    grid <- copulaNaturalFamilyGrid(candidates, supports, as.integer(maxModels))
    fitted <- copulaMapCandidates(nrow(grid), function(row)
      evaluateCombination(as.character(grid[row, ])))
  } else {
    ## Coordinate ascent: vary one eta coordinate at a time with the others
    ## held at their current family, K*d evaluations per sweep, repeating until
    ## no coordinate moves. Each comparison is still the exact likelihood-ratio
    ## identity, so this changes which combinations are visited, never how any
    ## of them is scored, and the retained optimum cannot be worse than the
    ## starting combination. Against the full product grid on twelve
    ## replicates the two agreed every time, with a log-likelihood difference
    ## of exactly zero.
    store <- new.env(parent = emptyenv())
    key <- function(families) paste(families, collapse = "|")
    current <- vapply(options, `[`, character(1), 1L)
    for (sweep in seq_len(max(1L, as.integer(maxSweeps)))) {
      moved <- FALSE
      for (j in seq_along(options)) {
        trials <- lapply(options[[j]], function(family) {
          trial <- current; trial[j] <- family; trial })
        pending <- which(vapply(trials,
          function(t) is.null(store[[key(t)]]), logical(1)))
        if (length(pending)) {
          got <- copulaMapCandidates(length(pending), function(k)
            evaluateCombination(trials[[pending[k]]]))
          for (k in seq_along(pending))
            assign(key(trials[[pending[k]]]), got[[k]], envir = store)
        }
        scores <- vapply(trials, function(t) {
          got <- store[[key(t)]]
          if (is.null(got) || !is.finite(got$validationBicAdvantage)) -Inf else
            got$validationBicAdvantage
        }, numeric(1))
        best <- options[[j]][which.max(scores)]
        if (!identical(best, current[j])) { current[j] <- best; moved <- TRUE }
      }
      if (!moved) break
    }
    fitted <- unname(as.list(store))
  }
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
      ## Reported, not gated by default. The classical guidance is that k
      ## above 0.7 makes an importance estimator unreliable, but the shipped
      ## positive-support candidates sit near that value in cases the method
      ## demonstrably gets right, and k is a lagging indicator of a tail that
      ## has not been sampled, so it is surfaced for the caller to judge.
      maximum_pareto_k = result$validation$maximumParetoK %||% NA_real_,
      eligible = result$convergence == 0L &&
        is.finite(result$validationBicAdvantage) &&
        result$validation$minimumEssFraction >= minimumEssFraction &&
        result$validation$mcse <= maximumMcse &&
        (!is.finite(result$validation$maximumParetoK %||% NA_real_) ||
          result$validation$maximumParetoK <= maximumParetoK))
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
