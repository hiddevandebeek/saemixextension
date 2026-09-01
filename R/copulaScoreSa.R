## Likelihood-score stochastic-approximation backend for a fixed
## Gaussian-copula FREM model. It updates one coherent parameter vector:
## population-location coefficients, every free native marginal parameter,
## every Gaussian partial correlation, and supported residual-error
## parameters. Normal, lognormal, Weibull, Laplace, Gamma-scale, location,
## dependence-angle, and residual components are analytical. Gamma-shape and
## declared custom-margin components use centered differences with a shrinking
## step covered by the paper's summable random-perturbation assumption.

copulaScoreBounds <- function(value, lower, upper) {
  value <- as.numeric(value)
  lower <- rep_len(as.numeric(lower), length(value))
  upper <- rep_len(as.numeric(upper), length(value))
  both <- is.finite(lower) & is.finite(upper)
  onlyLower <- is.finite(lower) & !is.finite(upper)
  onlyUpper <- !is.finite(lower) & is.finite(upper)
  free <- !is.finite(lower) & !is.finite(upper)
  wideLocation <- both & lower < 0 & upper > 0 & (upper - lower) > 1e6
  both[wideLocation] <- FALSE
  free <- free | wideLocation
  list(value = value, lower = lower, upper = upper, both = both,
    onlyLower = onlyLower, onlyUpper = onlyUpper, free = free)
}

copulaScoreToInternal <- function(value, lower, upper, epsilon = 1e-10) {
  bounds <- copulaScoreBounds(value, lower, upper)
  with(bounds, {
    answer <- numeric(length(value))
    if (any(both)) {
      fraction <- (value[both] - lower[both]) / (upper[both] - lower[both])
      answer[both] <- stats::qlogis(pmin(1 - epsilon,
        pmax(epsilon, fraction)))
    }
    if (any(onlyLower))
      answer[onlyLower] <- log(pmax(epsilon,
        value[onlyLower] - lower[onlyLower]))
    if (any(onlyUpper))
      answer[onlyUpper] <- log(pmax(epsilon,
        upper[onlyUpper] - value[onlyUpper]))
    answer[free] <- value[free]
    answer
  })
}

copulaScoreFromInternal <- function(value, lower, upper) {
  bounds <- copulaScoreBounds(value, lower, upper)
  with(bounds, {
    answer <- numeric(length(value))
    if (any(both))
      answer[both] <- lower[both] +
        (upper[both] - lower[both]) * stats::plogis(value[both])
    if (any(onlyLower))
      answer[onlyLower] <- lower[onlyLower] + exp(value[onlyLower])
    if (any(onlyUpper))
      answer[onlyUpper] <- upper[onlyUpper] - exp(value[onlyUpper])
    answer[free] <- value[free]
    answer
  })
}

copulaScoreInternalDerivative <- function(value, lower, upper) {
  bounds <- copulaScoreBounds(value, lower, upper)
  with(bounds, {
    answer <- numeric(length(value))
    if (any(both)) {
      probability <- stats::plogis(value[both])
      answer[both] <- (upper[both] - lower[both]) *
        probability * (1 - probability)
    }
    if (any(onlyLower)) answer[onlyLower] <- exp(value[onlyLower])
    if (any(onlyUpper)) answer[onlyUpper] <- -exp(value[onlyUpper])
    answer[free] <- 1
    answer
  })
}

copulaGaussianMarginNativeScore <- function(
    x, margin, z, zInfluence, numericalStep = 1e-4) {
  parameters <- margin$parameters; names <- names(parameters)
  gaussianScore <- z - zInfluence
  if (identical(margin$name, "normal")) {
    mean <- if ("mean" %in% names) unname(parameters["mean"]) else 0
    sd <- unname(parameters["sd"])
    score <- matrix(NA_real_, length(x), length(parameters),
      dimnames = list(NULL, names))
    if ("mean" %in% names) score[, "mean"] <- gaussianScore / sd
    score[, "sd"] <- (-1 + z * gaussianScore) / sd
    return(list(parameter = score, x = -gaussianScore / sd,
      numerical = character(), oneSided = FALSE))
  }
  if (identical(margin$name, "lognormal")) {
    sdlog <- unname(parameters["sdlog"])
    score <- cbind(meanlog = gaussianScore / sdlog,
      sdlog = (-1 + z * gaussianScore) / sdlog)
    return(list(parameter = score,
      x = -1 / x - gaussianScore / (sdlog * x),
      numerical = character(), oneSided = FALSE))
  }
  normalDensity <- stats::dnorm(z)
  if (any(!is.finite(normalDensity)) || any(normalDensity <= 1e-300))
    return(NULL)
  combine <- function(dlog, dF) dlog + zInfluence * dF / normalDensity
  density <- exp(margin$log_density(x, parameters))
  score <- matrix(NA_real_, length(x), length(parameters),
    dimnames = list(NULL, names)); numerical <- character()
  scoreX <- NULL; oneSided <- FALSE

  if (identical(margin$name, "gamma")) {
    shape <- unname(parameters["shape"]); scale <- unname(parameters["scale"])
    h <- numericalStep * max(1, abs(shape))
    lo <- max(margin$lower[match("shape", names)], shape - h)
    hi <- min(margin$upper[match("shape", names)], shape + h)
    oneSided <- oneSided || abs((hi - shape) - (shape - lo)) >
      1e-10 * max(1, h)
    dFshape <- (stats::pgamma(x, hi, scale = scale) -
      stats::pgamma(x, lo, scale = scale)) / (hi - lo)
    score[, "shape"] <- combine(log(x) - log(scale) - digamma(shape),
      dFshape)
    score[, "scale"] <- combine(x / scale^2 - shape / scale,
      -x * density / scale)
    scoreX <- (shape - 1) / x - 1 / scale +
      zInfluence * density / normalDensity
    numerical <- "gamma.shape.cdf"
  } else if (identical(margin$name, "weibull")) {
    shape <- unname(parameters["shape"]); scale <- unname(parameters["scale"])
    logRatio <- log(x / scale); power <- exp(shape * logRatio)
    survival <- exp(-power)
    score[, "shape"] <- combine(1 / shape + (1 - power) * logRatio,
      survival * power * logRatio)
    score[, "scale"] <- combine(shape * (power - 1) / scale,
      -survival * shape * power / scale)
    scoreX <- (shape - 1 - shape * power) / x +
      zInfluence * density / normalDensity
  } else if (identical(margin$name, "laplace")) {
    sd <- unname(parameters["sd"]); b <- sd / sqrt(2)
    F <- margin$cdf(x, parameters)
    dFsd <- ifelse(x < 0, -x * F / (b * sd),
      -x * (1 - F) / (b * sd))
    score[, "sd"] <- combine(-1 / sd + abs(x) / (b * sd), dFsd)
    scoreX <- -sign(x) / b + zInfluence * density / normalDensity
  }

  if (is.null(scoreX) || anyNA(score)) {
    ## Generic declared continuous margins retain the summable centered-
    ## difference route for each native parameter and the native coordinate.
    numerical <- c(numerical, paste0(margin$name, ".generic"))
    baseLog <- margin$log_density(x, parameters)
    baseF <- margin$cdf(x, parameters)
    for (j in seq_along(parameters)) if (anyNA(score[, j])) {
      h <- numericalStep * max(1, abs(parameters[j]))
      lo <- max(margin$lower[j], parameters[j] - h)
      hi <- min(margin$upper[j], parameters[j] + h)
      oneSided <- oneSided || abs((hi - parameters[j]) -
        (parameters[j] - lo)) > 1e-10 * max(1, h)
      plus <- minus <- parameters; plus[j] <- hi; minus[j] <- lo
      dlog <- (margin$log_density(x, plus) -
        margin$log_density(x, minus)) / (hi - lo)
      dF <- (margin$cdf(x, plus) - margin$cdf(x, minus)) / (hi - lo)
      score[, j] <- combine(dlog, dF)
    }
    if (is.null(scoreX)) {
      h <- numericalStep * pmax(1, abs(x))
      xp <- x + h; xm <- x - h
      supportValid <- is.finite(margin$log_density(xm, parameters))
      if (any(!supportValid)) oneSided <- TRUE
      xm[!supportValid] <- x[!supportValid]
      denominator <- xp - xm
      dlog <- (margin$log_density(xp, parameters) -
        margin$log_density(xm, parameters)) / denominator
      dF <- (margin$cdf(xp, parameters) -
        margin$cdf(xm, parameters)) / denominator
      scoreX <- combine(dlog, dF)
    }
  }
  if (any(!is.finite(score)) || any(!is.finite(scoreX))) return(NULL)
  list(parameter = score, x = scoreX, numerical = unique(numerical),
    oneSided = oneSided)
}

copulaGaussianFremScoreLayout <- function(
    margins, vine, d, dEta, X = NULL, locMap = NULL,
    beta0 = NULL, betaFree = NULL, withMu = TRUE,
    response = NULL) {
  marginLayout <- copulaMarginLayout(margins)
  correlationAngles <- copulaGaussianCorrelationAngles(
    copulaGaussianRvineCor(vine, d))
  edgePar <- correlationAngles
  edgeLower <- rep(1e-6, length(edgePar))
  edgeUpper <- rep(pi - 1e-6, length(edgePar))
  hasDesign <- !is.null(X) && !is.null(locMap) && !is.null(beta0)
  if (hasDesign) {
    X <- as.matrix(X); locMap <- as.matrix(locMap); beta0 <- as.numeric(beta0)
    if (ncol(X) != nrow(locMap) || ncol(locMap) != d ||
        length(beta0) != ncol(X)) stop("invalid population-score design")
    locationIndex <- which(rowSums(abs(locMap)) > 0)
    if (!is.null(betaFree))
      locationIndex <- intersect(locationIndex, as.integer(betaFree))
    location <- beta0[locationIndex]
  } else {
    locationIndex <- integer()
    location <- if (isTRUE(withMu)) rep(0, dEta) else numeric()
  }
  residual <- residualFree <- numeric()
  if (!is.null(response)) {
    if (!is.list(response) || is.null(response$pres) ||
        is.null(response$free)) stop("invalid population-score response block")
    if (length(response$y) != length(response$f) ||
        length(response$etype) != length(response$f) ||
        any(!is.finite(response$y)) || any(!is.finite(response$f)) ||
        anyNA(response$etype))
      stop("population-score response values, predictions, and types must align")
    residualFree <- as.integer(response$free)
    residual <- as.numeric(response$pres)[residualFree]
    if (!length(residual) || any(!is.finite(residual)) || any(residual <= 0))
      stop("population-score residual parameters must be finite and positive")
  }
  native <- c(location, marginLayout$par, edgePar, residual)
  lower <- c(rep(-Inf, length(location)), marginLayout$lower, edgeLower,
    rep(1e-10, length(residual)))
  upper <- c(rep( Inf, length(location)), marginLayout$upper, edgeUpper,
    rep(Inf, length(residual)))
  list(native = native, lower = lower, upper = upper,
    nLocation = length(location), nMargin = length(marginLayout$par),
    nEdge = length(edgePar), nResidual = length(residual),
    residualFree = residualFree, response = response,
    locationIndex = locationIndex,
    marginLayout = marginLayout,
    hasDesign = hasDesign, X = X, locMap = locMap, beta0 = beta0,
    d = d, dEta = dEta, vine0 = vine, margins0 = margins)
}

copulaGaussianFremScoreMaterialize <- function(internal, layout,
                                                buildVine = TRUE) {
  native <- copulaScoreFromInternal(internal, layout$lower, layout$upper)
  nLocation <- layout$nLocation; nMargin <- layout$nMargin
  location <- if (nLocation) native[seq_len(nLocation)] else numeric()
  marginParameters <- if (nMargin)
    native[nLocation + seq_len(nMargin)] else numeric()
  edgeParameters <- if (layout$nEdge)
    native[nLocation + nMargin + seq_len(layout$nEdge)] else numeric()
  residualParameters <- if (layout$nResidual)
    native[nLocation + nMargin + layout$nEdge +
      seq_len(layout$nResidual)] else numeric()
  ## The unconstrained-to-native map already enforces every declared bound.
  ## Avoid revalidating unchanged callback objects several times per score
  ## iteration; only their small native parameter vectors are replaced here.
  margins <- layout$margins0
  for (j in seq_along(margins)) {
    index <- layout$marginLayout$index[[j]]
    if (length(index)) {
      parameters <- margins[[j]]$parameters
      parameters[margins[[j]]$free] <- marginParameters[index]
      margins[[j]]$parameters <- parameters
    }
  }
  correlation <- if (layout$nEdge)
    copulaGaussianCorrelationFromAngles(edgeParameters, layout$d)$R else
    copulaGaussianRvineCor(layout$vine0, layout$d)
  vine <- if (!isTRUE(buildVine)) NULL else if (layout$nEdge)
    copulaGaussianRvineUpdateCor(layout$vine0, correlation) else layout$vine0
  if (layout$hasDesign) {
    beta <- layout$beta0
    if (nLocation) beta[layout$locationIndex] <- location
    delta <- numeric(layout$d)
  } else {
    beta <- layout$beta0; delta <- numeric(layout$d)
    if (nLocation) delta[seq_len(layout$dEta)] <- location
  }
  residual <- if (is.null(layout$response)) NULL else
    as.numeric(layout$response$pres)
  if (length(residualParameters))
    residual[layout$residualFree] <- residualParameters
  list(native = native, beta = beta, delta = delta, residual = residual,
    margins = margins, vine = vine, correlation = correlation)
}

copulaGaussianFremCompleteScoreInternal <- function(
    internal, layout, E, w, numericalStep) {
  candidate <- try(copulaGaussianFremScoreMaterialize(internal, layout,
    buildVine = FALSE),
    silent = TRUE)
  if (inherits(candidate, "try-error")) return(NULL)
  residual <- if (layout$hasDesign)
    E - copulaLocation(layout$X, candidate$beta, layout$locMap) else
    sweep(E, 2L, candidate$delta, "-")
  evaluated <- try(copulaGaussianFremEvaluateMargins(
    residual, candidate$margins), silent = TRUE)
  if (inherits(evaluated, "try-error") || !all(evaluated$valid)) return(NULL)
  R <- candidate$correlation
  Omega <- solve(R); zInfluence <- evaluated$z %*% (diag(layout$d) - Omega)
  parameterGradient <- numeric(layout$nMargin)
  scoreX <- matrix(NA_real_, nrow(E), layout$d)
  numerical <- character(); oneSided <- FALSE
  for (j in seq_len(layout$d)) {
    local <- copulaGaussianMarginNativeScore(residual[, j],
      candidate$margins[[j]], evaluated$z[, j], zInfluence[, j],
      numericalStep)
    if (is.null(local)) return(NULL)
    index <- layout$marginLayout$index[[j]]
    if (length(index))
      parameterGradient[index] <- colSums(w *
        local$parameter[, candidate$margins[[j]]$free, drop = FALSE])
    scoreX[, j] <- local$x; numerical <- c(numerical, local$numerical)
    oneSided <- oneSided || isTRUE(local$oneSided)
  }
  locationGradient <- numeric(layout$nLocation)
  if (layout$nLocation) {
    if (layout$hasDesign) for (k in seq_along(layout$locationIndex)) {
      l <- layout$locationIndex[k]
      derivative <- sweep(matrix(layout$locMap[l, ], nrow(E), layout$d,
        byrow = TRUE), 1L, layout$X[, l], "*")
      locationGradient[k] <- -sum(w * rowSums(scoreX * derivative))
    } else locationGradient <- -colSums(w *
      scoreX[, seq_len(layout$dEta), drop = FALSE])
  }
  S <- crossprod(evaluated$z, w * evaluated$z)
  angleStart <- layout$nLocation + layout$nMargin
  angles <- if (layout$nEdge)
    candidate$native[angleStart + seq_len(layout$nEdge)] else numeric()
  correlationGradient <- if (layout$nEdge)
    -copulaGaussianCorrelationAngleObjective(
      angles, S, layout$d)$gradient else numeric()
  residualGradient <- numeric(layout$nResidual)
  if (layout$nResidual) {
    responseSd <- error(layout$response$f, candidate$residual,
      layout$response$etype)
    responseType <- rep_len(as.integer(layout$response$etype),
      length(layout$response$f))
    difference <- layout$response$y - layout$response$f
    common <- difference^2 / responseSd^3 - 1 / responseSd
    full <- numeric(length(candidate$residual))
    for (type in sort(unique(responseType))) {
      rows <- responseType == type; ia <- 2L * type - 1L; ib <- 2L * type
      a <- candidate$residual[ia]; b <- candidate$residual[ib]
      full[ia] <- sum(common[rows] * a / responseSd[rows]) / nrow(E)
      full[ib] <- sum(common[rows] * b * layout$response$f[rows]^2 /
        responseSd[rows]) / nrow(E)
    }
    residualGradient <- full[layout$residualFree]
  }
  nativeGradient <- c(locationGradient, parameterGradient,
    correlationGradient, residualGradient)
  derivative <- copulaScoreInternalDerivative(internal,
    layout$lower, layout$upper)
  gradient <- nativeGradient * derivative
  if (length(gradient) != length(internal) || any(!is.finite(gradient)))
    return(NULL)
  list(gradient = gradient, nativeGradient = nativeGradient,
    numerical = unique(numerical), candidate = candidate,
    residual = residual, evaluated = evaluated, oneSided = oneSided)
}

copulaGaussianFremPopulationScoreStep <- function(
    E, w, margins0, vine0, d, dEta, gain,
    X = NULL, locMap = NULL, beta0 = NULL, betaFree = NULL,
    state = NULL, scoreScale = 0.05, finiteDifference = 1e-4,
    projection = 24, withMu = TRUE, adaptMetric = TRUE,
    average = FALSE, useAverage = FALSE, response = NULL,
    analyticScore = TRUE, categoricalUniform = NULL,
    referenceUniform = NULL) {
  E <- as.matrix(E); w <- as.numeric(w); w <- w / sum(w)
  if (ncol(E) != d || nrow(E) != length(w) ||
      anyNA(E[, seq_len(dEta), drop = FALSE]) ||
      any(!is.finite(E[!is.na(E)])) || any(!is.finite(w)) || any(w < 0) ||
      !copulaIsFullGaussianVine(vine0, d))
    stop("population score step requires finite eta and valid Gaussian-copula rows")
  if (length(gain) != 1L || !is.finite(gain) || gain <= 0 || gain > 1 ||
      length(scoreScale) != 1L || !is.finite(scoreScale) || scoreScale <= 0 ||
      length(finiteDifference) != 1L || !is.finite(finiteDifference) ||
      finiteDifference <= 0)
    stop("invalid population score gain or numerical control")
  layout <- copulaGaussianFremScoreLayout(margins0, vine0, d, dEta,
    X, locMap, beta0, betaFree, withMu, response)
  current <- copulaScoreToInternal(layout$native, layout$lower, layout$upper)
  discrete <- which(vapply(margins0, function(m)
    identical(m$type, "discrete"), logical(1)))
  if (is.null(referenceUniform))
    referenceUniform <- matrix(NA_real_, nrow(E), d) else {
    referenceUniform <- as.matrix(referenceUniform)
    if (any(dim(referenceUniform) != c(nrow(E), d)))
      stop("score-sa fixed-reference augmentation does not align with the MCMC batch")
  }
  hasReference <- any(is.finite(referenceUniform))
  completeContinuous <- !length(discrete) && !anyNA(E)
  completeAugmented <- !anyNA(E) && (completeContinuous ||
    (!is.null(categoricalUniform) &&
      all(is.finite(as.matrix(categoricalUniform)[, discrete, drop = FALSE]))))
  if (length(discrete)) {
    categoricalUniform <- as.matrix(categoricalUniform)
    if (any(dim(categoricalUniform) != c(nrow(E), d)))
      stop("score-sa categorical augmentation does not align with the MCMC batch")
  }
  materializeState <- function(candidate) {
    location <- if (layout$hasDesign)
      copulaLocation(layout$X, candidate$beta, layout$locMap) else
      matrix(candidate$delta, nrow(E), d, byrow = TRUE)
    residual <- E - location
    for (j in seq_len(d)) {
      rows <- is.finite(referenceUniform[, j])
      if (any(rows)) residual[rows, j] <- candidate$margins[[j]]$quantile(
        referenceUniform[rows, j], candidate$margins[[j]]$parameters)
    }
    list(residual = residual, absolute = location + residual)
  }
  objectiveFailure <- NULL
  objective <- function(internal) {
    candidate <- try(copulaGaussianFremScoreMaterialize(internal, layout,
      buildVine = !completeAugmented),
      silent = TRUE)
    if (inherits(candidate, "try-error")) {
      objectiveFailure <<- as.character(candidate); return(-Inf)
    }
    completeState <- materializeState(candidate)
    residual <- completeState$residual
    value <- try({
      if (completeAugmented || hasReference) {
        evaluated <- if (completeContinuous && !hasReference)
          copulaGaussianFremEvaluateMargins(residual, candidate$margins) else
          copulaGaussianFremAugmentedEvaluateMargins(residual,
            candidate$margins, categoricalUniform, referenceUniform)
        if (!all(evaluated$valid)) stop("invalid complete-data margin")
        correlation <- candidate$correlation
        logPrior <- copulaGaussianLogDensity(evaluated$z, correlation) -
          rowSums(stats::dnorm(evaluated$z, log = TRUE)) +
          rowSums(evaluated$logMargin)
        sum(w * logPrior)
      } else {
        sum(w * copulaGaussianFremLogPrior(residual,
          candidate$vine, candidate$margins, dEta, "joint"))
      }
    }, silent = TRUE)
    if (inherits(value, "try-error") || !is.finite(value)) {
      objectiveFailure <<- if (inherits(value, "try-error"))
        as.character(value) else "non-finite population objective"
      return(-Inf)
    }
    if (!is.null(layout$response)) {
      if (is.function(layout$response$evaluate) && hasReference) {
        responseLog <- try(layout$response$evaluate(
          completeState$absolute[, seq_len(layout$dEta), drop = FALSE],
          candidate$residual), silent = TRUE)
        if (inherits(responseLog, "try-error") ||
            any(!is.finite(responseLog))) {
          objectiveFailure <<- "moving-support response evaluation failed"
          return(-Inf)
        }
        return(value + sum(responseLog) / nrow(E))
      }
      responseSd <- try(error(layout$response$f, candidate$residual,
        layout$response$etype), silent = TRUE)
      if (inherits(responseSd, "try-error") ||
          any(!is.finite(responseSd)) || any(responseSd <= 0)) {
        objectiveFailure <<- paste0("invalid response standard deviation; residual=",
          paste(signif(candidate$residual, 4), collapse = "/"),
          ", f.range=", paste(signif(range(layout$response$f), 4), collapse = "/"),
          ", lengths(f/type)=", length(layout$response$f), "/",
          length(layout$response$etype), ", types=",
          paste(unique(layout$response$etype), collapse = "/"),
          ", sd.range=", if (inherits(responseSd, "try-error")) "error" else
            paste(signif(range(responseSd), 4), collapse = "/"))
        return(-Inf)
      }
      responseLog <- -.5 * ((layout$response$y - layout$response$f) /
        responseSd)^2 - log(responseSd) - .5 * log(2 * pi)
      if (any(!is.finite(responseLog))) {
        objectiveFailure <<- "non-finite response objective"; return(-Inf)
      }
      value <- value + sum(responseLog) / nrow(E)
    }
    value
  }
  validCandidate <- function(internal) {
    candidate <- try(copulaGaussianFremScoreMaterialize(internal, layout,
      buildVine = FALSE), silent = TRUE)
    if (inherits(candidate, "try-error")) return(FALSE)
    completeState <- materializeState(candidate)
    residual <- completeState$residual
    evaluated <- try(if (completeContinuous && !hasReference)
      copulaGaussianFremEvaluateMargins(residual, candidate$margins) else
      if (completeAugmented || hasReference)
        copulaGaussianFremAugmentedEvaluateMargins(residual,
          candidate$margins, categoricalUniform, referenceUniform) else
        stop("incomplete categorical state"), silent = TRUE)
    if (inherits(evaluated, "try-error") || !all(evaluated$valid)) return(FALSE)
    if (!is.null(layout$response)) {
      if (is.function(layout$response$evaluate) && hasReference) {
        responseLog <- try(layout$response$evaluate(
          completeState$absolute[, seq_len(layout$dEta), drop = FALSE],
          candidate$residual), silent = TRUE)
        if (inherits(responseLog, "try-error") ||
            any(!is.finite(responseLog))) return(FALSE)
        return(TRUE)
      }
      responseSd <- try(error(layout$response$f, candidate$residual,
        layout$response$etype), silent = TRUE)
      if (inherits(responseSd, "try-error") ||
          any(!is.finite(responseSd)) || any(responseSd <= 0)) return(FALSE)
    }
    TRUE
  }
  h <- finiteDifference * sqrt(gain)
  analyticEligible <- isTRUE(analyticScore) && completeContinuous &&
    !hasReference
  analytic <- if (analyticEligible)
    copulaGaussianFremCompleteScoreInternal(current, layout, E, w, h) else NULL
  numericalFallbackCount <- 0L
  if (!is.null(analytic)) {
    if (isTRUE(analytic$oneSided) && !isTRUE(adaptMetric))
      stop("score-sa requires centered numerical components after finite initialization")
    gradient <- analytic$gradient
  } else {
    gradient <- numeric(length(current))
    for (j in seq_along(current)) {
      plus <- minus <- current
      plus[j] <- plus[j] + h; minus[j] <- minus[j] - h
      fp <- objective(plus); fm <- objective(minus)
      if (is.finite(fp) && is.finite(fm)) {
        gradient[j] <- (fp - fm) / (2 * h)
      } else {
        numericalFallbackCount <- numericalFallbackCount + 1L
        if (!isTRUE(adaptMetric))
          stop("score-sa requires centered global differences after finite initialization; coordinate ",
            j, " had finite(+/0/-)=", is.finite(fp), "/",
            is.finite(objective(current)), "/", is.finite(fm),
            "; last failure: ", objectiveFailure %||% "unknown")
        f0 <- objective(current)
        if (is.finite(fp) && is.finite(f0)) gradient[j] <- (fp - f0) / h else
          if (is.finite(fm) && is.finite(f0)) gradient[j] <- (f0 - fm) / h else
            gradient[j] <- 0
      }
    }
  }
  if (is.null(state)) state <- list(iteration = 0L, projectionCount = 0L,
    average = current, averageCount = 0L,
    runningSquare = pmax(gradient^2, 1e-2), projectionCentre = current)
  for (field in c("backtrackCount", "postFreezeBacktrackCount",
      "noMoveCount", "postFreezeNoMoveCount", "postFreezeProjectionCount",
      "numericalFallbackCount"))
    if (is.null(state[[field]])) state[[field]] <- 0L
  state$numericalFallbackCount <- state$numericalFallbackCount +
    numericalFallbackCount
  if (is.null(state$runningSquare) ||
      length(state$runningSquare) != length(gradient))
    state$runningSquare <- pmax(gradient^2, 1e-2)
  metricGain <- if (isTRUE(adaptMetric)) (state$iteration + 10)^-.7 else 0
  if (metricGain > 0)
    state$runningSquare <- state$runningSquare + metricGain *
      (gradient^2 - state$runningSquare)
  preconditioner <- pmin(10, pmax(.1,
    1 / sqrt(pmax(state$runningSquare, 1e-8))))
  direction <- preconditioner * gradient
  proposal <- current + gain * scoreScale * direction
  if (is.null(state$projectionCentre) ||
      length(state$projectionCentre) != length(current))
    state$projectionCentre <- current
  projectionLower <- state$projectionCentre - projection
  projectionUpper <- state$projectionCentre + projection
  projected <- pmin(projectionUpper, pmax(projectionLower, proposal))
  projectionEvent <- any(abs(projected - proposal) > 0)
  state$projectionCount <- state$projectionCount + as.integer(projectionEvent)
  if (projectionEvent && !isTRUE(adaptMetric))
    state$postFreezeProjectionCount <- state$postFreezeProjectionCount + 1L
  ## Projection prevents parameter explosion; an invalid finite-precision PIT
  ## is handled by deterministic step halving, not clipping the density.
  localScale <- 1
  projectedValid <- if (is.null(analytic)) is.finite(objective(projected)) else
    validCandidate(projected)
  while (!isTRUE(projectedValid) && localScale > 2^-12) {
    localScale <- localScale / 2
    projected <- pmin(projectionUpper, pmax(projectionLower,
      current + localScale * gain * scoreScale * direction))
    projectedValid <- if (is.null(analytic)) is.finite(objective(projected)) else
      validCandidate(projected)
  }
  if (localScale < 1) {
    state$backtrackCount <- state$backtrackCount + 1L
    if (!isTRUE(adaptMetric))
      state$postFreezeBacktrackCount <- state$postFreezeBacktrackCount + 1L
  }
  if (!isTRUE(projectedValid)) {
    projected <- current
    state$noMoveCount <- state$noMoveCount + 1L
    if (!isTRUE(adaptMetric))
      state$postFreezeNoMoveCount <- state$postFreezeNoMoveCount + 1L
  }
  state$iteration <- state$iteration + 1L
  if (isTRUE(average)) {
    if (!isTRUE(state$averagingStarted)) {
      state$average <- projected; state$averageCount <- 1L
      state$averagingStarted <- TRUE
    } else {
      state$averageCount <- state$averageCount + 1L
      state$average <- state$average +
        (projected - state$average) / state$averageCount
    }
    if (is.null(state$scoreAverage)) {
      state$scoreAverage <- gradient; state$scoreAverageCount <- 1L
    } else {
      state$scoreAverageCount <- state$scoreAverageCount + 1L
      state$scoreAverage <- state$scoreAverage +
        (gradient - state$scoreAverage) / state$scoreAverageCount
    }
  }
  state$internal <- projected
  reportedInternal <- if (isTRUE(useAverage) &&
      isTRUE(state$averagingStarted)) state$average else projected
  reported <- copulaGaussianFremScoreMaterialize(reportedInternal, layout)
  ## A complete-data objective value is a diagnostic, not part of score-SA.
  ## Evaluate it only for a numerical-score step or the terminal averaged fit.
  reportedValue <- if (is.null(analytic) || isTRUE(useAverage))
    objective(reportedInternal) else NA_real_
  list(beta = reported$beta, margins = reported$margins,
    sd = copulaMarginScales(reported$margins), vine = reported$vine,
    delta = reported$delta, residual = reported$residual,
    value = reportedValue,
    score = gradient, scoreMax = max(abs(gradient)),
    scoreAverage = state$scoreAverage %||% rep(NA_real_, length(gradient)),
    scoreAverageMax = if (is.null(state$scoreAverage)) NA_real_ else
      max(abs(state$scoreAverage)),
    finiteDifference = h, gain = gain, scoreScale = scoreScale,
    scoreMethod = if (is.null(analytic)) "global-centered-difference" else
      "analytic-with-declared-local-numerical-components",
    numericalScoreComponents = if (is.null(analytic)) "all" else
      analytic$numerical,
    preconditionerMin = min(preconditioner),
    preconditionerMax = max(preconditioner), metricGain = metricGain,
    averaging = isTRUE(average), usedAverage = isTRUE(useAverage) &&
      isTRUE(state$averagingStarted),
    state = state, backend = "gaussian-copula-frem-score-sa",
    scoreTheory = list(
      route = "Delyon-Section-8.2-equation-74",
      fixedFingerprint = TRUE,
      coordinateSystem = "bounded internal coordinates; native = tau(internal)",
      augmentedOrPartiallyMarginalizedScore = TRUE,
      categoricalAugmentation = if (length(discrete))
        "fixed-support latent uniforms; exact rectangle recovered on integration" else
        "none",
      movingSupportAugmentation = if (hasReference)
        "fixed percentile coordinates; quantile and response paths differentiated" else
        "none",
      metricFrozen = !isTRUE(adaptMetric),
      metricBounds = c(0.1 * scoreScale, 10 * scoreScale),
      finiteDifferenceOrder = 2L,
      analyticExcept = if (is.null(analytic)) "all coordinates" else
        analytic$numerical,
      finiteDifferenceSchedule = "h_k = h0 * sqrt(gamma_k)",
      weightedTruncationBias =
        "O(sum gamma_k * h_k^2) = O(sum gamma_k^2)",
      projectionCount = state$projectionCount,
      backtrackCount = state$backtrackCount,
      noMoveCount = state$noMoveCount,
      postFreezeProjectionCount = state$postFreezeProjectionCount,
      postFreezeBacktrackCount = state$postFreezeBacktrackCount,
      postFreezeNoMoveCount = state$postFreezeNoMoveCount,
      numericalFallbackCount = state$numericalFallbackCount,
      runtimeConditionsObserved = !isTRUE(adaptMetric) &&
        state$postFreezeProjectionCount == 0L &&
        state$postFreezeBacktrackCount == 0L &&
        state$postFreezeNoMoveCount == 0L,
      caveat = paste0("Conditional on score smoothness/moments, controlled-",
        "MCMC Poisson-equation bounds, stable iterates, and negligible ",
        "floating-point error.")))
}

## Record the current controlled-MCMC batch for the score estimator.  Unlike
## the legacy common-Q route, this stores no retained particle history.
