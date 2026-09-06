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
    ## The error-type partition of the observations is a property of the data,
    ## not of the candidate, but the residual gradient was rebuilding it --
    ## rep_len, sort(unique(.)) and one comparison per type over every
    ## observation -- on every gradient evaluation, of which there are several
    ## per iteration. Store the partition with the layout instead.
    responseType <- rep_len(as.integer(response$etype), length(response$f))
    response$typeRows <- lapply(sort(unique(responseType)),
      function(type) which(responseType == type))
    names(response$typeRows) <- as.character(sort(unique(responseType)))
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
    internal, layout, E, w, numericalStep, skipMargins = integer(),
    referenceUniform = NULL, perRow = FALSE) {
  candidate <- try(copulaGaussianFremScoreMaterialize(internal, layout,
    buildVine = FALSE),
    silent = TRUE)
  if (inherits(candidate, "try-error")) return(NULL)
  residual <- if (layout$hasDesign)
    E - copulaLocation(layout$X, candidate$beta, layout$locMap) else
    sweep(E, 2L, candidate$delta, "-")
  if (length(skipMargins)) {
    referenceUniform <- as.matrix(referenceUniform)
    z <- logMargin <- matrix(NA_real_, nrow(E), layout$d)
    valid <- rep(TRUE, nrow(E))
    for (j in seq_len(layout$d)) if (j %in% skipMargins) {
      u <- referenceUniform[, j]
      ok <- is.finite(u) & u > 0 & u < 1
      valid <- valid & ok; z[ok, j] <- stats::qnorm(u[ok]); logMargin[ok, j] <- 0
    } else {
      local <- try(copulaGaussianFremEvaluateMargins(
        matrix(residual[, j], ncol = 1L), list(candidate$margins[[j]])),
        silent = TRUE)
      if (inherits(local, "try-error")) return(NULL)
      z[, j] <- local$z[, 1L]; logMargin[, j] <- local$logMargin[, 1L]
      valid <- valid & local$valid
    }
    evaluated <- list(z = z, logMargin = logMargin, valid = valid)
  } else evaluated <- try(copulaGaussianFremEvaluateMargins(
    residual, candidate$margins), silent = TRUE)
  if (inherits(evaluated, "try-error") || !all(evaluated$valid)) return(NULL)
  R <- candidate$correlation
  Omega <- solve(R); zInfluence <- evaluated$z %*% (diag(layout$d) - Omega)
  parameterGradient <- numeric(layout$nMargin)
  marginRows <- vector("list", layout$d)
  scoreX <- matrix(NA_real_, nrow(E), layout$d)
  numerical <- character(); oneSided <- FALSE
  for (j in seq_len(layout$d)) {
    if (j %in% skipMargins) {
      ## A fixed-reference coordinate has uniform density in its retained
      ## percentile and a parameter-invariant Gaussian score. Its complete
      ## population score is therefore zero here; any response-path and direct-
      ## row contribution is supplied by the declared hybrid coordinates.
      index <- layout$marginLayout$index[[j]]
      if (length(index)) parameterGradient[index] <- 0
      scoreX[, j] <- 0
      next
    }
    local <- copulaGaussianMarginNativeScore(residual[, j],
      candidate$margins[[j]], evaluated$z[, j], zInfluence[, j],
      numericalStep)
    if (is.null(local)) return(NULL)
    index <- layout$marginLayout$index[[j]]
    if (length(index)) {
      contribution <- local$parameter[, candidate$margins[[j]]$free,
        drop = FALSE]
      parameterGradient[index] <- colSums(w * contribution)
      if (isTRUE(perRow)) marginRows[[j]] <- contribution
    }
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
  responseSd <- common <- NULL
  if (layout$nResidual) {
    responseSd <- error(layout$response$f, candidate$residual,
      layout$response$etype)
    difference <- layout$response$y - layout$response$f
    common <- difference^2 / responseSd^3 - 1 / responseSd
    full <- numeric(length(candidate$residual))
    for (label in names(layout$response$typeRows)) {
      type <- as.integer(label)
      rows <- layout$response$typeRows[[label]]
      ia <- 2L * type - 1L; ib <- 2L * type
      a <- candidate$residual[ia]; b <- candidate$residual[ib]
      full[ia] <- sum(common[rows] * a / responseSd[rows]) / nrow(E)
      full[ib] <- sum(common[rows] * b * layout$response$f[rows]^2 /
        responseSd[rows]) / nrow(E)
    }
    residualGradient <- full[layout$residualFree]
  }
  nativeGradient <- c(locationGradient, parameterGradient,
    correlationGradient, residualGradient)
  ## The same score, kept per latent row instead of summed over them. The
  ## Fisher information is the covariance of the per-subject score, so it
  ## cannot be formed from the aggregate; and every block above is already a
  ## weighted sum over rows, so keeping the terms costs nothing but the
  ## allocation. The correlation block is the only one that is not written that
  ## way -- it goes through the scatter matrix S -- but its objective is linear
  ## in S, so the row contribution is the same linear map applied to that row's
  ## own outer product.
  rows <- if (!isTRUE(perRow)) NULL else {
    answer <- matrix(0, nrow(E), length(nativeGradient))
    at <- 0L
    if (layout$nLocation) {
      if (layout$hasDesign) for (k in seq_along(layout$locationIndex)) {
        l <- layout$locationIndex[k]
        derivativeMap <- sweep(matrix(layout$locMap[l, ], nrow(E), layout$d,
          byrow = TRUE), 1L, layout$X[, l], "*")
        answer[, at + k] <- -w * rowSums(scoreX * derivativeMap)
      } else for (k in seq_len(layout$nLocation))
        answer[, at + k] <- -w * scoreX[, k]
      at <- at + layout$nLocation
    }
    if (layout$nMargin) {
      for (j in seq_len(layout$d)) {
        columns <- layout$marginLayout$index[[j]]
        if (!length(columns)) next
        if (j %in% skipMargins) {
          answer[, at + columns] <- 0
        } else {
          local <- marginRows[[j]]
          answer[, at + columns] <- w * local
        }
      }
      at <- at + layout$nMargin
    }
    if (layout$nEdge) {
      ## Extract the linear map from S to each angle gradient once, then apply
      ## it to each row's own contribution to S.
      ## The angle gradient is affine in S, not linear: the term in Omega
      ## alone survives at S = 0. That constant belongs to the rows in the
      ## proportions they carry, which are their weights, since the weights
      ## sum to one.
      constant <- copulaGaussianCorrelationAngleObjective(angles,
        matrix(0, layout$d, layout$d), layout$d)$gradient
      answer[, at + seq_len(layout$nEdge)] <-
        answer[, at + seq_len(layout$nEdge)] - outer(w, constant)
      base <- copulaGaussianCorrelationAngleObjective(angles, S,
        layout$d)$gradient
      for (a in seq_len(layout$d)) for (b in seq_len(a)) {
        bump <- matrix(0, layout$d, layout$d)
        bump[a, b] <- bump[b, a] <- 1
        slope <- copulaGaussianCorrelationAngleObjective(angles, S + bump,
          layout$d)$gradient - base
        ## The bump raises both symmetric entries of S at once, and a row
        ## contributes its outer product to both, so the off-diagonal term is
        ## not doubled again here.
        contribution <- if (a == b) evaluated$z[, a]^2 else
          evaluated$z[, a] * evaluated$z[, b]
        answer[, at + seq_len(layout$nEdge)] <-
          answer[, at + seq_len(layout$nEdge)] -
          outer(w * contribution, slope)
      }
      at <- at + layout$nEdge
    }
    if (layout$nResidual) {
      per <- matrix(0, nrow(E), length(candidate$residual))
      for (label in names(layout$response$typeRows)) {
        type <- as.integer(label)
        these <- layout$response$typeRows[[label]]
        ia <- 2L * type - 1L; ib <- 2L * type
        a <- candidate$residual[ia]; b <- candidate$residual[ib]
        per[, ia] <- per[, ia] + as.numeric(rowsum(
          common[these] * a / responseSd[these],
          layout$response$row[these], reorder = FALSE))[
            match(seq_len(nrow(E)), sort(unique(layout$response$row[these])))]
        per[, ib] <- per[, ib] + as.numeric(rowsum(
          common[these] * b * layout$response$f[these]^2 / responseSd[these],
          layout$response$row[these], reorder = FALSE))[
            match(seq_len(nrow(E)), sort(unique(layout$response$row[these])))]
      }
      per[is.na(per)] <- 0
      answer[, at + seq_along(layout$residualFree)] <-
        per[, layout$residualFree, drop = FALSE] / nrow(E)
    }
    answer
  }
  derivative <- copulaScoreInternalDerivative(internal,
    layout$lower, layout$upper)
  gradient <- nativeGradient * derivative
  if (length(gradient) != length(internal) || any(!is.finite(gradient)))
    return(NULL)
  list(gradient = gradient, nativeGradient = nativeGradient,
    nativeRows = rows,
    numerical = unique(numerical), candidate = candidate,
    residual = residual, evaluated = evaluated, oneSided = oneSided)
}

## Hybrid score for fixed-reference continuous coordinates.  At the current
## parameter value E and Q_theta(U) coincide, so the ordinary complete score is
## valid for every coordinate that does not move Q_theta(U).  Only population
## locations feeding a referenced eta and parameters of referenced margins need
## the full pathwise response derivative.  Evaluate centered differences for
## those coordinates and retain the shared analytic margin, correlation and
## residual scores everywhere else.
copulaGaussianFremReferenceScoreInternal <- function(
    internal, layout, E, w, numericalStep, referenceUniform, objective,
    materializeState, perRow = FALSE, objectiveRows = NULL) {
  referenceUniform <- as.matrix(referenceUniform)
  referenced <- which(colSums(is.finite(referenceUniform)) > 0L)
  if (!length(referenced))
    return(copulaGaussianFremCompleteScoreInternal(
      internal, layout, E, w, numericalStep))
  base <- copulaGaussianFremCompleteScoreInternal(
    internal, layout, E, w, numericalStep, skipMargins = referenced,
    referenceUniform = referenceUniform, perRow = perRow)
  if (is.null(base)) return(NULL)
  ## The corrections below are applied in internal coordinates, so the rows
  ## are carried that way too.
  internalRows <- if (!isTRUE(perRow) || is.null(base$nativeRows)) NULL else
    sweep(base$nativeRows, 2L,
      copulaScoreInternalDerivative(internal, layout$lower, layout$upper), "*")

  numericalIndex <- integer(); responsePathIndex <- integer()
  movingEta <- intersect(referenced, seq_len(layout$dEta))
  if (length(movingEta) && is.function(layout$response$evaluate) &&
      layout$nLocation) {
    if (layout$hasDesign) {
      local <- which(vapply(layout$locationIndex, function(index)
        any(layout$locMap[index, movingEta, drop = FALSE] != 0), logical(1)))
    } else local <- intersect(movingEta, seq_len(layout$nLocation))
    responsePathIndex <- c(responsePathIndex, local)
  }
  for (j in referenced) {
    index <- layout$marginLayout$index[[j]]
    if (length(index)) {
      internalIndex <- layout$nLocation + index
      if (j %in% movingEta && all(is.finite(referenceUniform[, j])))
        responsePathIndex <- c(responsePathIndex, internalIndex) else
        numericalIndex <- c(numericalIndex, internalIndex)
    }
  }
  numericalIndex <- sort(unique(numericalIndex))
  responsePathIndex <- sort(unique(responsePathIndex))
  if (!length(c(numericalIndex, responsePathIndex))) {
    base$numerical <- unique(c(base$numerical, "fixed-reference.none"))
    return(base)
  }

  gradient <- base$gradient
  oneSided <- isTRUE(base$oneSided)
  ## A response-gradient callback differentiates the structural model once per
  ## referenced eta coordinate. Parameter-specific directions then use only
  ## quantile/materialization calls and a chain-rule dot product.
  if (length(responsePathIndex) && is.function(layout$response$gradient)) {
    candidate <- base$candidate
    complete <- materializeState(candidate)
    responseGradient <- try(layout$response$gradient(
      complete$absolute[, seq_len(layout$dEta), drop = FALSE],
      candidate$residual, movingEta, numericalStep), silent = TRUE)
    if (inherits(responseGradient, "try-error") ||
        any(dim(responseGradient) != c(nrow(E), layout$dEta)) ||
        any(!is.finite(responseGradient))) return(NULL)
    for (index in responsePathIndex) {
      plus <- minus <- internal
      plus[index] <- plus[index] + numericalStep
      minus[index] <- minus[index] - numericalStep
      plusCandidate <- try(copulaGaussianFremScoreMaterialize(
        plus, layout, buildVine = FALSE), silent = TRUE)
      minusCandidate <- try(copulaGaussianFremScoreMaterialize(
        minus, layout, buildVine = FALSE), silent = TRUE)
      if (inherits(plusCandidate, "try-error") ||
          inherits(minusCandidate, "try-error")) return(NULL)
      plusState <- materializeState(plusCandidate)
      minusState <- materializeState(minusCandidate)
      direction <- (plusState$absolute[, seq_len(layout$dEta), drop = FALSE] -
        minusState$absolute[, seq_len(layout$dEta), drop = FALSE]) /
        (2 * numericalStep)
      contribution <- rowSums(responseGradient * direction) / nrow(E)
      gradient[index] <- gradient[index] + sum(contribution)
      if (!is.null(internalRows))
        internalRows[, index] <- internalRows[, index] + contribution
    }
  } else numericalIndex <- sort(unique(c(numericalIndex,
    responsePathIndex)))
  for (index in numericalIndex) {
    plus <- minus <- internal
    plus[index] <- plus[index] + numericalStep
    minus[index] <- minus[index] - numericalStep
    fp <- objective(plus); fm <- objective(minus)
    if (!is.finite(fp) || !is.finite(fm)) return(NULL)
    gradient[index] <- (fp - fm) / (2 * numericalStep)
    ## The same two evaluations, not summed: a difference of per-row objectives
    ## is the per-row difference, so the rows come free with the aggregate.
    if (!is.null(internalRows) && is.function(objectiveRows)) {
      rowsPlus <- objectiveRows(plus); rowsMinus <- objectiveRows(minus)
      if (is.null(rowsPlus) || is.null(rowsMinus)) return(NULL)
      internalRows[, index] <- (rowsPlus - rowsMinus) / (2 * numericalStep)
    }
  }
  base$gradient <- gradient
  base$internalRows <- internalRows
  base$nativeGradient <- gradient / copulaScoreInternalDerivative(
    internal, layout$lower, layout$upper)
  base$numerical <- unique(c(base$numerical,
    if (length(responsePathIndex) && is.function(layout$response$gradient))
      "fixed-reference.chain-rule-path" else
        "fixed-reference.path-coordinates"))
  base$oneSided <- oneSided
  base
}

## Complete-data score for the parameter scale with the natural parameters
## held fixed. Every coordinate -- population locations through their typical
## values, native margin parameters, correlation angles and residual-error
## parameters -- is a centered difference of the per-row complete objective,
## so the per-subject rows that the Fisher information needs come from the
## same two evaluations as the aggregate. The differencing step is the
## gain-scaled h of the paper's numerical-error argument; the response term is
## evaluated at fixed predictions and so is constant in every coordinate but
## the residual ones.
copulaGaussianFremNaturalScoreInternal <- function(
    internal, layout, numericalStep, objectiveRows, perRow = FALSE) {
  candidate <- try(copulaGaussianFremScoreMaterialize(internal, layout,
    buildVine = FALSE), silent = TRUE)
  if (inherits(candidate, "try-error")) return(NULL)
  gradient <- numeric(length(internal)); rows <- NULL
  for (index in seq_along(internal)) {
    plus <- minus <- internal
    plus[index] <- plus[index] + numericalStep
    minus[index] <- minus[index] - numericalStep
    rowsPlus <- objectiveRows(plus); rowsMinus <- objectiveRows(minus)
    if (is.null(rowsPlus) || is.null(rowsMinus)) return(NULL)
    difference <- (rowsPlus - rowsMinus) / (2 * numericalStep)
    if (any(!is.finite(difference))) return(NULL)
    gradient[index] <- sum(difference)
    if (isTRUE(perRow)) {
      if (is.null(rows)) rows <- matrix(0, length(difference), length(internal))
      rows[, index] <- difference
    }
  }
  list(gradient = gradient,
    nativeGradient = gradient / copulaScoreInternalDerivative(internal,
      layout$lower, layout$upper),
    internalRows = rows, nativeRows = NULL,
    numerical = "natural-psi.centered-differences",
    candidate = candidate, oneSided = FALSE)
}

copulaGaussianFremPopulationScoreStep <- function(
    E, w, margins0, vine0, d, dEta, gain,
    X = NULL, locMap = NULL, beta0 = NULL, betaFree = NULL,
    state = NULL, scoreScale = 1, finiteDifference = 1e-4,
    projection = 24, withMu = TRUE, adaptMetric = TRUE,
    average = FALSE, useAverage = FALSE, response = NULL,
    analyticScore = TRUE, categoricalUniform = NULL,
    referenceUniform = NULL,
    populationScale = c("transformed-additive", "parameter"),
    transform = NULL, subject = NULL, preheating = FALSE,
    smoothing = FALSE, freeze = FALSE, deltaGain = gain) {
  populationScale <- match.arg(populationScale)
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
  layout$populationScale <- populationScale
  layout$transform <- if (is.null(transform)) rep(0L, dEta) else
    rep_len(as.integer(transform), dEta)
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
  ## Parameter-scale fit with the natural parameters held fixed: no reference
  ## percentile has been supplied for the eta coordinates, so the score is the
  ## ordinary complete-data score in psi (see copulaScoreBatchUpdate).
  naturalPsi <- identical(populationScale, "parameter") &&
    !any(is.finite(referenceUniform[, seq_len(dEta), drop = FALSE]))
  if (naturalPsi && hasReference)
    stop("parameter-scale psi augmentation cannot be combined with fixed-reference coordinates")
  if (naturalPsi && (!layout$hasDesign || length(discrete)))
    stop("parameter-scale psi augmentation requires a design-aware continuous batch")
  ## The fixed reference is constant for the whole step, so its Gaussian score
  ## is too. Forming it once lets the parameter-scale path go through
  ## inverse_score, which for families that have a closed form skips the
  ## uniform entirely, and removes a qnorm over every row from every objective
  ## and gradient evaluation.
  referenceScore <- if (hasReference) stats::qnorm(referenceUniform) else
    referenceUniform
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
    if (identical(populationScale, "parameter")) {
      typical <- copulaWorkingToNatural(
        location[, seq_len(dEta), drop = FALSE], layout$transform)
      if (naturalPsi) {
        ## psi is the latent variable and does not move with the candidate;
        ## only its typical value does.
        absolute <- E[, seq_len(dEta), drop = FALSE]
        psi <- copulaWorkingToNatural(absolute, layout$transform)
        residual <- E
        residual[, seq_len(dEta)] <- absolute -
          location[, seq_len(dEta), drop = FALSE]
        return(list(residual = residual, absolute = absolute,
          typical = typical, natural = psi,
          conditioning = if (d > dEta) E[, dEta + seq_len(d - dEta),
            drop = FALSE] else matrix(numeric(), nrow(E), 0L)))
      }
      if (any(!is.finite(referenceUniform[, seq_len(dEta), drop = FALSE])))
        stop("parameter-scale score materialization requires parameter reference uniforms")
      psi <- copulaNaturalMarginsFromScore(referenceScore[, seq_len(dEta),
        drop = FALSE], typical, candidate$margins[seq_len(dEta)])
      absolute <- copulaNaturalToWorking(psi, layout$transform)
      residual <- E
      residual[, seq_len(dEta)] <- absolute -
        location[, seq_len(dEta), drop = FALSE]
      return(list(residual = residual, absolute = absolute,
        typical = typical, natural = psi,
        conditioning = if (d > dEta) E[, dEta + seq_len(d - dEta),
          drop = FALSE] else matrix(numeric(), nrow(E), 0L)))
    }
    residual <- E - location
    for (j in seq_len(d)) {
      rows <- is.finite(referenceUniform[, j])
      if (any(rows)) residual[rows, j] <- candidate$margins[[j]]$quantile(
        referenceUniform[rows, j], candidate$margins[[j]]$parameters)
    }
    list(residual = residual, absolute = location + residual)
  }
  ## Complete-data population terms for the psi augmentation: natural margins
  ## evaluated at the fixed psi against the candidate's typical values and
  ## parameters, observed covariates against their margins, and the working-
  ## coordinate Jacobian, which does not depend on the candidate but keeps the
  ## objective equal to the complete log density rather than to it up to a
  ## constant.
  naturalJacobian <- if (naturalPsi) rowSums(copulaWorkingLogJacobian(
    E[, seq_len(dEta), drop = FALSE], layout$transform)) else NULL
  evaluateNatural <- function(candidate, completeState) {
    parameter <- copulaNaturalMarginsEvaluate(completeState$natural,
      completeState$typical, candidate$margins[seq_len(dEta)])
    if (d == dEta) return(parameter)
    covariate <- copulaGaussianFremEvaluateMargins(
      E[, dEta + seq_len(d - dEta), drop = FALSE],
      candidate$margins[dEta + seq_len(d - dEta)])
    list(z = cbind(parameter$z, covariate$z),
      logMargin = cbind(parameter$logMargin, covariate$logMargin),
      valid = parameter$valid & covariate$valid)
  }
  ## Response log density at the fixed predictions, per observation. Only the
  ## residual-error parameters move it; it is what the non-reference objective
  ## already adds, kept per row here so the psi-route score can be formed per
  ## subject.
  fixedResponseLog <- function(candidate) {
    responseSd <- try(error(layout$response$f, candidate$residual,
      layout$response$etype), silent = TRUE)
    if (inherits(responseSd, "try-error") ||
        any(!is.finite(responseSd)) || any(responseSd <= 0)) return(NULL)
    -.5 * ((layout$response$y - layout$response$f) / responseSd)^2 -
      log(responseSd) - .5 * log(2 * pi)
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
      if (naturalPsi) {
        evaluated <- evaluateNatural(candidate, completeState)
        if (!all(evaluated$valid)) stop("invalid complete-data margin")
        sum(w * (copulaGaussianCopulaLogDensity(evaluated$z,
          candidate$correlation) + rowSums(evaluated$logMargin) +
          naturalJacobian))
      } else if (completeAugmented || hasReference) {
        evaluated <- if (completeContinuous && !hasReference)
          copulaGaussianFremEvaluateMargins(residual, candidate$margins) else
          copulaGaussianFremAugmentedEvaluateMargins(residual,
            candidate$margins, categoricalUniform, referenceUniform,
            referenceScore)
        if (!all(evaluated$valid)) stop("invalid complete-data margin")
        sum(w * (copulaGaussianCopulaLogDensity(evaluated$z,
          candidate$correlation) + rowSums(evaluated$logMargin)))
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
  ## The same objective, kept per latent row. It is a weighted sum over rows
  ## plus a response term that is a sum over observations, and each observation
  ## belongs to one row, so no extra evaluation is needed to know what each row
  ## contributed -- only the restraint not to add them up.
  objectiveRows <- function(internal) {
    candidate <- try(copulaGaussianFremScoreMaterialize(internal, layout,
      buildVine = !completeAugmented), silent = TRUE)
    if (inherits(candidate, "try-error")) return(NULL)
    completeState <- materializeState(candidate)
    residual <- completeState$residual
    value <- try({
      evaluated <- if (naturalPsi)
        evaluateNatural(candidate, completeState) else
        if (completeContinuous && !hasReference)
        copulaGaussianFremEvaluateMargins(residual, candidate$margins) else
        copulaGaussianFremAugmentedEvaluateMargins(residual,
          candidate$margins, categoricalUniform, referenceUniform,
          referenceScore)
      if (!all(evaluated$valid)) stop("invalid complete-data margin")
      w * (copulaGaussianCopulaLogDensity(evaluated$z,
        candidate$correlation) + rowSums(evaluated$logMargin) +
        (if (naturalPsi) naturalJacobian else 0))
    }, silent = TRUE)
    if (inherits(value, "try-error") || any(!is.finite(value))) return(NULL)
    if (!is.null(layout$response) && is.function(layout$response$evaluate) &&
        hasReference) {
      responseLog <- try(layout$response$evaluate(
        completeState$absolute[, seq_len(layout$dEta), drop = FALSE],
        candidate$residual), silent = TRUE)
      if (inherits(responseLog, "try-error") || any(!is.finite(responseLog)))
        return(NULL)
      value <- value + as.numeric(rowsum(responseLog, layout$response$row,
        reorder = TRUE)) / nrow(E)
    } else if (!is.null(layout$response) && naturalPsi) {
      responseLog <- fixedResponseLog(candidate)
      if (is.null(responseLog) || any(!is.finite(responseLog))) return(NULL)
      value <- value + as.numeric(rowsum(responseLog, layout$response$row,
        reorder = TRUE)) / nrow(E)
    }
    value
  }
  validCandidate <- function(internal) {
    candidate <- try(copulaGaussianFremScoreMaterialize(internal, layout,
      buildVine = FALSE), silent = TRUE)
    if (inherits(candidate, "try-error")) return(FALSE)
    completeState <- materializeState(candidate)
    residual <- completeState$residual
    evaluated <- try(if (naturalPsi)
      evaluateNatural(candidate, completeState) else
      if (completeContinuous && !hasReference)
      copulaGaussianFremEvaluateMargins(residual, candidate$margins) else
      if (completeAugmented || hasReference)
        copulaGaussianFremAugmentedEvaluateMargins(residual,
          candidate$margins, categoricalUniform, referenceUniform,
          referenceScore) else
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
  metricTrace <- isTRUE(getOption("saemix.fisherTrace", FALSE))
  ## The information is built on every fit, because it is what the standard
  ## errors are made of and it costs only the per-subject score rows the score
  ## already computes: measured overhead under one per cent, with the estimates
  ## bit-identical either way. It does not steer the recursion.
  fisherInformation <- metricTrace ||
    isTRUE(getOption("saemix.fisherInformation", TRUE))
  analyticEligible <- isTRUE(analyticScore) && completeContinuous
  ## Per-subject scores are only formed on the complete analytic path, where
  ## every block decomposes exactly. Their weighted sum has to reproduce the
  ## aggregate they were split out of, and that is checked rather than assumed
  ## whenever they are asked for.
  wantRows <- isTRUE(getOption("saemix.scorePerRow", FALSE)) ||
    (fisherInformation && !is.null(subject))
  analytic <- if (!analyticEligible) NULL else if (naturalPsi)
    copulaGaussianFremNaturalScoreInternal(current, layout, h,
      objectiveRows, perRow = wantRows) else if (hasReference)
    copulaGaussianFremReferenceScoreInternal(current, layout, E, w, h,
      referenceUniform, objective, materializeState,
      perRow = wantRows, objectiveRows = objectiveRows) else
    copulaGaussianFremCompleteScoreInternal(current, layout, E, w, h,
      perRow = wantRows)
  if (wantRows && !is.null(analytic) && !is.null(analytic$internalRows)) {
    reconstructed <- colSums(analytic$internalRows)
    scale <- pmax(1, abs(analytic$gradient))
    if (max(abs(reconstructed - analytic$gradient) / scale) > 1e-8)
      stop("per-row score does not sum to the aggregate score (reference path)")
  }
  if (wantRows && !is.null(analytic) && !is.null(analytic$nativeRows) &&
      is.null(analytic$internalRows)) {
    reconstructed <- colSums(analytic$nativeRows)
    scale <- pmax(1, abs(analytic$nativeGradient))
    if (max(abs(reconstructed - analytic$nativeGradient) / scale) > 1e-8)
      stop("per-row score does not sum to the aggregate score
",
        "  blocks: location ", layout$nLocation, ", margin ", layout$nMargin,
        ", edge ", layout$nEdge, ", residual ", layout$nResidual, "
",
        "  aggregate:    ",
        paste(sprintf("%.6g", analytic$nativeGradient), collapse = " "), "
",
        "  from rows:    ",
        paste(sprintf("%.6g", reconstructed), collapse = " "))
  }
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
    projectionCentre = current)
  for (field in c("backtrackCount", "postFreezeBacktrackCount",
      "noMoveCount", "postFreezeNoMoveCount", "postFreezeProjectionCount",
      "numericalFallbackCount"))
    if (is.null(state[[field]])) state[[field]] <- 0L
  state$numericalFallbackCount <- state$numericalFallbackCount +
    numericalFallbackCount
  fisherReady <- FALSE
  if (is.null(state$metricUpdateCount)) state$metricUpdateCount <- 0L

  ## The score covariance, accumulated across iterations rather than estimated
  ## from one batch. Delattre & Kuhn's information estimate needs the posterior
  ## mean of each subject's complete-data score, and a single batch estimates
  ## that badly: the posterior draws mix slowly, with an effective sample
  ## fraction of a few percent that does not improve when the chain is made
  ## longer. Averaging along the recursion instead, as Fisher-SGD does, spends
  ## draws the fit is producing anyway and sees a new MCMC state every
  ## iteration.
  ##
  ##   Delta_i^k = (1 - gamma) Delta_i^(k-1) + gamma s_i^k,
  ##   I_k       = (1/N) sum_i Delta_i^k (Delta_i^k)'
  ##
  ## I_k is a sum of outer products, so it is positive semi-definite by
  ## construction and needs no ridge to be inverted, unlike the clipped
  ## per-coordinate scaling it replaces.
  ## The preconditioner, following Gu and Kong (1998).
  ##
  ## Their recursion is
  ##
  ##   Gamma_k = Gamma_{k-1} + gamma_k ( Ibar(theta_{k-1}, X_k) - Gamma_{k-1} )
  ##   theta_k = theta_{k-1} + gamma_k Gamma_k^-1 Hbar(theta_{k-1}, X_k)
  ##
  ## and the content is in what Gamma averages. It is not the score. Their
  ## Lemma 1, a generalisation of Louis' missing-information theorem, asks for
  ## a function whose posterior expectation is the derivative of the mean
  ## field, and gives
  ##
  ##   I(theta, x) = -grad^2 log f(theta, x) - s(theta, x) s(theta, x)'.
  ##
  ## That is an information: it varies slowly with theta, so averaging it
  ## along a trajectory that is still moving is legitimate. Averaging the
  ## scores and squaring them instead -- which is what this did first -- is
  ## not, because a score is large far from the optimum and vanishes at it, so
  ## the metric came out inflated exactly where the recursion still had
  ## ground to cover, and the steps were smallest exactly where they needed to
  ## be largest. Measured against plain saemix on identical data that showed
  ## as a correlation of 0.34 against 0.46, and no gain schedule repaired it:
  ## raising the gain from a sum of 1.8 to 7.3 moved it only to 0.40, because
  ## the fault was in the matrix and not in the step length.
  ##
  ## The matrix itself is built from the per-subject scores alone, so it
  ## costs nothing beyond the score the recursion already computes.
  ## Tracing is separated from stepping on purpose. Judged only by the fits it
  ## drives, a broken metric is indistinguishable from a metric that is fine
  ## but mistuned, because once the recursion has been pushed somewhere absurd
  ## every quantity measured there is absurd too. With the trace on and the
  ## metric off, the pieces are measured along a trajectory that is known to
  ## converge, which is the only place their agreement means anything.
  if (fisherInformation && !is.null(subject) &&
      !is.null(analytic) &&
      (!is.null(analytic$internalRows) || !is.null(analytic$nativeRows))) {
    subject <- as.integer(subject)
    ## The recursion steps in its own internal coordinates, so the metric is
    ## built there too.
    internalRows <- if (!is.null(analytic$internalRows))
      analytic$internalRows else sweep(analytic$nativeRows, 2L,
      copulaScoreInternalDerivative(current, layout$lower, layout$upper), "*")
    perSubject <- rowsum(internalRows, subject, reorder = TRUE)
    subjects <- nrow(perSubject)
    ## Per-subject scores: the rows carry weights summing to one over the whole
    ## batch, so scaling by the number of subjects turns a share into a score.
    score <- perSubject * subjects

    ## Baey, Delattre, Kuhn, Leger and Lemler (2023), Algorithm 1, with the
    ## Fisher information estimator of Delattre and Kuhn (2023). Written in
    ## their notation:
    ##
    ##   Delta_i^k = (1 - gamma_k) Delta_i^{k-1} + gamma_k s_i^k
    ##   Istar_k   = (1/N) sum_i Delta_i^k (Delta_i^k)'
    ##   theta_k   = theta_{k-1} + gamma_k I_k^{-1} v_k
    ##
    ## Delta follows Baey et al. exactly, driven by the parameter gain. That
    ## coupling is not incidental: the heating phase takes full Newton steps at
    ## gamma = 1, and a full Newton step against an information that lags the
    ## iterate leaves the feasible region.
    state$deltaSubject <- if (is.null(state$deltaSubject) ||
        !identical(dim(state$deltaSubject), dim(score))) score else
      (1 - deltaGain) * state$deltaSubject + deltaGain * score
    starInformation <- crossprod(state$deltaSubject) / subjects
    state$fisherCount <- (state$fisherCount %||% 0L) + 1L

    ## The same estimator, debiased, for reporting only.
    ##
    ## Istar is an outer product of a smoothed score, so it estimates
    ##
    ##   E[Delta Delta'] = mu mu' + Var(Delta),   mu = E[s_i | y_i],
    ##
    ## and the second term does not vanish at a finite gain. Delattre and Kuhn
    ## are consistent because Var(Delta) -> 0, but the size of that term at any
    ## finite k is set by the posterior variance of the score, and that varies
    ## enormously between coordinates: a location score is linear in the latent
    ## variable, a spread score quadratic, but a shape score is cubic -- for the
    ## generalized gamma it is exactly -He_3(w)/6 -- so its posterior variance
    ## dwarfs the others when individual data are thin. Measured on the shape
    ## coordinate the inflation was large enough to put the reported information
    ## ABOVE the complete-data information, which is impossible: with the etas
    ## observed exactly the per-subject information for Q is E[He_3^2]/36 = 1/6,
    ## so SE(Q) can never fall below sqrt(6/n), and it did.
    ##
    ## Rather than estimate Var(Delta) and subtract it, smooth two copies of
    ## Delta on alternating iterations. Their Monte Carlo noises are then
    ## independent, so the cross product is unbiased for the part that matters,
    ##
    ##   E[Delta^A (Delta^B)'] = mu mu' ,
    ##
    ## with no correction term to get wrong. Each copy sees half the iterations
    ## and is individually noisier; the product is not.
    ##
    ## The copies alternate in BLOCKS, not on consecutive iterations. That is
    ## the whole difficulty: the noises cancel only to the extent that they are
    ## independent, and a Markov chain's successive states are not. Alternating
    ## every iteration removes only the fraction of the noise that decorrelates
    ## in a single step, which measurably did nothing here -- the chain mixes
    ## slowly enough that its effective sample fraction is a few per cent. The
    ## block length has to exceed the chain's autocorrelation time for the
    ## argument to hold, and 50 is chosen to sit comfortably beyond it while
    ## still giving each copy many blocks over a run of 1500 iterations.
    ##
    ## This is deliberately NOT fed back into the step. The preconditioner keeps
    ## using Istar, so the iterate sequence, and every convergence property
    ## resting on the spectrum of I_k, is bit-for-bit unchanged. A cross product
    ## is also not guaranteed positive definite in finite samples, which is
    ## harmless in a report and would not be in a metric.
    block <- 50L
    half <- if ((state$fisherCount %/% block) %% 2L) "deltaSubjectA" else
      "deltaSubjectB"
    state[[half]] <- if (is.null(state[[half]]) ||
        !identical(dim(state[[half]]), dim(score))) score else
      (1 - deltaGain) * state[[half]] + deltaGain * score
    if (!is.null(state$deltaSubjectA) && !is.null(state$deltaSubjectB) &&
        identical(dim(state$deltaSubjectA), dim(state$deltaSubjectB))) {
      cross <- crossprod(state$deltaSubjectA, state$deltaSubjectB) / subjects
      state$fisherCross <- (cross + t(cross)) / 2
    }

    ## Their stabilisation, section 3.4.2, as their code writes it:
    ##
    ##   fisher_info_mat = factor * fisher_info_mat + (1 - factor) * eye(p)
    ##
    ## so the ridge is toward the identity, `r_k = 1`. The paper offers
    ## `r_k = max(1, tr Istar_k)` as an alternative "to further avoid
    ## instabilities"; the reference implementation uses the plain identity, so
    ## that is what is used here.
    fisher <- if (isTRUE(preheating))
      gain * starInformation + (1 - gain) * diag(1, nrow(starInformation)) else
      starInformation
    ## A permanent ridge, which is what makes Assumption 3.3.3 a construction
    ## rather than an assumption.
    ##
    ## Baey et al. assume the spectrum of I_k stays in [mu_m, mu_M]. The upper
    ## bound is free here: I_k is an average of outer products of scores that
    ## are bounded on the compact internal-coordinate box, so
    ## mu_M <= sup ||s||^2. The lower bound is not free, because I_k is a sum
    ## of N rank-one terms and is singular in any direction the per-subject
    ## scores have not explored.
    ##
    ## Trying to derive it from a spanning condition on the limiting scores is
    ## circular: it uses convergence to establish the hypothesis convergence
    ## needs. Adding a fixed ridge instead gives mu_m = epsilon outright, and
    ## costs nothing that matters, because a positive definite metric cannot
    ## move the zeros of the gradient -- the stationary set of the mean ODE is
    ## unchanged. The spanning condition is then needed only for the standard
    ## errors to be finite, which is identifiability at the limit and belongs
    ## there rather than in the convergence argument.
    ##
    ## The ridge is *relative*, a fraction of the mean eigenvalue, not an
    ## absolute constant. This matters: the reason for preferring an
    ## information metric over a diagonal one is affine invariance, and adding
    ## a fixed multiple of the identity destroys precisely that -- 1e-6 means
    ## different things under different parameter scalings. A ridge of
    ## delta * tr(I)/p is equivariant and keeps the property the metric was
    ## chosen for.
    ##
    ## The size is not cosmetic. At 1e-6 a two-compartment model with three
    ## observations per subject -- four random effects, so fewer observations
    ## than random effects, and each subject's contribution to I_k rank
    ## deficient by construction -- fails outright with "natural parameter
    ## could not be mapped to a fixed reference interior". At 1e-3 and above it
    ## fits. So the ridge is a working regularisation, not a token to satisfy
    ## an assumption, and 1e-3 is the default: it rescues that case while
    ## perturbing a well-conditioned one by a few per cent, the mean eigenvalue
    ## being of order ten against a smallest eigenvalue of 0.33.
    ##
    ## Two honest limitations, which belong in the paper rather than in a
    ## footnote. First, it licenses the qualitative claim, almost-sure
    ## convergence, but not the quantitative one: Baey's bound carries 1/mu_m,
    ## so a tiny ridge makes that bound numerically vacuous. The operative
    ## constant in practice is the observed smallest eigenvalue, which is 0.33
    ## to 1.2 on the worked example, not the ridge. Second, it does not protect
    ## against the failure mode it appears to: bounding the spectrum below by
    ## mu_m bounds the step above by 1/mu_m, so a genuinely null direction is
    ## still permitted an enormous step. Protection would need a large ridge,
    ## which would shrink the metric toward plain gradient descent. This picks
    ## the faithful end of that trade-off deliberately.
    ##
    ## The reported information is left unridged; this affects only the step.
    ## Both floors, because each alone trades one condition for another. A
    ## purely relative ridge is equivariant but gives mu_m = delta tr(I)/p,
    ## which is not a fixed constant and so needs inf_k tr(I_k) > 0 -- mild,
    ## since tr(I) is the mean squared per-subject score and vanishes only if
    ## every subject's posterior mean score does, which does not happen at an
    ## interior stationary point, but still a condition. A purely absolute
    ## ridge is unconditional but not equivariant. The maximum of the two is
    ## unconditional *and* equivariant wherever there is a scale to be
    ## equivariant to; in the degenerate regime where tr(I) collapses there is
    ## no such scale, and the absolute floor takes over.
    ridge <- max(.cop$scoreMetricRidgeAbsolute %||% 1e-10,
      (.cop$scoreMetricRidge %||% 1e-3) * sum(diag(fisher)) / nrow(fisher))
    fisher <- fisher + ridge * diag(1, nrow(fisher))
    fisher <- (fisher + t(fisher)) / 2
    ## What the standard errors are built from is the *average of the
    ## information over the converging phase*, not its final value.
    ##
    ## This is what their reference implementation does and the paper's
    ## algorithm box does not show: `estim()` runs until the heating phase has
    ## ended and a convergence criterion is met, then takes ten thousand
    ## further iterations and returns the mean of `fisher_info_mat` over them,
    ## alongside the mean of theta. Reporting the last matrix instead leaves
    ## all of Delta's sampling noise in it, which inflates the information and
    ## shrinks every interval.
    ## Averaging addresses the VARIANCE of the estimate; it does not touch the
    ## bias, since the mean of many inflated matrices is still inflated. The
    ## cross product below is what removes the bias, and it is averaged the same
    ## way for the same reason.
    state$fisher <- if (!isTRUE(smoothing)) starInformation else {
      state$fisherAverageCount <- (state$fisherAverageCount %||% 0L) + 1L
      w <- 1 / state$fisherAverageCount
      if (is.null(state$fisherAverage) ||
          !identical(dim(state$fisherAverage), dim(starInformation)))
        state$fisherAverage <- starInformation else
        state$fisherAverage <- (1 - w) * state$fisherAverage +
          w * starInformation
      state$fisherAverage
    }
    if (!is.null(state$fisherCross)) {
      state$fisherDebiased <- if (!isTRUE(smoothing)) state$fisherCross else {
        w <- 1 / state$fisherAverageCount
        if (is.null(state$fisherCrossAverage) ||
            !identical(dim(state$fisherCrossAverage), dim(state$fisherCross)))
          state$fisherCrossAverage <- state$fisherCross else
          state$fisherCrossAverage <- (1 - w) * state$fisherCrossAverage +
            w * state$fisherCross
        state$fisherCrossAverage
      }
    }

    if (metricTrace) {
      smallest <- function(m) {
        e <- try(min(eigen(m, symmetric = TRUE, only.values = TRUE)$values),
          silent = TRUE)
        if (inherits(e, "try-error")) NA_real_ else e
      }
      state$louisPieces <- list(outerProduct = crossprod(score) / subjects,
        meanOuter = starInformation, metric = fisher, subjects = subjects,
        internal = current)
      state$louisTrace <- rbind(state$louisTrace, c(
        iteration = state$iteration, metric = smallest(fisher),
        empirical = smallest(starInformation), maxScore = max(abs(score))))
    }

    ## Assumption 3.3.3 asks for the spectrum of I_k to stay in [mu_m, mu_M].
    ## The pre-heating ridge supplies that early. Later, Istar_k is a sum of N
    ## rank-one terms in p dimensions and is non-singular once enough subjects
    ## have contributed, but a failed solve must not become a wild step, so an
    ## uninvertible matrix leaves the previous metric in place rather than
    ## being repaired into a different one.
    inverse <- try(chol2inv(chol(fisher)), silent = TRUE)
    if (!inherits(inverse, "try-error") && all(is.finite(inverse))) {
      state$inverseInformation <- inverse
      state$metricUpdateCount <- state$metricUpdateCount + 1L
      fisherReady <- TRUE
    }
  }
  if (!fisherReady && (is.null(state$inverseInformation) ||
      any(dim(state$inverseInformation) !=
          c(length(gradient), length(gradient))))) {
    ## Only before the first usable information estimate exists.
    state$inverseInformation <- diag(1, length(gradient))
    state$metricUpdateCount <- state$metricUpdateCount + 1L
  }
  metricEigen <- eigen(state$inverseInformation, symmetric = TRUE,
    only.values = TRUE)$values
  direction <- as.numeric(state$inverseInformation %*% gradient)
  proposal <- if (isTRUE(freeze)) current else
    current + gain * scoreScale * direction
  ## The truncation set is a box of half-width `projection` around the
  ## starting internal coordinates. That width is dimensionless -- it is
  ## generous for log, logit and angle coordinates -- but an unbounded
  ## location coordinate (the mean of a Normal covariate margin, a typical
  ## value under an identity transform) lives on the data's own scale, and a
  ## start more than `projection` units from the maximum leaves the maximum
  ## outside the box. Projected stochastic approximation then converges to a
  ## boundary point of the projected ODE, which is not a stationary point of
  ## the likelihood (Kushner and Clark 1978, Section 5.3): measured on
  ## 2026-09-05, a run started at twice the generating clearance pinned the
  ## CRP mean at start minus 24 for 2,400 iterations while the score stayed
  ## at 0.13. The remedy of Chen (2002, Section 2.3) and Andrieu, Moulines and
  ## Priouret (2005) is a nested sequence of truncation sets that expands
  ## whenever the boundary is hit: under the coercivity assumption that the
  ## convergence argument already makes, the boundary is hit finitely often
  ## almost surely, after which the recursion coincides with the unprojected
  ## one inside a compact set. Each hit doubles the half-width of the
  ## coordinates that were clipped; the number of expansions is reported, and
  ## a fit whose post-adaptation projection count is zero never touched the
  ## boundary at all. (Their construction also reinitialises the iterate on
  ## each expansion; that step is not taken here, so the guarantee is the
  ## weaker one stated in the manuscript's supplement.)
  if (is.null(state$projectionCentre) ||
      length(state$projectionCentre) != length(current)) {
    state$projectionCentre <- current
    state$projectionWidth <- rep(projection, length(current))
    state$projectionExpansions <- 0L
  }
  if (is.null(state$projectionWidth) ||
      length(state$projectionWidth) != length(current))
    state$projectionWidth <- rep(projection, length(current))
  projectionLower <- state$projectionCentre - state$projectionWidth
  projectionUpper <- state$projectionCentre + state$projectionWidth
  projected <- pmin(projectionUpper, pmax(projectionLower, proposal))
  clipped <- abs(projected - proposal) > 0
  projectionEvent <- any(clipped)
  state$projectionCount <- state$projectionCount + as.integer(projectionEvent)
  if (projectionEvent && !isTRUE(adaptMetric))
    state$postFreezeProjectionCount <- state$postFreezeProjectionCount + 1L
  if (projectionEvent) {
    state$projectionWidth[clipped] <- 2 * state$projectionWidth[clipped]
    state$projectionExpansions <- (state$projectionExpansions %||% 0L) + 1L
  }
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
  reportedState <- materializeState(reported)
  reportedSd <- if (identical(populationScale, "parameter"))
    apply(reportedState$residual[, seq_len(dEta), drop = FALSE], 2L,
      stats::sd) else copulaMarginScales(reported$margins)
  list(beta = reported$beta, margins = reported$margins,
    sd = reportedSd, vine = reported$vine,
    delta = reported$delta, residual = reported$residual,
    value = reportedValue,
    score = gradient, scoreMax = max(abs(gradient)),
    preconditionedStep = direction,
    scoreAverage = state$scoreAverage %||% rep(NA_real_, length(gradient)),
    scoreAverageMax = if (is.null(state$scoreAverage)) NA_real_ else
      max(abs(state$scoreAverage)),
    finiteDifference = h, gain = gain, scoreScale = scoreScale,
    scoreMethod = if (is.null(analytic)) "global-centered-difference" else
      if (naturalPsi) "natural-psi-centered-difference-score" else
      if (hasReference) "hybrid-fixed-reference-path-score" else
        "analytic-with-declared-local-numerical-components",
    numericalScoreComponents = if (is.null(analytic)) "all" else
      analytic$numerical,
    fisherCovariance = local({
      ## The recursion's own Fisher information, turned into the covariance of
      ## the estimate. It is accumulated in internal coordinates, so the delta
      ## method carries it back to the ones a reader can interpret: the score
      ## in internal coordinates is the native score times the derivative of
      ## the map, so the information is D I D and the covariance is
      ## D^-1 I^-1 D^-1, divided by the number of subjects.
      ##
      ## This is the estimator of Delattre and Kuhn, but averaged along the
      ## recursion rather than over one batch of posterior draws, which matters
      ## because those draws mix slowly -- an effective sample fraction of a
      ## few percent that does not improve with a longer chain. A fit visits a
      ## new MCMC state every iteration and so accumulates far more of them.
      ## Reported from the plain estimate. The debiased cross product is
      ## computed and carried in the state for comparison but is NOT used here,
      ## and that is a deliberate refusal rather than an oversight.
      ##
      ## Both estimates are known to be wrong, in opposite directions, and
      ## neither is yet trustworthy enough to report. The plain one is provably
      ## too large: it puts the reported information above the complete-data
      ## information, which no data can. The cross product comes out five to
      ## seventeen times smaller on every coordinate and indefinite, which
      ## cannot be right either -- an inflation that large would have wrecked
      ## the ordinary parameters' coverage, and it does not. The likely cause is
      ## that two copies separated by a block are centred on different theta
      ## while the iterate is still moving, which attenuates their product.
      ##
      ## Until that is resolved the honest choice is to leave the shipped
      ## behaviour alone rather than swap one unvalidated estimate for another.
      ## What the shape coordinate needs is stated in the handoff: more MCMC per
      ## iteration moves the reported error to roughly the right size.
      chosen <- state$fisher
      if (is.null(chosen) || is.null(state$deltaSubject)) NULL else {
        inverse <- try(chol2inv(chol(chosen)), silent = TRUE)
        if (inherits(inverse, "try-error")) NULL else {
          scale <- copulaScoreInternalDerivative(state$internal, layout$lower,
            layout$upper)
          answer <- sweep(sweep(inverse, 1L, scale, "*"), 2L, scale, "*") /
            nrow(state$deltaSubject)
          labels <- c(
            if (layout$nLocation) paste0("location.",
              seq_len(layout$nLocation)) else character(),
            if (layout$nMargin) unlist(lapply(seq_len(layout$d),
              function(j) if (!length(layout$marginLayout$index[[j]]))
                character() else paste0("margin.", j, ".",
                  seq_along(layout$marginLayout$index[[j]])))) else character(),
            if (layout$nEdge) paste0("correlation.",
              seq_len(layout$nEdge)) else character(),
            if (layout$nResidual) paste0("residual.",
              layout$residualFree) else character())
          if (length(labels) == nrow(answer))
            dimnames(answer) <- list(labels, labels)
          answer
        }
      }
    }),
    preconditionerMin = min(metricEigen),
    preconditionerMax = max(metricEigen), metricGain = gain,
    metricUpdateCount = state$metricUpdateCount,
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
        if (naturalPsi)
        "natural parameters held fixed; ordinary complete-data score" else
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
      projectionExpansions = state$projectionExpansions %||% 0L,
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
