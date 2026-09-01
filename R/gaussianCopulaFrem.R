## Exact collapsed likelihood helpers for Gaussian-copula FREM.
##
## The fitted population law remains joint in parameter random effects and
## covariates.  Observed covariates are conditioned on analytically and missing
## covariates are integrated out.  Thus only the PK/PD random effects have to be
## carried in the SAEM Markov state when covariates are error-free.

copulaIsFullGaussianVine <- function(vine, d = as.integer(dim(vine)["dim"])) {
  if (!inherits(vine, "vinecop_dist") || length(d) != 1L || is.na(d) || d < 1L)
    return(FALSE)
  full <- inherits(vine$structure, "rvine_structure") &&
    identical(as.integer(vine$structure$d), as.integer(d)) &&
    identical(as.integer(vine$structure$trunc_lvl), as.integer(d - 1L))
  if (!full) return(FALSE)
  flat <- copulaPadFlat(vine, d)
  all(vapply(flat, function(edge)
    edge$family %in% c("gaussian", "indep"), logical(1)))
}

copulaGaussianLogDensity <- function(z, correlation) {
  z <- as.matrix(z)
  if (!ncol(z)) return(rep(0, nrow(z)))
  correlation <- as.matrix(correlation)
  U <- try(chol(correlation), silent = TRUE)
  if (inherits(U, "try-error"))
    stop("Gaussian-copula FREM correlation matrix is not positive definite")
  standardized <- forwardsolve(t(U), t(z))
  -.5 * (ncol(z) * log(2 * pi) + 2 * sum(log(diag(U))) +
           colSums(standardized^2))
}

copulaGaussianRectangleProbability <- function(lower, upper, mean, covariance) {
  lower <- as.numeric(lower); upper <- as.numeric(upper)
  mean <- as.numeric(mean); covariance <- as.matrix(covariance)
  d <- length(lower)
  if (!d) return(1)
  if (length(upper) != d || length(mean) != d ||
      any(dim(covariance) != c(d, d)) || any(lower >= upper)) return(0)
  if (d == 1L) {
    sd <- sqrt(covariance[1L, 1L])
    return(max(0, stats::pnorm((upper - mean) / sd) -
      stats::pnorm((lower - mean) / sd)))
  }
  if (!requireNamespace("mvtnorm", quietly = TRUE))
    stop("mixed categorical Gaussian-copula likelihood requires package 'mvtnorm'")
  as.numeric(mvtnorm::pmvnorm(lower = lower, upper = upper, mean = mean,
    sigma = covariance, algorithm = mvtnorm::Miwa(steps = 128L)))
}

## Exact Gaussian-copula density/mass for mixed continuous and ordered
## categorical coordinates. Discrete coordinates contribute a conditional
## latent-normal rectangle probability, not c(F) times a probability mass.
copulaGaussianFremMixedLogDensity <- function(x, margins, correlation) {
  x <- as.matrix(x); d <- ncol(x)
  if (length(margins) != d || any(dim(correlation) != c(d, d)))
    stop("mixed Gaussian-copula dimensions are incompatible")
  continuous <- which(vapply(margins, function(m)
    identical(m$type, "continuous"), logical(1)))
  discrete <- setdiff(seq_len(d), continuous)
  zc <- matrix(numeric(), nrow(x), 0L)
  logContinuous <- numeric(nrow(x)); validContinuous <- rep(TRUE, nrow(x))
  if (length(continuous)) {
    evaluated <- copulaGaussianFremEvaluateMargins(
      x[, continuous, drop = FALSE], margins[continuous])
    validContinuous <- evaluated$valid
    zc <- evaluated$z
    logContinuous[] <- -Inf
    if (any(validContinuous)) {
      rows <- which(validContinuous)
      localZ <- zc[rows, , drop = FALSE]
      logContinuous[rows] <- copulaGaussianLogDensity(localZ,
        correlation[continuous, continuous, drop = FALSE]) -
        rowSums(stats::dnorm(localZ, log = TRUE)) +
        rowSums(evaluated$logMargin[rows, , drop = FALSE])
    }
  }
  if (!length(discrete)) return(logContinuous)

  if (length(continuous)) {
    solveContinuous <- solve(correlation[continuous, continuous, drop = FALSE])
    regression <- correlation[discrete, continuous, drop = FALSE] %*%
      solveContinuous
    conditionalMean <- matrix(NA_real_, nrow(x), length(discrete))
    if (any(validContinuous))
      conditionalMean[validContinuous, ] <-
        zc[validContinuous, , drop = FALSE] %*% t(regression)
    conditionalCovariance <- correlation[discrete, discrete, drop = FALSE] -
      regression %*% correlation[continuous, discrete, drop = FALSE]
  } else {
    conditionalMean <- matrix(0, nrow(x), length(discrete))
    conditionalCovariance <- correlation[discrete, discrete, drop = FALSE]
  }
  conditionalCovariance <- (conditionalCovariance + t(conditionalCovariance)) / 2
  lower <- upper <- matrix(NA_real_, nrow(x), length(discrete))
  for (k in seq_along(discrete)) {
    j <- discrete[k]; margin <- margins[[j]]
    upper[, k] <- stats::qnorm(margin$cdf(x[, j], margin$parameters))
    lower[, k] <- stats::qnorm(margin$cdf_left(x[, j], margin$parameters))
  }
  result <- rep(-Inf, nrow(x))
  for (i in which(validContinuous)) {
    probability <- copulaGaussianRectangleProbability(lower[i, ], upper[i, ],
      conditionalMean[i, ], conditionalCovariance)
    if (is.finite(probability) && probability > 0)
      result[i] <- logContinuous[i] + log(probability)
  }
  result
}

## Fixed-support augmentation for categorical coordinates.  If category k has
## Gaussian-score interval (a_k,b_k], introduce v in (0,1) and set
## u = F(k-) + v {F(k)-F(k-)}.  The augmented contribution is log P(k), so
## integrating v recovers the exact rectangle mass without putting a
## parameter-dependent threshold score in the Markov state.
copulaGaussianFremAugmentedEvaluateMargins <- function(
    x, margins, categoricalUniform = NULL, referenceUniform = NULL) {
  x <- as.matrix(x); d <- ncol(x)
  if (length(margins) != d)
    stop("augmented Gaussian-copula margin dimension mismatch")
  if (is.null(categoricalUniform))
    categoricalUniform <- matrix(NA_real_, nrow(x), d) else {
    categoricalUniform <- as.matrix(categoricalUniform)
    if (any(dim(categoricalUniform) != c(nrow(x), d)))
      stop("categorical augmentation must align with all population coordinates")
  }
  if (is.null(referenceUniform))
    referenceUniform <- matrix(NA_real_, nrow(x), d) else {
    referenceUniform <- as.matrix(referenceUniform)
    if (any(dim(referenceUniform) != c(nrow(x), d)))
      stop("fixed-reference augmentation must align with all population coordinates")
  }
  continuous <- which(vapply(margins, function(m)
    identical(m$type, "continuous"), logical(1)))
  discrete <- setdiff(seq_len(d), continuous)
  z <- matrix(NA_real_, nrow(x), d)
  logMargin <- matrix(NA_real_, nrow(x), d)
  valid <- rep(TRUE, nrow(x))
  for (j in continuous) {
    reference <- is.finite(referenceUniform[, j])
    if (any(reference)) {
      u <- referenceUniform[reference, j]
      ok <- u > 0 & u < 1
      rows <- which(reference)
      if (any(ok)) {
        z[rows[ok], j] <- stats::qnorm(u[ok])
        ## Density with respect to the fixed reference du is uniform.
        logMargin[rows[ok], j] <- 0
      }
      valid[rows] <- valid[rows] & ok
    }
    direct <- !reference
    if (any(direct)) {
      evaluated <- copulaGaussianFremEvaluateMargins(
        matrix(x[direct, j], ncol = 1L), list(margins[[j]]))
      z[direct, j] <- evaluated$z[, 1L]
      logMargin[direct, j] <- evaluated$logMargin[, 1L]
      valid[direct] <- valid[direct] & evaluated$valid
    }
  }
  for (j in discrete) {
    margin <- margins[[j]]; value <- x[, j]
    lo <- margin$cdf_left(value, margin$parameters)
    hi <- margin$cdf(value, margin$parameters)
    v <- categoricalUniform[, j]
    ok <- is.finite(lo) & is.finite(hi) & hi > lo &
      is.finite(v) & v > 0 & v < 1
    u <- rep(NA_real_, nrow(x))
    u[ok] <- lo[ok] + v[ok] * (hi[ok] - lo[ok])
    ok <- ok & is.finite(u) & u > 0 & u < 1
    if (any(ok)) {
      z[ok, j] <- stats::qnorm(u[ok])
      logMargin[ok, j] <- log(hi[ok] - lo[ok])
    }
    valid <- valid & ok
  }
  list(z = z, logMargin = logMargin, valid = valid,
    categorical = discrete,
    fixedReference = which(colSums(is.finite(referenceUniform)) > 0L))
}

## Draw the complete conditioning augmentation given a current eta draw and
## the observed covariates.  Observed discrete coordinates are truncated;
## missing continuous and discrete coordinates are sampled jointly.  Rejection
## sampling is exact and fails closed when a rare rectangle exceeds the cap.
copulaGaussianFremAugmentMixedConditioning <- function(
    eta, conditioning, vine, margins, dEta, rejectionMax = 100000L) {
  eta <- as.matrix(eta); conditioning <- as.matrix(conditioning)
  dEta <- as.integer(dEta); d <- length(margins)
  if (ncol(eta) != dEta || nrow(conditioning) != nrow(eta) ||
      ncol(conditioning) != d - dEta ||
      !copulaIsFullGaussianVine(vine, d))
    stop("invalid mixed-conditioning augmentation")
  if (any(vapply(margins[seq_len(dEta)], function(m)
      !identical(m$type, "continuous"), logical(1))))
    stop("mixed-conditioning augmentation requires continuous eta margins")
  R <- copulaGaussianRvineCor(vine, d)
  etaEvaluated <- copulaGaussianFremEvaluateMargins(
    eta, margins[seq_len(dEta)])
  if (!all(etaEvaluated$valid))
    stop("eta lies outside its declared margin support")
  answer <- conditioning
  uniform <- matrix(NA_real_, nrow(eta), d)
  conditioningGlobal <- dEta + seq_len(d - dEta)
  discreteGlobal <- conditioningGlobal[vapply(margins[conditioningGlobal],
    function(m) identical(m$type, "discrete"), logical(1))]

  for (i in seq_len(nrow(eta))) {
    observed <- !is.na(conditioning[i, ])
    observedContinuousLocal <- which(observed & vapply(
      margins[conditioningGlobal], function(m)
        identical(m$type, "continuous"), logical(1)))
    observedContinuousGlobal <- dEta + observedContinuousLocal
    knownGlobal <- c(seq_len(dEta), observedContinuousGlobal)
    knownZ <- etaEvaluated$z[i, ]
    if (length(observedContinuousGlobal)) {
      evaluated <- copulaGaussianFremEvaluateMargins(
        conditioning[i, observedContinuousLocal, drop = FALSE],
        margins[observedContinuousGlobal])
      if (!all(evaluated$valid))
        stop("observed continuous covariate lies outside its margin support")
      knownZ <- c(knownZ, evaluated$z[1L, ])
    }
    unknownGlobal <- setdiff(conditioningGlobal, observedContinuousGlobal)
    if (length(unknownGlobal)) {
      regression <- R[unknownGlobal, knownGlobal, drop = FALSE] %*%
        solve(R[knownGlobal, knownGlobal, drop = FALSE])
      meanUnknown <- as.numeric(regression %*% knownZ)
      covarianceUnknown <- R[unknownGlobal, unknownGlobal, drop = FALSE] -
        regression %*% R[knownGlobal, unknownGlobal, drop = FALSE]
      covarianceUnknown <- (covarianceUnknown + t(covarianceUnknown)) / 2
      lower <- rep(-Inf, length(unknownGlobal))
      upper <- rep(Inf, length(unknownGlobal))
      for (k in seq_along(unknownGlobal)) {
        j <- unknownGlobal[k]; local <- j - dEta
        if (observed[local] && identical(margins[[j]]$type, "discrete")) {
          value <- conditioning[i, local]; margin <- margins[[j]]
          lower[k] <- stats::qnorm(margin$cdf_left(value, margin$parameters))
          upper[k] <- stats::qnorm(margin$cdf(value, margin$parameters))
        }
      }
      U <- chol(covarianceUnknown); accepted <- FALSE
      for (attempt in seq_len(as.integer(rejectionMax))) {
        zUnknown <- meanUnknown +
          as.numeric(stats::rnorm(length(unknownGlobal)) %*% U)
        if (all(zUnknown > lower & zUnknown <= upper)) {
          accepted <- TRUE; break
        }
      }
      if (!accepted)
        stop("mixed categorical augmentation rejection sampler exhausted its cap")
      for (k in seq_along(unknownGlobal)) {
        j <- unknownGlobal[k]; local <- j - dEta
        if (is.na(conditioning[i, local]))
          answer[i, local] <- margins[[j]]$quantile(
            stats::pnorm(zUnknown[k]), margins[[j]]$parameters)
      }
      latent <- setNames(as.list(zUnknown), as.character(unknownGlobal))
    } else latent <- list()

    for (j in discreteGlobal) {
      local <- j - dEta
      z <- latent[[as.character(j)]]
      if (is.null(z))
        stop("internal error: categorical latent score was not augmented")
      margin <- margins[[j]]; value <- answer[i, local]
      lo <- margin$cdf_left(value, margin$parameters)
      hi <- margin$cdf(value, margin$parameters)
      v <- (stats::pnorm(z) - lo) / (hi - lo)
      if (!is.finite(v) || v <= 0 || v >= 1)
        stop("categorical latent-uniform augmentation left its fixed support")
      uniform[i, j] <- v
    }
  }
  if (anyNA(answer) || any(!is.finite(suppressWarnings(as.numeric(answer)))))
    stop("mixed-conditioning augmentation did not complete every covariate")
  list(conditioning = answer, categoricalUniform = uniform,
    method = "exact-fixed-support-latent-uniform")
}

copulaGaussianFremEvaluateMargins <- function(x, margins) {
  x <- as.matrix(x)
  if (ncol(x) != length(margins))
    stop("Gaussian-copula FREM margin dimension mismatch")
  z <- matrix(NA_real_, nrow(x), ncol(x))
  logMargin <- matrix(NA_real_, nrow(x), ncol(x))
  valid <- rep(TRUE, nrow(x))
  for (j in seq_len(ncol(x))) {
    margin <- margins[[j]]
    if (!identical(margin$type, "continuous"))
      stop("Gaussian-copula FREM currently requires continuous margins")
    values <- x[, j]
    if (identical(margin$name, "normal")) {
      mean <- if ("mean" %in% names(margin$parameters))
        unname(margin$parameters["mean"]) else 0
      sd <- unname(margin$parameters["sd"])
      localZ <- (values - mean) / sd
      density <- stats::dnorm(localZ, log = TRUE) - log(sd)
      ok <- is.finite(localZ) & is.finite(density)
    } else if (identical(margin$name, "lognormal")) {
      meanlog <- unname(margin$parameters["meanlog"])
      sdlog <- unname(margin$parameters["sdlog"])
      ok <- is.finite(values) & values > 0
      localZ <- rep(NA_real_, length(values))
      localZ[ok] <- (log(values[ok]) - meanlog) / sdlog
      density <- rep(-Inf, length(values))
      density[ok] <- stats::dnorm(localZ[ok], log = TRUE) -
        log(sdlog) - log(values[ok])
      ok <- ok & is.finite(localZ) & is.finite(density)
    } else {
      density <- margin$log_density(values, margin$parameters)
      probability <- margin$cdf(values, margin$parameters)
      ok <- is.finite(density) & is.finite(probability) &
        probability > 0 & probability < 1
      localZ <- rep(NA_real_, length(values))
      localZ[ok] <- stats::qnorm(probability[ok])
    }
    valid <- valid & ok
    if (any(ok)) {
      z[ok, j] <- localZ[ok]
      logMargin[ok, j] <- density[ok]
    }
  }
  list(z = z, logMargin = logMargin, valid = valid)
}

## Observed-data population density for rows (eta, covariates). NA is permitted
## only in covariate columns.  Missing covariates are marginalized exactly.
copulaGaussianFremLogPrior <- function(E, vine, margins, dEta,
                                       likelihoodTarget = c("joint", "conditional")) {
  likelihoodTarget <- match.arg(likelihoodTarget)
  E <- as.matrix(E)
  d <- ncol(E); dEta <- as.integer(dEta)
  if (length(margins) != d || dEta < 1L || dEta > d)
    stop("Gaussian-copula FREM dimensions are incompatible")
  if (!copulaIsFullGaussianVine(vine, d))
    stop("Gaussian-copula FREM requires a full all-Gaussian R-vine")
  if (anyNA(E[, seq_len(dEta), drop = FALSE]))
    stop("parameter random effects cannot be missing")
  conditioning <- if (dEta < d) E[, dEta + seq_len(d - dEta), drop = FALSE] else
    matrix(numeric(), nrow(E), 0L)
  if (any(is.infinite(conditioning)))
    stop("conditioning covariates may be missing but cannot be infinite")

  R <- copulaGaussianRvineCor(vine, d)
  result <- rep(-Inf, nrow(E))
  observedPattern <- if (ncol(conditioning))
    apply(!is.na(conditioning), 1L, paste0, collapse = "") else rep("", nrow(E))

  for (pattern in unique(observedPattern)) {
    rows <- which(observedPattern == pattern)
    observedCov <- if (!ncol(conditioning)) integer() else
      which(!is.na(conditioning[rows[1L], ]))
    index <- c(seq_len(dEta), dEta + observedCov)
    logJoint <- copulaGaussianFremMixedLogDensity(
      E[rows, index, drop = FALSE], margins[index],
      R[index, index, drop = FALSE])
    valid <- is.finite(logJoint); okRows <- rows[valid]
    if (!length(okRows)) next
    logJoint <- logJoint[valid]

    if (identical(likelihoodTarget, "conditional") && length(observedCov)) {
      covIndex <- dEta + observedCov
      logCov <- copulaGaussianFremMixedLogDensity(
        E[okRows, covIndex, drop = FALSE], margins[covIndex],
        R[covIndex, covIndex, drop = FALSE])
      logJoint <- logJoint - logCov
    }
    result[okRows] <- logJoint
  }
  attr(result, "pit_clipped") <- 0L
  attr(result, "density_floored") <- 0L
  result
}

## Marginal density of whichever conditioning coordinates are observed.  This
## is the p(c_obs) term needed to turn an exact conditional response likelihood
## into the joint FREM likelihood; missing coordinates are integrated out by
## selecting the corresponding Gaussian-score submatrix.
copulaGaussianFremConditioningLogDensity <- function(
    conditioning, vine, margins, dEta) {
  conditioning <- as.matrix(conditioning); dEta <- as.integer(dEta)
  if (!ncol(conditioning)) return(rep(0, nrow(conditioning)))
  d <- dEta + ncol(conditioning)
  if (length(margins) != d || !copulaIsFullGaussianVine(vine, d))
    stop("invalid Gaussian FREM conditioning density specification")
  R <- copulaGaussianRvineCor(vine, d); result <- rep(-Inf, nrow(conditioning))
  pattern <- apply(!is.na(conditioning), 1L, paste0, collapse = "")
  for (key in unique(pattern)) {
    rows <- which(pattern == key)
    observed <- which(!is.na(conditioning[rows[1L], ]))
    if (!length(observed)) { result[rows] <- 0; next }
    index <- dEta + observed
    density <- copulaGaussianFremMixedLogDensity(
      conditioning[rows, observed, drop = FALSE], margins[index],
      R[index, index, drop = FALSE])
    valid <- is.finite(density); result[rows[valid]] <- density[valid]
  }
  result
}

## Conditional normal-score parameters for z_eta | observed covariates. Missing
## covariates are marginalized, so each missingness pattern gets its own Schur
## complement. The return value is deliberately pattern-indexed for reuse.
copulaGaussianFremConditional <- function(conditioning, vine, margins, dEta) {
  conditioning <- as.matrix(conditioning)
  dEta <- as.integer(dEta); d <- dEta + ncol(conditioning)
  if (length(margins) != d || !copulaIsFullGaussianVine(vine, d))
    stop("invalid Gaussian-copula FREM conditional specification")
  if (any(is.infinite(conditioning)))
    stop("conditioning covariates may be missing but cannot be infinite")
  R <- copulaGaussianRvineCor(vine, d)
  etaIndex <- seq_len(dEta)
  mean <- matrix(0, nrow(conditioning), dEta)
  covariance <- vector("list", nrow(conditioning))
  pattern <- apply(!is.na(conditioning), 1L, paste0, collapse = "")
  patternState <- vector("list", length(unique(pattern)))
  names(patternState) <- unique(pattern)

  for (key in names(patternState)) {
    rows <- which(pattern == key)
    observed <- which(!is.na(conditioning[rows[1L], ]))
    if (!length(observed)) {
      conditionalCovariance <- R[etaIndex, etaIndex, drop = FALSE]
      conditionalMean <- matrix(0, length(rows), dEta)
    } else {
      observedIndex <- dEta + observed
      evaluated <- copulaGaussianFremEvaluateMargins(
        conditioning[rows, observed, drop = FALSE], margins[observedIndex])
      if (!all(evaluated$valid))
        stop("an observed conditioning value lies outside its declared margin support")
      solveObserved <- solve(R[observedIndex, observedIndex, drop = FALSE])
      regression <- R[etaIndex, observedIndex, drop = FALSE] %*% solveObserved
      conditionalMean <- evaluated$z %*% t(regression)
      conditionalCovariance <- R[etaIndex, etaIndex, drop = FALSE] -
        regression %*% R[observedIndex, etaIndex, drop = FALSE]
      conditionalCovariance <- (conditionalCovariance +
                                  t(conditionalCovariance)) / 2
    }
    mean[rows, ] <- conditionalMean
    covariance[rows] <- rep(list(conditionalCovariance), length(rows))
    patternState[[key]] <- list(rows = rows, observed = observed,
                                covariance = conditionalCovariance)
  }
  list(mean = mean, covariance = covariance, pattern = pattern,
       patternState = patternState)
}

## Exact latent augmentation for missing continuous Gaussian-copula FREM
## covariates. Given an accepted eta draw, sample the missing covariate normal
## scores from p(z_mis | z_eta, z_cov,obs), then map back through their declared
## margins. The resulting (eta,c_complete) row is a draw from the same joint
## posterior used by SAEM. For identity/log-Gaussianizable margins this restores
## the complete Gaussian sufficient statistic and its closed-form M-step rather
## than invoking a nested missing-pattern EM inside the M-step.
copulaGaussianFremImputeMissingConditioning <- function(
    eta, conditioning, vine, margins, dEta) {
  eta <- as.matrix(eta); conditioning <- as.matrix(conditioning)
  dEta <- as.integer(dEta); d <- dEta + ncol(conditioning)
  if (ncol(eta) != dEta || nrow(eta) != nrow(conditioning) ||
      length(margins) != d || !copulaIsFullGaussianVine(vine, d))
    stop("invalid missing-covariate augmentation dimensions")
  if (!anyNA(conditioning)) return(conditioning)
  if (any(vapply(margins, function(m)
      !identical(m$type, "continuous"), logical(1))))
    stop("missing-covariate Gaussian augmentation currently requires continuous margins")
  etaEvaluated <- copulaGaussianFremEvaluateMargins(
    eta, margins[seq_len(dEta)])
  if (!all(etaEvaluated$valid))
    stop("eta draw lies outside its declared margin during covariate augmentation")
  R <- copulaGaussianRvineCor(vine, d); answer <- conditioning
  pattern <- apply(is.na(conditioning), 1L, paste0, collapse = "")
  for (key in unique(pattern)) {
    rows <- which(pattern == key)
    missingLocal <- which(is.na(conditioning[rows[1L], ]))
    if (!length(missingLocal)) next
    observedLocal <- setdiff(seq_len(ncol(conditioning)), missingLocal)
    knownGlobal <- c(seq_len(dEta), dEta + observedLocal)
    missingGlobal <- dEta + missingLocal
    zKnown <- etaEvaluated$z[rows, , drop = FALSE]
    if (length(observedLocal)) {
      observed <- copulaGaussianFremEvaluateMargins(
        conditioning[rows, observedLocal, drop = FALSE],
        margins[dEta + observedLocal])
      if (!all(observed$valid))
        stop("observed covariate lies outside its margin during augmentation")
      zKnown <- cbind(zKnown, observed$z)
    }
    regression <- R[missingGlobal, knownGlobal, drop = FALSE] %*%
      solve(R[knownGlobal, knownGlobal, drop = FALSE])
    conditionalMean <- zKnown %*% t(regression)
    conditionalCovariance <- R[missingGlobal, missingGlobal, drop = FALSE] -
      regression %*% R[knownGlobal, missingGlobal, drop = FALSE]
    conditionalCovariance <- (conditionalCovariance +
      t(conditionalCovariance)) / 2
    U <- chol(conditionalCovariance)
    zMissing <- conditionalMean +
      matrix(stats::rnorm(length(rows) * length(missingLocal)),
        nrow = length(rows)) %*% U
    for (k in seq_along(missingLocal)) {
      margin <- margins[[missingGlobal[k]]]
      answer[rows, missingLocal[k]] <- margin$quantile(
        stats::pnorm(zMissing[, k]), margin$parameters)
    }
  }
  if (any(!is.finite(answer)))
    stop("missing-covariate augmentation produced a non-finite value")
  answer
}

## Build the exact negative log density of eta | observed covariates once per
## E-step.  The covariate PITs, Schur complements and Cholesky factors are
## invariant across all random-walk proposals in that E-step.  Using this
## conditional density in an MH ratio is identical to using the joint density:
## the omitted p(c_obs) is constant for each subject/chain row and cancels.
copulaGaussianFremConditionalUetaEvaluator <- function(
    conditioning, vine, margins, dEta, conditional = NULL) {
  conditioning <- as.matrix(conditioning); dEta <- as.integer(dEta)
  if (any(vapply(margins[dEta + seq_len(length(margins) - dEta)],
      function(m) identical(m$type, "discrete"), logical(1))))
    return(copulaGaussianFremCategoricalKernel(
      conditioning, vine, margins, dEta)$negative)
  if (!copulaIsFullGaussianVine(vine, length(margins)) ||
      ncol(conditioning) != length(margins) - dEta)
    stop("invalid conditional-prior evaluator specification")
  if (is.null(conditional)) conditional <- copulaGaussianFremConditional(
    conditioning, vine, margins, dEta)
  etaMargins <- margins[seq_len(dEta)]
  allNormal <- all(vapply(etaMargins, function(m)
    identical(m$name, "normal") && isTRUE(m$centered) &&
      identical(names(m$parameters), "sd"), logical(1)))
  patternState <- conditional$patternState

  if (allNormal) {
    scale <- vapply(etaMargins, function(m)
      unname(m$parameters[["sd"]]), numeric(1))
    meanEta <- sweep(conditional$mean, 2L, scale, "*")
    for (key in names(patternState)) {
      state <- patternState[[key]]
      covariance <- state$covariance * tcrossprod(scale)
      U <- chol(covariance)
      state$U <- U
      state$logNormalizer <- .5 * dEta * log(2 * pi) + sum(log(diag(U)))
      patternState[[key]] <- state
    }
    return(function(eta) {
      eta <- as.matrix(eta)
      if (nrow(eta) != nrow(conditioning) || ncol(eta) != dEta)
        stop("conditional eta evaluator received incompatible dimensions")
      answer <- numeric(nrow(eta))
      for (state in patternState) {
        rows <- state$rows
        difference <- eta[rows, , drop = FALSE] -
          meanEta[rows, , drop = FALSE]
        standardized <- forwardsolve(t(state$U), t(difference))
        answer[rows] <- state$logNormalizer + .5 * colSums(standardized^2)
      }
      answer
    })
  }

  ## Arbitrary continuous eta margins retain their exact Jacobian correction
  ## from native eta coordinates to Gaussian scores.
  function(eta) {
    eta <- as.matrix(eta)
    if (nrow(eta) != nrow(conditioning) || ncol(eta) != dEta)
      stop("conditional eta evaluator received incompatible dimensions")
    evaluated <- copulaGaussianFremEvaluateMargins(eta, etaMargins)
    if (!all(evaluated$valid))
      stop("conditional eta proposal lies outside its declared margin support")
    answer <- numeric(nrow(eta))
    correction <- rowSums(evaluated$logMargin) -
      rowSums(stats::dnorm(evaluated$z, log = TRUE))
    for (state in patternState) {
      rows <- state$rows
      centred <- evaluated$z[rows, , drop = FALSE] -
        conditional$mean[rows, , drop = FALSE]
      answer[rows] <- -(copulaGaussianLogDensity(
        centred, state$covariance) + correction[rows])
    }
    answer
  }
}

copulaGaussianFremConditionalKernel <- function(
    conditioning, vine, margins, dEta) {
  conditioning <- as.matrix(conditioning); dEta <- as.integer(dEta)
  if (any(vapply(margins[dEta + seq_len(length(margins) - dEta)],
      function(m) identical(m$type, "discrete"), logical(1))))
    return(copulaGaussianFremCategoricalKernel(
      conditioning, vine, margins, dEta))
  conditional <- copulaGaussianFremConditional(
    conditioning, vine, margins, dEta)
  patternState <- conditional$patternState
  for (key in names(patternState)) {
    patternState[[key]]$chol <- chol(patternState[[key]]$covariance)
  }
  etaMargins <- margins[seq_len(dEta)]
  allNormal <- all(vapply(etaMargins, function(m)
    identical(m$name, "normal") && isTRUE(m$centered) &&
      identical(names(m$parameters), "sd"), logical(1)))
  scale <- if (allNormal) vapply(etaMargins, function(m)
    unname(m$parameters[["sd"]]), numeric(1)) else NULL
  negative <- copulaGaussianFremConditionalUetaEvaluator(
    conditioning, vine, margins, dEta, conditional)
  random <- function() {
    z <- conditional$mean
    for (state in patternState) {
      rows <- state$rows
      noise <- matrix(stats::rnorm(length(rows) * dEta), ncol = dEta)
      z[rows, ] <- z[rows, , drop = FALSE] + noise %*% state$chol
    }
    if (allNormal) sweep(z, 2L, scale, "*") else
      copulaMarginsQuantile(stats::pnorm(z), etaMargins)
  }
  list(negative = negative, random = random,
    conditional = conditional, patternState = patternState)
}

copulaGaussianFremCategoricalKernel <- function(
    conditioning, vine, margins, dEta, rejectionMax = 100000L) {
  conditioning <- as.matrix(conditioning); dEta <- as.integer(dEta)
  d <- length(margins); dConditioning <- d - dEta
  if (ncol(conditioning) != dConditioning ||
      !copulaIsFullGaussianVine(vine, d) ||
      any(vapply(margins[seq_len(dEta)], function(m)
        !identical(m$type, "continuous"), logical(1))))
    stop("invalid categorical Gaussian-copula FREM kernel")
  R <- copulaGaussianRvineCor(vine, d); etaIndex <- seq_len(dEta)
  etaMargins <- margins[etaIndex]
  allNormal <- all(vapply(etaMargins, function(m)
    identical(m$name, "normal") && isTRUE(m$centered) &&
      identical(names(m$parameters), "sd"), logical(1)))
  etaScale <- if (allNormal) vapply(etaMargins, function(m)
    unname(m$parameters[["sd"]]), numeric(1)) else NULL
  pattern <- apply(!is.na(conditioning), 1L, paste0, collapse = "")
  states <- vector("list", length(unique(pattern))); names(states) <- unique(pattern)

  for (key in names(states)) {
    rows <- which(pattern == key)
    observedLocal <- which(!is.na(conditioning[rows[1L], ]))
    observedGlobal <- dEta + observedLocal
    discreteGlobal <- observedGlobal[vapply(margins[observedGlobal], function(m)
      identical(m$type, "discrete"), logical(1))]
    continuousGlobal <- setdiff(observedGlobal, discreteGlobal)
    zContinuous <- matrix(numeric(), length(rows), 0L)
    if (length(continuousGlobal)) {
      local <- continuousGlobal - dEta
      evaluated <- copulaGaussianFremEvaluateMargins(
        conditioning[rows, local, drop = FALSE], margins[continuousGlobal])
      if (!all(evaluated$valid))
        stop("an observed continuous covariate lies outside its margin support")
      zContinuous <- evaluated$z
    }
    jointIndex <- c(etaIndex, discreteGlobal)
    if (length(continuousGlobal)) {
      inverseContinuous <- solve(R[continuousGlobal, continuousGlobal, drop = FALSE])
      regression <- R[jointIndex, continuousGlobal, drop = FALSE] %*%
        inverseContinuous
      meanJoint <- zContinuous %*% t(regression)
      covarianceJoint <- R[jointIndex, jointIndex, drop = FALSE] -
        regression %*% R[continuousGlobal, jointIndex, drop = FALSE]
    } else {
      meanJoint <- matrix(0, length(rows), length(jointIndex))
      covarianceJoint <- R[jointIndex, jointIndex, drop = FALSE]
    }
    covarianceJoint <- (covarianceJoint + t(covarianceJoint)) / 2
    nd <- length(discreteGlobal); discretePosition <- dEta + seq_len(nd)
    lower <- upper <- matrix(numeric(), length(rows), nd)
    if (nd) for (k in seq_len(nd)) {
      j <- discreteGlobal[k]; local <- j - dEta; margin <- margins[[j]]
      value <- conditioning[rows, local]
      upper[, k] <- stats::qnorm(margin$cdf(value, margin$parameters))
      lower[, k] <- stats::qnorm(margin$cdf_left(value, margin$parameters))
      if (any(lower[, k] >= upper[, k]))
        stop("categorical observation has zero declared probability")
    }
    covEE <- covarianceJoint[etaIndex, etaIndex, drop = FALSE]
    state <- list(rows = rows, discrete = discreteGlobal,
      meanE = meanJoint[, etaIndex, drop = FALSE], covEE = covEE,
      cholEE = chol(covEE), lower = lower, upper = upper)
    if (nd) {
      covDD <- covarianceJoint[discretePosition, discretePosition, drop = FALSE]
      covED <- covarianceJoint[etaIndex, discretePosition, drop = FALSE]
      state$meanD <- meanJoint[, discretePosition, drop = FALSE]
      state$covDD <- covDD; state$cholDD <- chol(covDD)
      state$regressionEonD <- covED %*% solve(covDD)
      state$covEgivenD <- covEE - state$regressionEonD %*% t(covED)
      state$covEgivenD <- (state$covEgivenD + t(state$covEgivenD)) / 2
      state$cholEgivenD <- chol(state$covEgivenD)
      state$regressionDonE <- t(covED) %*% solve(covEE)
      state$covDgivenE <- covDD - state$regressionDonE %*% covED
      state$covDgivenE <- (state$covDgivenE + t(state$covDgivenE)) / 2
    }
    states[[key]] <- state
  }

  negative <- function(eta) {
    eta <- as.matrix(eta)
    if (nrow(eta) != nrow(conditioning) || ncol(eta) != dEta)
      stop("categorical conditional evaluator received incompatible eta")
    evaluated <- copulaGaussianFremEvaluateMargins(eta, etaMargins)
    if (!all(evaluated$valid)) return(rep(Inf, nrow(eta)))
    correction <- rowSums(evaluated$logMargin) -
      rowSums(stats::dnorm(evaluated$z, log = TRUE))
    result <- rep(Inf, nrow(eta))
    for (state in states) {
      rows <- state$rows; zeta <- evaluated$z[rows, , drop = FALSE]
      logEta <- copulaGaussianLogDensity(
        zeta - state$meanE, state$covEE) + correction[rows]
      if (!length(state$discrete)) { result[rows] <- -logEta; next }
      meanDgivenE <- state$meanD +
        (zeta - state$meanE) %*% t(state$regressionDonE)
      probability <- vapply(seq_along(rows), function(k)
        copulaGaussianRectangleProbability(state$lower[k, ], state$upper[k, ],
          meanDgivenE[k, ], state$covDgivenE), numeric(1))
      valid <- probability > 0 & is.finite(probability)
      result[rows[valid]] <- -(logEta[valid] + log(probability[valid]))
    }
    result
  }

  random <- function() {
    zeta <- matrix(NA_real_, nrow(conditioning), dEta)
    for (state in states) {
      rows <- state$rows; nr <- length(rows); nd <- length(state$discrete)
      if (!nd) {
        zeta[rows, ] <- state$meanE +
          matrix(stats::rnorm(nr * dEta), ncol = dEta) %*% state$cholEE
        next
      }
      zd <- matrix(NA_real_, nr, nd)
      for (k in seq_len(nr)) {
        if (nd == 1L) {
          sd <- sqrt(state$covDD[1L, 1L]); mu <- state$meanD[k, 1L]
          pa <- stats::pnorm((state$lower[k, 1L] - mu) / sd)
          pb <- stats::pnorm((state$upper[k, 1L] - mu) / sd)
          if (!is.finite(pa + pb) || pb <= pa)
            stop("categorical truncation interval has zero numerical probability")
          zd[k, 1L] <- mu + sd * stats::qnorm(pa + stats::runif(1L) * (pb - pa))
        } else {
          accepted <- FALSE
          for (attempt in seq_len(rejectionMax)) {
            proposal <- state$meanD[k, ] +
              as.numeric(stats::rnorm(nd) %*% state$cholDD)
            if (all(proposal > state$lower[k, ] & proposal <= state$upper[k, ])) {
              zd[k, ] <- proposal; accepted <- TRUE; break
            }
          }
          if (!accepted)
            stop("categorical truncated-normal rejection sampler exhausted its cap")
        }
      }
      meanEgivenD <- state$meanE +
        (zd - state$meanD) %*% t(state$regressionEonD)
      zeta[rows, ] <- meanEgivenD +
        matrix(stats::rnorm(nr * dEta), ncol = dEta) %*% state$cholEgivenD
    }
    if (allNormal) sweep(zeta, 2L, etaScale, "*") else
      copulaMarginsQuantile(stats::pnorm(zeta), etaMargins)
  }
  list(negative = negative, random = random, states = states,
    categorical = TRUE)
}

copulaGaussianFremRandEta <- function(conditioning, vine, margins, dEta) {
  copulaGaussianFremConditionalKernel(
    conditioning, vine, margins, dEta)$random()
}
