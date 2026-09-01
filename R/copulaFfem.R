## Exact translation of a fitted Gaussian-copula FREM population model to an
## FFEM representation for any selected covariate subset.

copulaFremToFfem <- function(object, covariates = NULL,
                             representation = c("auto", "linear",
                               "distributional")) {
  representation <- match.arg(representation)
  state <- if (inherits(object, "SaemixObject")) copulaGet(object) else object
  if (!is.list(state) || is.null(state$vine) || is.null(state$margins))
    stop("object must be a fitted SaemixObject or copula state")
  dEta <- as.integer(state$dEta %||% 0L)
  dConditioning <- as.integer(state$dConditioning %||%
    (length(state$margins) - dEta))
  if (dEta < 1L || dConditioning < 1L ||
      length(state$margins) != dEta + dConditioning ||
      !copulaIsFullGaussianVine(state$vine, dEta + dConditioning))
    stop("translation requires a fitted full Gaussian-copula FREM model")

  covariateNames <- state$conditioningName
  if (is.null(covariateNames) || length(covariateNames) != dConditioning)
    covariateNames <- colnames(state$conditioning)
  if (is.null(covariateNames) || length(covariateNames) != dConditioning)
    covariateNames <- paste0("covariate", seq_len(dConditioning))
  etaNames <- state$variableName
  if (is.null(etaNames) || length(etaNames) != dEta)
    etaNames <- paste0("eta", seq_len(dEta))

  if (is.null(covariates)) selected <- seq_len(dConditioning) else
    if (is.character(covariates)) {
      selected <- match(covariates, covariateNames)
      if (anyNA(selected)) stop("unknown covariate name: ",
        paste(covariates[is.na(selected)], collapse = ", "))
    } else if (is.logical(covariates)) {
      if (length(covariates) != dConditioning)
        stop("logical covariate selector has incompatible length")
      selected <- which(covariates)
    } else selected <- as.integer(covariates)
  if (!length(selected) || anyNA(selected) || any(selected < 1L) ||
      any(selected > dConditioning) || anyDuplicated(selected))
    stop("covariates must select one or more unique fitted covariates")

  R <- copulaGaussianRvineCor(state$vine, dEta + dConditioning)
  etaIndex <- seq_len(dEta)
  covariateIndex <- dEta + selected
  inverseCovariate <- solve(R[covariateIndex, covariateIndex, drop = FALSE])
  betaScore <- R[etaIndex, covariateIndex, drop = FALSE] %*% inverseCovariate
  conditionalCorrelation <- R[etaIndex, etaIndex, drop = FALSE] -
    betaScore %*% R[covariateIndex, etaIndex, drop = FALSE]
  conditionalCorrelation <- (conditionalCorrelation +
    t(conditionalCorrelation)) / 2
  dimnames(betaScore) <- list(etaNames, covariateNames[selected])
  dimnames(conditionalCorrelation) <- list(etaNames, etaNames)
  etaMargins <- state$margins[etaIndex]
  selectedMargins <- state$margins[covariateIndex]
  normalEta <- all(vapply(etaMargins, function(m)
    identical(m$name, "normal") && isTRUE(m$centered) &&
      identical(names(m$parameters), "sd"), logical(1)))
  continuousCovariates <- all(vapply(selectedMargins, function(m)
    identical(m$type, "continuous"), logical(1)))
  linearAvailable <- normalEta && continuousCovariates
  if (identical(representation, "linear") && !linearAvailable)
    stop(paste0("a conventional FFEM requires Normal parameter margins and ",
      "continuous selected covariates; use representation='distributional'"))
  if (identical(representation, "auto")) representation <-
    if (linearAvailable) "linear" else "distributional"

  etaScale <- copulaMarginScales(etaMargins)
  coefficient <- if (normalEta)
    diag(etaScale, dEta) %*% betaScore else NULL
  omega <- if (normalEta)
    diag(etaScale, dEta) %*% conditionalCorrelation %*%
      diag(etaScale, dEta) else NULL
  if (!is.null(coefficient)) {
    dimnames(coefficient) <- list(etaNames, covariateNames[selected])
    dimnames(omega) <- list(etaNames, etaNames)
  }

  normalCovariates <- continuousCovariates && all(vapply(selectedMargins,
    function(m) identical(m$name, "normal") &&
      all(c("mean", "sd") %in% names(m$parameters)), logical(1)))
  native <- NULL
  if (normalEta && normalCovariates) {
    centre <- vapply(selectedMargins, function(m)
      unname(m$parameters[["mean"]]), numeric(1))
    scale <- vapply(selectedMargins, function(m)
      unname(m$parameters[["sd"]]), numeric(1))
    betaNative <- coefficient %*% diag(1 / scale, length(scale))
    dimnames(betaNative) <- list(etaNames, covariateNames[selected])
    names(centre) <- names(scale) <- covariateNames[selected]
    native <- list(coefficient = betaNative, centre = centre, scale = scale)
  }

  answer <- list(representation = representation,
    linearAvailable = linearAvailable, selected = selected,
    covariates = covariateNames[selected], allCovariates = covariateNames,
    eta = etaNames, betaScore = betaScore,
    conditionalCorrelation = conditionalCorrelation,
    coefficient = coefficient, omega = omega, native = native,
    state = list(vine = state$vine, margins = state$margins,
      dEta = dEta, dConditioning = dConditioning,
      conditioningName = covariateNames, variableName = etaNames))
  class(answer) <- c("saemixCopulaFfem", "list")
  answer
}

copulaFfemConditioning <- function(object, newdata) {
  if (!inherits(object, "saemixCopulaFfem"))
    stop("object must be produced by copulaFremToFfem")
  if (is.vector(newdata) && !is.list(newdata)) {
    if (length(object$selected) == 1L && is.null(names(newdata)))
      newdata <- matrix(newdata, ncol = 1L) else {
      if (length(newdata) != length(object$selected) || is.null(names(newdata)))
        stop(paste0("a vector for multiple selected covariates must contain ",
          "one named value per covariate"))
      newdata <- matrix(newdata, nrow = 1L,
        dimnames = list(NULL, names(newdata)))
    }
  } else newdata <- as.matrix(newdata)
  storage.mode(newdata) <- "double"
  if (!nrow(newdata) || ncol(newdata) != length(object$selected) ||
      any(is.na(newdata)) || any(!is.finite(newdata)))
    stop("newdata must contain one finite column per selected covariate")
  if (!is.null(colnames(newdata))) {
    order <- match(object$covariates, colnames(newdata))
    if (anyNA(order)) stop("newdata does not contain every selected covariate")
    newdata <- newdata[, order, drop = FALSE]
  }
  full <- matrix(NA_real_, nrow(newdata), object$state$dConditioning,
    dimnames = list(NULL, object$allCovariates))
  full[, object$selected] <- newdata
  full
}

copulaFfemTransform <- function(object, newdata) {
  full <- copulaFfemConditioning(object, newdata)
  margins <- object$state$margins
  dEta <- object$state$dEta
  answer <- matrix(NA_real_, nrow(full), length(object$selected),
    dimnames = list(NULL, object$covariates))
  for (k in seq_along(object$selected)) {
    margin <- margins[[dEta + object$selected[k]]]
    if (!identical(margin$type, "continuous"))
      stop("ordered categorical covariates do not have one deterministic FFEM transform")
    u <- margin$cdf(full[, object$selected[k]], margin$parameters)
    answer[, k] <- stats::qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12))
  }
  answer
}

copulaFfemLocation <- function(object, newdata) {
  if (!isTRUE(object$linearAvailable))
    stop("a linear FFEM location is unavailable; use copulaFfemSimulate")
  copulaFfemTransform(object, newdata) %*% t(object$coefficient)
}

copulaFfemQuantile <- function(object, newdata,
                               probabilities = c(.05, .5, .95)) {
  probabilities <- as.numeric(probabilities)
  if (!length(probabilities) || any(!is.finite(probabilities)) ||
      any(probabilities <= 0 | probabilities >= 1))
    stop("probabilities must lie strictly between zero and one")
  full <- copulaFfemConditioning(object, newdata)
  selectedMargins <- object$state$margins[
    object$state$dEta + object$selected]
  if (any(vapply(selectedMargins, function(m)
      !identical(m$type, "continuous"), logical(1))))
    stop("analytic FFEM quantiles require continuous selected covariates")
  conditional <- copulaGaussianFremConditional(full, object$state$vine,
    object$state$margins, object$state$dEta)
  rows <- vector("list", nrow(full) * length(object$eta) * length(probabilities))
  cursor <- 0L
  for (i in seq_len(nrow(full))) for (j in seq_along(object$eta)) {
    z <- conditional$mean[i, j] +
      sqrt(conditional$covariance[[i]][j, j]) * stats::qnorm(probabilities)
    margin <- object$state$margins[[j]]
    value <- margin$quantile(stats::pnorm(z), margin$parameters)
    for (k in seq_along(probabilities)) {
      cursor <- cursor + 1L
      rows[[cursor]] <- data.frame(row = i, eta = object$eta[j],
        probability = probabilities[k], value = value[k])
    }
  }
  do.call(rbind, rows)
}

copulaFfemSimulate <- function(object, newdata, n = 1L, seed = NULL) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n < 1L)
    stop("n must be a positive integer")
  full <- copulaFfemConditioning(object, newdata)
  sourceRow <- rep(seq_len(nrow(full)), each = n)
  full <- full[rep(seq_len(nrow(full)), each = n), , drop = FALSE]
  if (!is.null(seed)) set.seed(as.integer(seed))
  answer <- copulaGaussianFremRandEta(full, object$state$vine,
    object$state$margins, object$state$dEta)
  colnames(answer) <- object$eta
  attr(answer, "sourceRow") <- sourceRow
  answer
}

print.saemixCopulaFfem <- function(x, ...) {
  cat("Gaussian-copula FREM to FFEM translation\n")
  cat("  representation:", x$representation, "\n")
  cat("  covariates:", paste(x$covariates, collapse = ", "), "\n")
  cat("  parameter random effects:", paste(x$eta, collapse = ", "), "\n")
  if (isTRUE(x$linearAvailable)) {
    cat("  coefficient matrix on Gaussian-score covariate scale:\n")
    print(x$coefficient)
    cat("  unexplained parameter covariance:\n")
    print(x$omega)
  } else cat("  exact output: conditional quantiles/simulation\n")
  invisible(x)
}
