## Observed-data Fisher information from the score function.
##
## Delattre & Kuhn (2023, arXiv:1909.06094) estimate the Fisher information of
## a latent-variable model from first derivatives alone:
##
##   I = (1/N) sum_i E[ grad log f_i | y_i ] E[ grad log f_i | y_i ]'
##
## with f_i the complete-data density of subject i. The inner expectations are
## the same posterior expectations Fisher's identity already uses for the
## score, so the estimator needs no second derivatives, no extra likelihood
## evaluations beyond the posterior draws, and no new MCMC machinery. It is
## unbiased and consistent for the Fisher information at the true parameter.
##
## Two things follow from it. Inverted, it is the asymptotic covariance of the
## estimate, which is what a flexible-margin fit otherwise has no route to:
## the likelihood-score recursion reports a point estimate and nothing else.
## And it is the natural preconditioner for that recursion, in place of a
## per-coordinate scaling that cannot see the correlation between, say, a
## marginal spread and the copula correlation that competes with it.
##
## This is the post-hoc estimator, evaluated once at the fitted parameter. It
## is deliberately separate from the recursion: it can be checked against a
## standard fit's own standard errors before anything in the estimation loop
## depends on it.
##
## The complete-data log density of one subject splits into the population
## term and the response term,
##
##   log f_i = log p(eta_i; theta) + sum_j log p(y_ij | phi_i = mu(theta) + eta_i),
##
## and both are evaluated exactly. The derivative in theta is taken by central
## differences of that per-subject vector: the objective is a sum over
## subjects, so one pair of evaluations per coordinate yields every subject's
## derivative at once, which is why the per-subject scores cost no more than
## the aggregate one.

## The parameter vector this works in: population locations for the eta
## coordinates, the free parameters of the eta margins, the Gaussian
## correlation angles, and the free residual-error parameters. It is the
## natural parameterisation rather than the recursion's internal one, so the
## covariance needs no delta method to be interpretable.
copulaFisherParameterVector <- function(state, locations, residual,
                                        residualFree) {
  marginLayout <- copulaMarginLayout(state$margins[seq_len(state$dEta)])
  angles <- copulaGaussianCorrelationAngles(
    copulaGaussianRvineCor(state$vine, state$d))
  value <- c(locations, marginLayout$par, angles, residual[residualFree])
  names(value) <- c(
    paste0("location.", seq_along(locations)),
    if (length(marginLayout$par)) paste0("margin.",
      seq_along(marginLayout$par)) else character(),
    if (length(angles)) paste0("angle.", seq_along(angles)) else character(),
    if (length(residualFree)) paste0("residual.", residualFree) else
      character())
  list(value = value, marginLayout = marginLayout,
    nLocation = length(locations), nMargin = length(marginLayout$par),
    nAngle = length(angles), residualFree = residualFree)
}

## The per-subject complete-data log density, and the parameter vector it is a
## function of, built from a fitted object. Exposed on its own so that the
## quantity everything else rests on can be differentiated and checked against
## an independent implementation without running the estimator.
copulaFisherPerSubjectLogDensity <- function(object) {
  copulaFisherContext(object)$perSubject
}

copulaFisherContext <- function(object) {
  if (!inherits(object, "SaemixObject"))
    stop("object must be a fitted SaemixObject")
  state <- copulaGet(object)
  if (!copulaIsFullGaussianVine(state$vine, state$d))
    stop("score Fisher information requires a Gaussian-copula population")
  index <- as.integer(state$etaIndex)
  dEta <- length(index)
  transform <- as.integer(object["model"]["transform.par"][index])
  meanPhi <- object["results"]["mean.phi"][, index, drop = FALSE]
  if (any(abs(sweep(meanPhi, 2L, meanPhi[1L, ], "-")) > 1e-10))
    stop("score Fisher information does not yet cover covariate models on eta")
  locations <- as.numeric(meanPhi[1L, ])
  residual <- as.numeric(object["results"]["respar"])
  errorModel <- as.character(object["model"]["error.model"])
  residualFree <- copulaScoreResidualIndices(errorModel)

  data <- object["data"]["data"]
  subject <- match(data[, object["data"]["name.group"]],
    unique(data[, object["data"]["name.group"]]))
  observations <- as.numeric(data[, object["data"]["name.response"]])
  predictors <- data[, object["data"]["name.predictors"], drop = FALSE]
  ytype <- if (nchar(object["data"]["name.ytype"]))
    data[, object["data"]["name.ytype"]] else rep(1L, nrow(data))
  structural <- object["model"]["model"]
  n <- nrow(meanPhi)
  conditioning <- if ((state$dConditioning %||% 0L) > 0L)
    as.matrix(state$conditioning) else matrix(numeric(), n, 0L)
  natural <- identical(state$populationScale, "parameter")

  parameters <- copulaFisherParameterVector(state, locations, residual,
    residualFree)
  theta <- parameters$value

  ## Per-subject complete-data log density at one draw of the random effects.
  perSubject <- function(value, eta) {
    cursor <- 0L
    location <- if (parameters$nLocation) {
      cursor <- parameters$nLocation
      value[seq_len(cursor)]
    } else locations
    marginParameters <- if (parameters$nMargin) {
      got <- value[cursor + seq_len(parameters$nMargin)]
      cursor <- cursor + parameters$nMargin
      got
    } else numeric()
    angles <- if (parameters$nAngle) {
      got <- value[cursor + seq_len(parameters$nAngle)]
      cursor <- cursor + parameters$nAngle
      got
    } else numeric()
    thisResidual <- residual
    if (length(residualFree))
      thisResidual[residualFree] <- value[cursor + seq_along(residualFree)]
    margins <- state$margins
    margins[seq_len(dEta)] <- copulaMarginsWithParameters(
      state$margins[seq_len(dEta)], parameters$marginLayout, marginParameters)
    correlation <- if (parameters$nAngle)
      copulaGaussianCorrelationFromAngles(angles, state$d)$R else
      copulaGaussianRvineCor(state$vine, state$d)
    vine <- copulaGaussianRvineUpdateCor(state$vine, correlation)
    predictor <- matrix(location, n, dEta, byrow = TRUE)
    population <- if (natural)
      copulaNaturalFremLogPrior(eta, conditioning, vine, margins, dEta,
        predictor, transform, "joint") else
      copulaGaussianFremLogPrior(cbind(eta, conditioning), vine, margins,
        dEta, "joint")
    phi <- predictor + eta
    psi <- copulaWorkingToNatural(phi, transform)
    prediction <- structural(psi, subject, predictors)
    spread <- error(prediction, thisResidual, ytype)
    if (any(!is.finite(spread)) || any(spread <= 0))
      return(rep(-Inf, n))
    response <- -.5 * ((observations - prediction) / spread)^2 - log(spread) -
      .5 * log(2 * pi)
    population + as.numeric(rowsum(response, subject, reorder = FALSE))
  }

  list(perSubject = perSubject, theta = theta, parameters = parameters,
    subjects = n, state = state)
}

## How many draws this needs is not a detail. For a population location the
## complete-data score is large -- given the random effect the data pin the
## individual parameter down -- while the observed-data score is small, because
## the random effect absorbs a shift in the location. Fisher's identity turns
## the first into the second by averaging over the posterior, so the entire
## cancellation lives in that average, and an under-dispersed sample cancels
## too little: the score stays too large, the information comes out too big and
## the standard error too small. Measured on a fit whose location could not be
## known better than 0.0138, the reported standard error ran
##
##   200 draws 0.0075, 400 0.0101, 800 0.0112, 1600 0.0126,
##
## against an empirical spread of 0.0142, while the margin and correlation
## coordinates moved in the fourth decimal. So the default is generous, and the
## returned `stability` says whether it was generous enough: it compares the
## standard errors from the full sample against those from half of it, and a
## value much above one means the average has not settled.
copulaScoreFisherInformation <- function(object, draws = 1000L,
    max.iter = NULL, seed = 935001L, step = 1e-4, check = TRUE, thin = 1L) {
  ## Whether the posterior average has settled is a question about the length
  ## of the chain behind the draws, not about how many of them were kept, so a
  ## split of one sample cannot answer it -- both halves inherit the same
  ## under-dispersion. It is answered by running the whole thing again on a
  ## shorter chain and seeing whether the answer moved.
  if (isTRUE(check)) {
    coarse <- copulaScoreFisherInformation(object,
      draws = max(4L, as.integer(draws / 2L)),
      max.iter = if (is.null(max.iter)) NULL else
        max(1L, as.integer(max.iter / 2L)),
      seed = seed, step = step, check = FALSE, thin = thin)
    full <- copulaScoreFisherInformation(object, draws = draws,
      max.iter = max.iter, seed = seed, step = step, check = FALSE,
      thin = thin)
    stability <- if (is.null(full$standardErrors) ||
        is.null(coarse$standardErrors)) NULL else
      full$standardErrors / coarse$standardErrors
    if (!is.null(stability) && max(abs(stability - 1), na.rm = TRUE) > .1)
      warning("score Fisher information has not settled: halving the chain ",
        "moves a standard error by ",
        sprintf("%.0f%%", 100 * max(abs(stability - 1), na.rm = TRUE)),
        "; increase draws and max.iter", call. = FALSE)
    full$stability <- stability
    return(full)
  }
  context <- copulaFisherContext(object)
  perSubject <- context$perSubject
  theta <- context$theta
  n <- context$subjects

  pool <- copulaPosteriorEtaDraws(object, draws, max.iter, as.integer(seed),
    thin = thin)
  samples <- dim(pool)[3L]
  if (samples < 4L)
    stop("score Fisher information needs at least four posterior draws")
  ## The draws are split in two independent halves, because the naive estimator
  ## is biased at a finite number of them. Delta_i is estimated, not known, so
  ##
  ##   E[ Delta_i^ Delta_i^' ] = Delta_i Delta_i' + Var(Delta_i^),
  ##
  ## and squaring an estimate inflates the information by its own Monte Carlo
  ## variance, which understates every standard error that comes out of the
  ## inverse. Halves A and B are independent given the data, so the product of
  ## one with the other has the cross term and not the variance. The naive form
  ## is returned alongside it: the gap between them is what the bias was.
  half <- rep(c(1L, 2L), length.out = samples)
  scores <- matrix(0, n, length(theta))
  halves <- list(matrix(0, n, length(theta)), matrix(0, n, length(theta)))
  squares <- matrix(0, n, length(theta))
  for (r in seq_len(samples)) {
    eta <- pool[, , r, drop = FALSE]
    dim(eta) <- dim(pool)[1:2]
    for (k in seq_along(theta)) {
      h <- step * max(1, abs(theta[k]))
      up <- down <- theta; up[k] <- up[k] + h; down[k] <- down[k] - h
      difference <- (perSubject(up, eta) - perSubject(down, eta)) / (2 * h)
      if (any(!is.finite(difference)))
        stop("non-finite per-subject score at coordinate ", k)
      scores[, k] <- scores[, k] + difference / samples
      squares[, k] <- squares[, k] + difference^2 / samples
      halves[[half[r]]][, k] <- halves[[half[r]]][, k] +
        difference / sum(half == half[r])
    }
  }
  ## Monte Carlo standard error of each per-subject score, from the spread of
  ## the draws that produced it.
  ## pmax() takes its attributes from the first argument, so the matrix has to
  ## come first or the dimensions are silently dropped.
  standardError <- sqrt(pmax(squares - scores^2, 0) / samples)

  naive <- crossprod(scores) / n
  cross <- crossprod(halves[[1L]], halves[[2L]]) / n
  information <- (cross + t(cross)) / 2
  covarianceOf <- function(matrix) {
    decomposition <- try(chol(matrix), silent = TRUE)
    if (inherits(decomposition, "try-error")) NULL else
      chol2inv(decomposition) / n
  }
  covariance <- covarianceOf(information)
  dimnames(information) <- dimnames(naive) <- list(names(theta), names(theta))
  list(information = information, naiveInformation = naive,
    covariance = covariance,
    standardErrors = if (is.null(covariance)) NULL else
      stats::setNames(sqrt(pmax(0, diag(covariance))), names(theta)),
    naiveStandardErrors = local({
      got <- covarianceOf(naive)
      if (is.null(got)) NULL else
        stats::setNames(sqrt(pmax(0, diag(got))), names(theta))
    }),
    scores = scores, scoreStandardError = standardError, parameters = theta,
    meanScore = colMeans(scores), draws = samples,
    method = "delattre-kuhn-score-covariance, split-half")
}
