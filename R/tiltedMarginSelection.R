## Marginal selection through a single nested family.
##
## The screening in naturalMarginSelection.R compares named parametric
## families: it needs a candidate list, it fits each candidate by a
## quasi-Newton search over the bridge objective, and the answer can only ever
## be one of the shapes on the list. This selector replaces the list with one
## ladder. For each eta coordinate the marginal is an exponentially tilted base
## family of degree K = 0, 1, ..., which is lognormal (or Normal) at K = 0 and
## nested in K, and the fit at each degree is the maximum-likelihood fit to the
## pooled conditional distribution of that coordinate -- cheap and reliable
## because the family is a regular exponential family with an analytic score.
##
## The construction has three parts, and each of them is exact:
##
## 1. Fisher's identity makes the pooled posterior the right target. The
##    derivative of the marginal log-likelihood in the margin parameters is
##    sum_i E_{pi_0i}[ d log f(psi_ij) ], so a stationary point of the marginal
##    likelihood in the margin alone is a maximiser of the pooled posterior
##    log-likelihood. Fitting the margin to the posterior draws is therefore an
##    exact step of the marginal problem, not an approximation of it, and it is
##    the sense in which "choose the marginal from the individual conditional
##    distributions" is a theorem rather than a heuristic.
##
## 2. Within a degree the fit has an analytic gradient and, at a fixed
##    location, an exactly concave objective whose stationarity condition reads
##      (1/N) sum_i E_{pi_0i}[ L_k(z_ij) ]  =  E_beta[ L_k ]
##    -- the fitted marginal reproduces the pooled posterior moments of the
##    basis. The anchor ties the location to the coefficients, which adds a
##    rank-one term to the gradient and leaves the objective concave only near
##    its maximum; the fit starts from beta = 0 every time, so there is still no
##    starting value to choose and no restart.
##
## 3. Across degrees the comparison is the exact posterior likelihood-ratio
##    identity already used by the family screen, penalised by BIC. The ladder
##    is nested, so the penalty is the only thing that can stop it, and the
##    degree it stops at is a genuine model choice rather than a shortlist
##    membership.
##
## The cost is d * (maxDegree + 1) exact evaluations per sweep against K^d
## combinations for the product grid, and each evaluation fits one margin from
## a fixed starting point rather than searching a family it was handed.

## One rung of one coordinate's ladder.
##
## The margin fit depends only on the coordinate, the degree and the pooled
## draws -- never on what the other coordinates are doing -- so it is cached
## across sweeps rather than repeated. Only the validation metric decides
## anything; the training one is a diagnostic, and computing it for every trial
## doubled the bridge work, so it is left to the caller to ask for once the
## ladder has settled.
copulaTiltedLadderFit <- function(degree, coordinate, families, margins,
                                  supports, train, validation, incumbent,
                                  fits = NULL, withTraining = FALSE) {
  slot <- paste(coordinate, degree, sep = "|")
  fit <- if (!is.null(fits)) fits[[slot]] else NULL
  if (is.null(fit)) {
    fit <- copulaNaturalMarginTiltedFit(train$natural[, coordinate],
      train$typical[, coordinate], degree, supports[coordinate])
    if (!is.null(fits)) assign(slot, fit, envir = fits)
  }
  margins[[coordinate]] <- fit$margin
  training <- if (withTraining)
    copulaNaturalBridgeMetric(margins, train, incumbent) else NULL
  validationMetric <- copulaNaturalBridgeMetric(margins, validation, incumbent)
  incumbentParameters <- sum(vapply(incumbent$margins[seq_len(incumbent$dEta)],
    function(margin) sum(margin$free), integer(1)))
  candidateParameters <- sum(vapply(margins,
    function(margin) sum(margin$free), integer(1)))
  deltaParameters <- candidateParameters - incumbentParameters
  families[coordinate] <- paste0("tilted", degree)
  list(families = families, margins = margins, degree = degree,
    coordinate = coordinate, convergence = 0L,
    newtonIterations = fit$newtonIterations,
    training = training, validation = validationMetric,
    deltaParameters = deltaParameters,
    validationBicAdvantage = 2 * validationMetric$deltaLogLik -
      deltaParameters * log(train$n))
}

copulaSelectTiltedMargins <- function(object, supports, maxDegree = 4L,
    posteriorDraws = 500L, max.iter = NULL, seed = 935001L, maxSweeps = 3L,
    minimumEssFraction = .01, minimumPosteriorEssFraction = .05,
    maximumMcse = .5) {
  if (!inherits(object, "SaemixObject"))
    stop("object must be a fitted SaemixObject")
  state <- copulaGet(object)
  if (!state$populationScale %in% c("transformed-additive", "parameter") ||
      !copulaIsFullGaussianVine(state$vine, state$d))
    stop("tilted-margin selection requires a Gaussian-copula incumbent")
  supports <- rep_len(as.character(supports), state$dEta)
  if (any(!supports %in% c("real", "positive")))
    stop("tilted margins support real or positive coordinates")
  maxDegree <- as.integer(maxDegree)
  if (maxDegree < 0L || maxDegree > 6L)
    stop("maxDegree must be between 0 and 6")
  trainingDraws <- copulaPosteriorEtaDraws(object, posteriorDraws, max.iter,
    as.integer(seed))
  validationDraws <- copulaPosteriorEtaDraws(object, posteriorDraws, max.iter,
    as.integer(seed) + 1L)
  train <- copulaNaturalBridgePrepare(
    copulaNaturalPosteriorData(object, trainingDraws), state)
  validation <- copulaNaturalBridgePrepare(
    copulaNaturalPosteriorData(object, validationDraws), state)

  ## Start every coordinate at the base of its own ladder, then raise one
  ## coordinate at a time. Scored trials are cached by the full family vector,
  ## so a coordinate that is revisited without having moved costs nothing, and
  ## the margin fits behind them are cached by coordinate and degree, since a
  ## fit never depends on what the other coordinates are doing.
  store <- new.env(parent = emptyenv())
  fits <- new.env(parent = emptyenv())
  key <- function(families) paste(families, collapse = "|")
  current <- rep("tilted0", state$dEta)
  margins <- lapply(seq_len(state$dEta), function(j) {
    base <- copulaNaturalMarginTiltedFit(train$natural[, j],
      train$typical[, j], 0L, supports[j])
    assign(paste(j, 0L, sep = "|"), base, envir = fits)
    base$margin
  })
  for (sweep in seq_len(max(1L, as.integer(maxSweeps)))) {
    moved <- FALSE
    for (j in seq_len(state$dEta)) {
      trials <- lapply(0:maxDegree, function(degree) {
        trial <- current; trial[j] <- paste0("tilted", degree); trial })
      pending <- which(vapply(trials,
        function(t) is.null(store[[key(t)]]), logical(1)))
      if (length(pending)) {
        got <- copulaMapCandidates(length(pending), function(k)
          copulaTiltedLadderFit(pending[k] - 1L, j, current, margins, supports,
            train, validation, state, fits = fits))
        ## copulaMapCandidates may have run these in separate processes, in
        ## which case the fits they made did not come back with them.
        for (k in seq_along(pending)) {
          slot <- paste(j, pending[k] - 1L, sep = "|")
          if (is.null(fits[[slot]]))
            assign(slot, list(margin = got[[k]]$margins[[j]],
              newtonIterations = got[[k]]$newtonIterations), envir = fits)
        }
        for (k in seq_along(pending))
          assign(key(trials[[pending[k]]]), got[[k]], envir = store)
      }
      scores <- vapply(trials, function(t) {
        got <- store[[key(t)]]
        if (is.null(got) || !is.finite(got$validationBicAdvantage) ||
            got$validation$minimumEssFraction < minimumEssFraction ||
            got$validation$mcse > maximumMcse) -Inf else
          got$validationBicAdvantage
      }, numeric(1))
      if (all(!is.finite(scores)))
        stop("no tilt degree passed the overlap diagnostics")
      best <- trials[[which.max(scores)]]
      if (!identical(best[j], current[j])) {
        current <- best; margins <- store[[key(best)]]$margins; moved <- TRUE
      }
    }
    if (!moved) break
  }
  fitted <- unname(as.list(store))
  table <- do.call(rbind, lapply(seq_along(fitted), function(i) {
    result <- fitted[[i]]
    data.frame(model = i, families = paste(result$families, collapse = "/"),
      coordinate = result$coordinate, degree = result$degree,
      newton_iterations = result$newtonIterations,
      delta_parameters = result$deltaParameters,
      training_delta_loglik = if (is.null(result$training)) NA_real_ else
        result$training$deltaLogLik,
      validation_delta_loglik = result$validation$deltaLogLik,
      validation_bic_advantage = result$validationBicAdvantage,
      validation_mcse = result$validation$mcse,
      minimum_ess_fraction = result$validation$minimumEssFraction,
      maximum_pareto_k = result$validation$maximumParetoK %||% NA_real_)
  }))
  selected <- which(table$families == key(current))[1L]
  ## The training metric is a diagnostic -- it says whether the pool the margin
  ## was fitted on and the pool it was scored on agree -- so it is computed for
  ## the retained ladder only, not for every rung of every coordinate.
  trainingMetric <- copulaNaturalBridgeMetric(margins, train, state)
  table$training_delta_loglik[selected] <- trainingMetric$deltaLogLik
  trainDiag <- attr(trainingDraws, "diagnostics")
  validationDiag <- attr(validationDraws, "diagnostics")
  mixing <- trainDiag$minimumMcmcEssFraction >= minimumPosteriorEssFraction &&
    validationDiag$minimumMcmcEssFraction >= minimumPosteriorEssFraction
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
  answer <- list(selected = selected, table = table, families = current,
    degrees = as.integer(sub("tilted", "", current)), margins = margins,
    population = population, posteriorDrawsPerPool = posteriorDraws,
    diagnostics = list(posteriorMixingAdequate = mixing,
      fixedDuringFinalFit = TRUE, selectionOutsideScoreSa = TRUE))
  class(answer) <- c("saemixTiltedMarginSelection", "list")
  answer
}

print.saemixTiltedMarginSelection <- function(x, ...) {
  cat("Tilted-margin degree selection\n")
  cat("  selected degrees:", paste(x$degrees, collapse = ", "), "\n")
  cat("  posterior draws per pool:", x$posteriorDrawsPerPool, "\n\n")
  print(x$table[, c("families", "coordinate", "degree", "delta_parameters",
    "validation_delta_loglik", "validation_bic_advantage",
    "minimum_ess_fraction")], row.names = FALSE)
  invisible(x)
}
