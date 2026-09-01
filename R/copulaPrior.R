## Runtime state for the fixed Gaussian-copula score model.

.cop <- new.env(parent = emptyenv())
.trc <- new.env(parent = emptyenv())

copulaClear <- function() rm(list = ls(.cop), envir = .cop)

copulaPopulation <- function(vine, margins = NULL, sd = NULL,
                             scale = c("auto", "transformed-additive", "parameter"),
                             conditioning = NULL,
                             populationAlgorithm = "score-sa",
                             scoreScale = "auto", scoreFiniteDifference = 1e-4,
                             scoreProjection = 24, scoreGainScale = 0.2,
                             scoreGainPower = 0.8, scoreGainOffset = 30,
                             scoreBurn = 50L, ...) {
  scale <- match.arg(scale)
  structure(list(arguments = c(list(vine = vine, margins = margins, sd = sd,
    populationScale = scale, conditioning = conditioning,
    populationAlgorithm = populationAlgorithm, scoreScale = scoreScale,
    scoreFiniteDifference = scoreFiniteDifference,
    scoreProjection = scoreProjection, scoreGainScale = scoreGainScale,
    scoreGainPower = scoreGainPower, scoreGainOffset = scoreGainOffset,
    scoreBurn = scoreBurn), list(...))), class = "saemixPopulation")
}

copulaUsePopulation <- function(population) {
  if (!inherits(population, "saemixPopulation") ||
      !is.list(population$arguments))
    stop("population must be created by copulaPopulation()")
  do.call(copulaSet, population$arguments)
  invisible(TRUE)
}

copulaRestoreState <- function(state) {
  copulaClear()
  for (name in names(state)) assign(name, state[[name]], envir = .cop)
  invisible(TRUE)
}

copulaVineForMargins <- function(vine, margins) {
  types <- ifelse(vapply(margins, `[[`, character(1), "type") == "discrete",
    "d", "c")
  if (identical(as.character(vine$var_types), types)) return(vine)
  suppressWarnings(rvinecopulib::vinecop_dist(
    vine$pair_copulas, vine$structure, var_types = types))
}

copulaSet <- function(vine, margins = NULL, sd = NULL,
                      conditioning = NULL,
                      populationScale = c("auto", "transformed-additive", "parameter"),
                      populationAlgorithm = "score-sa",
                      likelihoodTarget = "joint", familySet = NULL,
                      selectStructure = FALSE, refitEvery = 1L,
                      activeFrom = 1L, fitFrom = 1L,
                      warmStartOnActivate = FALSE, guard = FALSE,
                      numericalPolicy = "exact", scoreScale = "auto",
                      scoreFiniteDifference = 1e-4, scoreProjection = 24,
                      scoreGainScale = 0.2, scoreGainPower = 0.8,
                      scoreGainOffset = 30, scoreBurn = 50L, ...) {
  previous <- as.list(.cop)
  committed <- FALSE
  on.exit(if (!committed) copulaRestoreState(previous), add = TRUE)
  populationScale <- match.arg(populationScale)

  dots <- list(...)
  if (length(dots))
    stop("unsupported population controls: ", paste(names(dots), collapse = ", "))
  if (!identical(populationAlgorithm, "score-sa"))
    stop("only the likelihood-score population estimator is supported")
  if (!identical(likelihoodTarget, "joint") ||
      !is.null(familySet) || isTRUE(selectStructure) ||
      as.integer(refitEvery) != 1L || as.integer(activeFrom) != 1L ||
      as.integer(fitFrom) != 1L)
    stop("score-sa requires one fixed joint model from the first iteration")
  if (!identical(numericalPolicy, "exact"))
    stop("score-sa requires exact density evaluation")
  if (!inherits(vine, "vinecop_dist"))
    stop("vine must inherit from vinecop_dist")

  d <- as.integer(vine$structure$d)
  if (is.null(conditioning)) {
    conditioningValues <- NULL
    conditioningNames <- character()
    dConditioning <- 0L
  } else {
    if (is.list(conditioning) && !is.data.frame(conditioning)) {
      conditioningNames <- conditioning$variableName %||% character()
      conditioningValues <- conditioning$values
    } else {
      conditioningNames <- character()
      conditioningValues <- conditioning
    }
    conditioningValues <- as.matrix(conditioningValues)
    storage.mode(conditioningValues) <- "double"
    if (!nrow(conditioningValues) || !ncol(conditioningValues) ||
        any(is.infinite(conditioningValues)))
      stop("conditioning must be a non-empty numeric matrix; NA is allowed")
    dConditioning <- ncol(conditioningValues)
    if (dConditioning >= d || any(colSums(!is.na(conditioningValues)) == 0L))
      stop("each conditioning coordinate needs observed values and one latent coordinate")
    if (!length(conditioningNames)) conditioningNames <- colnames(conditioningValues)
    if (is.null(conditioningNames) || length(conditioningNames) != dConditioning)
      conditioningNames <- paste0("conditioning", seq_len(dConditioning))
  }
  dEta <- d - dConditioning

  naturalMargins <- !is.null(margins) && all(vapply(margins,
    inherits, logical(1), "saemix_natural_parameter_margin"))
  etaMargins <- !is.null(margins) && all(vapply(margins,
    inherits, logical(1), "saemix_copula_margin"))
  if (identical(populationScale, "auto")) populationScale <-
    if (naturalMargins) "parameter" else "transformed-additive"
  if (is.null(margins)) {
    if (is.null(sd) || length(sd) != d || any(!is.finite(sd)) || any(sd <= 0))
      stop("supply one margin per coordinate or positive Gaussian scales")
    margins <- lapply(sd, copulaMarginNormal)
  } else if (identical(populationScale, "parameter")) {
    if (!naturalMargins || !is.null(sd) || length(margins) != d)
      stop("parameter-scale fitting requires one natural margin per coordinate")
    if (dConditioning > 0L)
      stop("parameter-scale conditioning margins are not yet enabled")
    lapply(margins, copulaNaturalMarginValidate)
  } else {
    if (!is.null(sd) || !is.list(margins) || length(margins) != d)
      stop("supply exactly one margin per vine coordinate")
    if (!etaMargins)
      stop("transformed-additive fitting requires eta-margin objects")
    lapply(margins, copulaMarginValidate)
  }
  if (identical(populationScale, "transformed-additive")) {
    etaMargins <- margins[seq_len(dEta)]
    if (any(!vapply(etaMargins, `[[`, logical(1), "centered")) ||
        any(vapply(etaMargins, function(m)
          any(m$free & m$roles == "location"), logical(1))) ||
        any(vapply(etaMargins, `[[`, character(1), "type") != "continuous") ||
        any(!vapply(etaMargins, `[[`, logical(1), "scale_is_sd")))
      stop("eta margins must be centered continuous laws with finite SD scales")
  }

  vine <- copulaVineForMargins(vine, margins)
  if (!copulaIsFullGaussianVine(vine, d))
    stop("score-sa requires a fixed full Gaussian R-vine")
  scoreScaleRequested <- scoreScale
  scoreScaleAutomatic <- is.character(scoreScale) &&
    length(scoreScale) == 1L && identical(tolower(scoreScale), "auto")
  if (is.character(scoreScale) && !scoreScaleAutomatic)
    stop("scoreScale must be 'auto' or one positive numeric value")
  if (scoreScaleAutomatic) scoreScale <- 1
  numericScalar <- function(x, positive = FALSE) length(x) == 1L &&
    is.finite(x) && (!positive || x > 0)
  if (!numericScalar(scoreScale, TRUE) ||
      !numericScalar(scoreFiniteDifference, TRUE) ||
      !numericScalar(scoreProjection, TRUE) ||
      !numericScalar(scoreGainScale, TRUE) ||
      !numericScalar(scoreGainPower, TRUE) || scoreGainPower <= .75 ||
      scoreGainPower > 1 || !numericScalar(scoreGainOffset) ||
      scoreGainOffset < 0 || length(scoreBurn) != 1L ||
      is.na(scoreBurn) || scoreBurn < 0)
    stop("invalid score stochastic-approximation controls")

  etaScale <- if (identical(populationScale, "parameter"))
    rep(.3, dEta) else copulaMarginScales(margins)[seq_len(dEta)]
  etaCorrelation <- copulaGaussianRvineCor(vine, d)[seq_len(dEta),
    seq_len(dEta), drop = FALSE]
  proposalOmega <- diag(etaScale, dEta) %*% etaCorrelation %*%
    diag(etaScale, dEta)

  copulaClear()
  values <- list(vine = vine, margins = margins,
    populationScale = populationScale, populationAlgorithm = "score-sa",
    likelihoodTarget = "joint", numericalPolicy = "exact",
    d = d, dEta = dEta, dConditioning = dConditioning,
    conditioning = conditioningValues,
    conditioningName = as.character(conditioningNames),
    sd = if (identical(populationScale, "parameter")) rep(NA_real_, d) else
      copulaMarginScales(margins), proposalOmega = proposalOmega,
    scoreScale = as.numeric(scoreScale),
    scoreScaleRequested = scoreScaleRequested,
    scoreScaleSource = if (scoreScaleAutomatic)
      "unit-after-finite-diagonal-score-metric" else "numeric-override",
    scoreFiniteDifference = as.numeric(scoreFiniteDifference),
    scoreProjection = as.numeric(scoreProjection),
    scoreGainScale = as.numeric(scoreGainScale),
    scoreGainPower = as.numeric(scoreGainPower),
    scoreGainOffset = as.numeric(scoreGainOffset), scoreBurn = as.integer(scoreBurn),
    scoreState = NULL, mode = "joint", modelFrozen = TRUE,
    activeFrom = 1L, activated = TRUE, augmentMissingGaussian = TRUE,
    warmStartOnActivate = isTRUE(warmStartOnActivate), guard = isTRUE(guard),
    estimatedMargins = TRUE, estimatedVine = TRUE,
    locationRankChecked = FALSE, trace = list(), lastJoint = NULL,
    timing = list(estep = 0, batchUpdate = 0, iteration = 0))
  for (name in names(values)) assign(name, values[[name]], envir = .cop)
  .cop$fingerprint <- copulaFingerprint(vine, d, margins)
  committed <- TRUE
  invisible(TRUE)
}

copulaActive <- function() !is.null(.cop$vine)
copulaActiveAt <- function(kiter) copulaActive()
copulaStepSize <- function(kiter, stockGamma) stockGamma
copulaActivateFromGaussian <- function(kiter, omega) invisible(FALSE)

copulaRecordPhase <- function(kiter, gamma) {
  .cop$lastGamma <- gamma
  .cop$phase <- if (gamma >= 1 - sqrt(.Machine$double.eps))
    "constant" else "decreasing"
  invisible(.cop$phase)
}

copulaProfileEvent <- function(...) invisible(NULL)

withSeed <- function(seed, expression) {
  hadSeed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (hadSeed) previous <- get(".Random.seed", envir = globalenv())
  on.exit(if (hadSeed) assign(".Random.seed", previous, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv()), add = TRUE)
  set.seed(seed)
  force(expression)
}

copulaTakeBeta <- function() {
  value <- .cop$betaJoint
  .cop$betaJoint <- NULL
  value
}

copulaTakeResidual <- function() {
  value <- .cop$residualJoint
  .cop$residualJoint <- NULL
  value
}

copulaAcceptedBeta <- function(current, fitted = NULL, locMap = NULL,
                               betaFree = NULL) {
  current <- as.numeric(current)
  if (is.null(fitted) || length(fitted) != length(current) ||
      any(!is.finite(fitted)) || is.null(locMap)) return(current)
  active <- which(rowSums(abs(as.matrix(locMap))) > 0)
  if (!is.null(betaFree)) active <- intersect(active, as.integer(betaFree))
  current[active] <- fitted[active]
  current
}

copulaGet <- function(object = NULL) {
  if (is.null(object)) return(as.list(.cop))
  state <- attr(object, "saemix.copula", exact = TRUE)
  if (is.null(state)) stop("object does not contain a fitted copula state")
  state
}

copulaSnapshot <- function(state = copulaGet(), etaIndex = NULL,
                           variableName = NULL) {
  keep <- c("vine", "margins", "sd", "d", "dEta", "dConditioning",
    "conditioning", "conditioningName", "fingerprint", "modelFrozen",
    "lastJoint", "npar.margin", "npar.vine", "estimatedMargins",
    "estimatedVine", "timing", "trace", "populationScale", "proposalOmega",
    "likelihoodTarget", "numericalPolicy", "populationAlgorithm",
    "scoreScale", "scoreScaleRequested", "scoreScaleSource",
    "scoreFiniteDifference", "scoreProjection",
    "scoreGainScale", "scoreGainPower", "scoreGainOffset", "scoreBurn",
    "scoreState", "rwProposalFrozen", "rwProposalFreezeIteration",
    "rwBlockSizeFrozen")
  out <- state[intersect(keep, names(state))]
  out$version <- 2L
  out$etaIndex <- if (is.null(etaIndex)) seq_len(state$dEta) else
    as.integer(etaIndex)
  out$variableName <- variableName
  out$fingerprint <- copulaFingerprint(state$vine, state$d, state$margins)
  class(out) <- c("saemixCopulaSnapshot", "list")
  out
}

copulaPostprocessGuard <- function(object, feature) {
  if (is.null(attr(object, "saemix.copula", exact = TRUE)))
    return(invisible(TRUE))
  stop(feature, " assumes Gaussian random effects; use the copula-aware method",
    call. = FALSE)
}

copulaPadFlat <- function(vine, d) {
  edges <- unlist(vine$pair_copulas, recursive = FALSE)
  full <- d * (d - 1L) / 2L
  if (length(edges) < full)
    edges <- c(edges, rep(list(rvinecopulib::bicop_dist("indep")),
      full - length(edges)))
  edges
}

copulaFingerprint <- function(vine, d, margins = NULL) {
  edges <- copulaPadFlat(vine, d)
  answer <- list(structure = vine$structure,
    truncLvl = as.integer(vine$structure$trunc_lvl),
    varTypes = as.character(vine$var_types),
    family = vapply(edges, `[[`, character(1), "family"),
    rotation = vapply(edges, `[[`, numeric(1), "rotation"),
    npar = vapply(edges, function(edge) length(edge$parameters), integer(1)))
  if (!is.null(margins)) {
    answer$marginName <- unname(vapply(margins, `[[`, character(1), "name"))
    answer$marginType <- unname(vapply(margins, `[[`, character(1), "type"))
    answer$marginFree <- unname(lapply(margins, `[[`, "free"))
    answer$marginRoles <- unname(lapply(margins, `[[`, "roles"))
    answer$marginClass <- unname(vapply(margins,
      function(margin) class(margin)[1L], character(1)))
  }
  answer
}

copulaAssertFrozen <- function(vine = .cop$vine) {
  current <- copulaFingerprint(vine, .cop$d, .cop$margins)
  if (!identical(current, .cop$fingerprint))
    stop("fixed population model changed during estimation", call. = FALSE)
  invisible(TRUE)
}

copulaLocation <- function(X, beta, locMap) sweep(X, 2L, beta, "*") %*% locMap

copulaUeta <- function(eta) {
  eta <- as.matrix(eta)
  if (ncol(eta) != .cop$d)
    stop("population state and sampled coordinates have different dimensions")
  -copulaGaussianFremLogPrior(eta, .cop$vine, .cop$margins, .cop$dEta,
    "joint")
}

copulaRandEta <- function(n, conditioning = NULL) {
  if (.cop$dConditioning > 0L) {
    conditioning <- as.matrix(conditioning)
    if (nrow(conditioning) != n || ncol(conditioning) != .cop$dConditioning)
      stop("conditioning rows do not match the requested population draws")
    return(copulaGaussianFremRandEta(
      conditioning, .cop$vine, .cop$margins, .cop$dEta))
  }
  uniforms <- rvinecopulib::rvinecop(n, .cop$vine)
  copulaMarginsQuantile(uniforms, .cop$margins)
}

populationRandPhi <- function(object, predictor) {
  state <- copulaGet(object)
  predictor <- as.matrix(predictor)
  if (state$dConditioning > 0L)
    stop("joint FREM simulation requires explicit conditioning covariates")
  uniforms <- rvinecopulib::rvinecop(nrow(predictor), state$vine)
  predictor + copulaMarginsQuantile(uniforms, state$margins)
}

copulaOmega <- function() {
  correlation <- copulaGaussianRvineCor(.cop$vine, .cop$d)
  scale <- if (identical(.cop$populationScale, "parameter"))
    as.numeric(.cop$sd[seq_len(.cop$dEta)]) else
    copulaMarginScales(.cop$margins)[seq_len(.cop$dEta)]
  if (any(!is.finite(scale)) || any(scale <= 0))
    scale <- sqrt(diag(.cop$proposalOmega))[seq_len(.cop$dEta)]
  correlation <- correlation[seq_len(.cop$dEta), seq_len(.cop$dEta),
    drop = FALSE]
  diag(scale, .cop$dEta) %*% correlation %*% diag(scale, .cop$dEta)
}

copulaIsNestedGaussian <- function(vine, margins) {
  copulaIsFullGaussianVine(vine, length(margins)) &&
    all(vapply(margins, function(margin) identical(margin$name, "normal"),
      logical(1)))
}

copulaNestedGaussianLogDensity <- function(eta, vine, margins) {
  if (!copulaIsNestedGaussian(vine, margins)) return(NULL)
  scale <- vapply(margins, function(margin)
    margin$parameters[["sd"]], numeric(1))
  correlation <- copulaGaussianRvineCor(vine, length(margins))
  covariance <- diag(scale) %*% correlation %*% diag(scale)
  copulaGaussianLogDensity(as.matrix(eta), covariance)
}

copulaLogPrior <- function(E, vine, sdv = NULL, cores = 1L, margins = NULL,
                           numericalPolicy = "exact", ...) {
  if (!identical(numericalPolicy, "exact"))
    stop("only exact likelihood evaluation is supported")
  if (is.null(margins)) margins <- lapply(sdv, copulaMarginNormal)
  copulaGaussianFremLogPrior(as.matrix(E), vine, margins,
    ncol(as.matrix(E)), "joint")
}

saemixTraceReset <- function() {
  .trc$L <- list()
  options(saemixTrace = TRUE)
}
saemixTraceGet <- function() .trc$L
.saemixTracePush <- function(kiter, betas, omdiag, pres, gamma, statrese) {
  .trc$L[[length(.trc$L) + 1L]] <- list(kiter = kiter, betas = betas,
    omdiag = omdiag, pres = pres, gamma = gamma, statrese = statrese)
}
