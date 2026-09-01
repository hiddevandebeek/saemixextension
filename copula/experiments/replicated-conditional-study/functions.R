study_truth <- function() list(V = 10, CL = 3, etaV = .22, etaCL = .35,
  residual = .12, gammaShape = 2, gammaScale = 20,
  albuminMean = 40, albuminSd = 5)

study_correlation <- function() matrix(c(
  1, .20, .10, .25,
  .20, 1, .55, .10,
  .10, .55, 1, .20,
  .25, .10, .20, 1), 4L, 4L, byrow = TRUE)

study_structure <- function() rvinecopulib::cvine_structure(c(3, 1, 4, 2))

study_pk <- function(psi, id, xidep) {
  V <- psi[id, 1L]
  CL <- psi[id, 2L]
  xidep[, 1L] / V * exp(-(CL / V) * xidep[, 2L])
}

study_model <- function(fixed = c(V = 9, CL = 2.7),
                        omega = diag(c(.22^2, .35^2)), residual = .12) {
  errorInitial <- if (length(residual) == 1L) c(0, residual) else residual
  saemixModel(model = study_pk, modeltype = "structural",
    description = "Replicated Gamma-tail conditional prediction study",
    psi0 = matrix(fixed, 1L, dimnames = list(NULL, c("V", "CL"))),
    transform.par = c(1, 1), covariance.model = matrix(1, 2L, 2L),
    omega.init = omega, error.model = "proportional",
    error.init = errorInitial, verbose = FALSE)
}

study_residual <- function(fit) {
  value <- as.numeric(fit@results@respar[fit@results@indx.res])
  if (length(value) != 1L || !is.finite(value) || value <= 0)
    stop("could not identify the fitted proportional residual SD")
  value
}

study_conditional_draw <- function(covariate, vine, margins, zRandom,
                                   dEta = 2L, covariateIndex = 3L) {
  Rfull <- copulaGaussianRvineCor(vine, length(margins))
  index <- c(seq_len(dEta), covariateIndex)
  R <- Rfull[index, index, drop = FALSE]
  u <- margins[[covariateIndex]]$cdf(covariate,
    margins[[covariateIndex]]$parameters)
  zc <- qnorm(pmin(pmax(u, 1e-12), 1 - 1e-12))
  D <- diag(copulaMarginScales(margins)[seq_len(dEta)], dEta)
  regression <- D %*% R[seq_len(dEta), dEta + 1L, drop = FALSE]
  covariance <- D %*% (R[seq_len(dEta), seq_len(dEta), drop = FALSE] -
    tcrossprod(R[seq_len(dEta), dEta + 1L])) %*% D
  outer(zc, as.numeric(regression)) + zRandom %*% chol(covariance)
}

study_response_quantiles <- function(eta, typical, residual, times,
                                     residualRandom) {
  V <- typical[1L] * exp(eta[, 1L])
  CL <- typical[2L] * exp(eta[, 2L])
  response <- vapply(seq_along(times), function(j) {
    prediction <- 100 / V * exp(-(CL / V) * times[j])
    pmax(prediction * (1 + residual * residualRandom[, j]), 1e-10)
  }, numeric(nrow(eta)))
  apply(response, 2L, quantile, probs = c(.1, .5, .9), names = FALSE)
}

study_metric <- function(fitted, truth) {
  difference <- log(pmax(fitted, 1e-10)) - log(pmax(truth, 1e-10))
  sqrt(mean(difference^2))
}

study_regularize_correlation <- function(value, minimum = .05) {
  value <- (value + t(value)) / 2
  diag(value) <- 1
  smallest <- min(eigen(value, symmetric = TRUE, only.values = TRUE)$values)
  if (smallest < minimum) {
    weight <- (minimum - smallest) / (1 - smallest)
    value <- (1 - weight) * value + weight * diag(nrow(value))
  }
  cov2cor(value)
}

study_state_summary <- function(state) {
  last <- state$lastJoint
  list(
    margins = vapply(state$margins, `[[`, character(1), "name"),
    scales = copulaMarginScales(state$margins),
    correlation = copulaGaussianRvineCor(state$vine, length(state$margins)),
    projectionCount = last$projectionCount %||% NA_integer_,
    postFreezeBacktrackCount = last$postFreezeBacktrackCount %||% NA_integer_,
    postFreezeNoMoveCount = last$postFreezeNoMoveCount %||% NA_integer_,
    runtimeConditionsObserved = last$scoreTheory$runtimeConditionsObserved %||%
      NA)
}

study_run_one <- function(replicate, outputRoot, nSubjects = 300L,
                          nValidation = 30000L) {
  file <- file.path(outputRoot, sprintf("replicate_%03d.rds", replicate))
  if (file.exists(file)) return(file)
  started <- Sys.time()
  seed <- 880000L + 1000L * replicate
  truth <- study_truth()
  Rtruth <- study_correlation()
  structure <- study_structure()
  truthVine <- copulaGaussianRvineFromCor(Rtruth, structure)
  truthMargins <- list(copulaMarginNormal(truth$etaV),
    copulaMarginNormal(truth$etaCL),
    copulaMarginCovariateGamma(truth$gammaShape, truth$gammaScale),
    copulaMarginCovariateNormal(truth$albuminMean, truth$albuminSd))
  times <- c(.25, .5, 1, 2, 4, 8, 12)

  set.seed(seed)
  joint <- copulaMarginsQuantile(rvinecopulib::rvinecop(nSubjects, truthVine),
    truthMargins)
  colnames(joint) <- c("eta_V", "eta_CL", "CRP", "ALB")
  psi <- cbind(V = truth$V * exp(joint[, 1L]),
    CL = truth$CL * exp(joint[, 2L]))
  data <- data.frame(id = rep(seq_len(nSubjects), each = length(times)),
    dose = 100, time = rep(times, nSubjects))
  prediction <- study_pk(psi, data$id, cbind(data$dose, data$time))
  data$y <- pmax(prediction *
    (1 + truth$residual * rnorm(nrow(data))), 1e-8)
  conditioning <- joint[, c("CRP", "ALB"), drop = FALSE]
  sxdata <- saemixData(name.data = data, header = TRUE, name.group = "id",
    name.predictors = c("dose", "time"), name.response = "y", verbose = FALSE)

  control <- list(seed = seed + 100L, save = FALSE, save.graphs = FALSE,
    print = FALSE, displayProgress = FALSE, warnings = FALSE,
    nbiter.saemix = c(500L, 1000L), nbiter.mcmc = c(2L, 2L, 2L, 0L),
    ll.is = FALSE, fim = FALSE, map = FALSE)
  stockControl <- control
  stockControl$seed <- seed + 50L
  copulaClear()
  stockTime <- system.time(stockFit <- saemix(study_model(), sxdata,
    stockControl))["elapsed"]
  stockFixed <- as.numeric(stockFit@results@fixed.effects)[1:2]
  names(stockFixed) <- c("V", "CL")
  stockOmega <- stockFit@results@omega[1:2, 1:2, drop = FALSE]
  stockResidual <- study_residual(stockFit)

  standardMargins <- list(
    copulaMarginNormal(sqrt(stockOmega[1L, 1L])),
    copulaMarginNormal(sqrt(stockOmega[2L, 2L])),
    copulaFitCovariateMargin(conditioning[, 1L], "normal"),
    copulaFitCovariateMargin(conditioning[, 2L], "normal"))
  etaEbe <- stockFit@results@cond.mean.phi[, 1:2, drop = FALSE] -
    matrix(log(stockFixed), nrow = nSubjects, ncol = 2L, byrow = TRUE)
  covariateScore <- sapply(seq_len(2L), function(j) {
    margin <- standardMargins[[j + 2L]]
    qnorm(pmin(pmax(margin$cdf(conditioning[, j], margin$parameters),
      1e-10), 1 - 1e-10))
  })
  R0 <- study_regularize_correlation(cor(cbind(etaEbe, covariateScore)))
  standardPopulation <- copulaPopulation(
    copulaGaussianRvineFromCor(R0, structure), margins = standardMargins,
    scale = "transformed-additive",
    conditioning = list(values = conditioning,
      variableName = c("CRP", "ALB")),
    populationAlgorithm = "score-sa", scoreScale = .025,
    scoreGainScale = .2, scoreGainPower = .8, scoreGainOffset = 30,
    scoreBurn = 50L, scoreFiniteDifference = 1e-4, scoreProjection = 24)
  copulaClear()
  standardTime <- system.time(standardFit <- saemix(
    study_model(stockFixed, stockOmega, stockResidual), sxdata,
    control, population = standardPopulation))["elapsed"]
  standardState <- copulaGet(standardFit)

  fixedStart <- as.numeric(standardFit@results@fixed.effects)[1:2]
  names(fixedStart) <- c("V", "CL")
  omegaStart <- standardFit@results@omega[1:2, 1:2, drop = FALSE]
  residualStart <- study_residual(standardFit)
  flexibleMargins <- list(
    copulaMarginNormal(copulaMarginScales(standardState$margins)[1L]),
    copulaMarginNormal(copulaMarginScales(standardState$margins)[2L]),
    copulaFitCovariateMargin(conditioning[, 1L], "gamma"),
    copulaFitCovariateMargin(conditioning[, 2L], "normal"))
  flexiblePopulation <- copulaPopulation(
    copulaGaussianRvineFromCor(
      copulaGaussianRvineCor(standardState$vine, 4L), structure),
    margins = flexibleMargins, scale = "transformed-additive",
    conditioning = list(values = conditioning,
      variableName = c("CRP", "ALB")),
    populationAlgorithm = "score-sa", scoreScale = .025,
    scoreGainScale = .2, scoreGainPower = .8, scoreGainOffset = 30,
    scoreBurn = 50L, scoreFiniteDifference = 1e-4, scoreProjection = 24)
  flexibleControl <- control
  flexibleControl$seed <- seed + 200L
  copulaClear()
  flexibleTime <- system.time(flexibleFit <- saemix(
    study_model(fixedStart, omegaStart, residualStart), sxdata,
    flexibleControl, population = flexiblePopulation))["elapsed"]
  flexibleState <- copulaGet(flexibleFit)

  truthState <- list(vine = truthVine, margins = truthMargins,
    dEta = 2L, dConditioning = 2L,
    conditioningName = c("CRP", "ALB"), variableName = c("V", "CL"))
  ffem <- list(truth = copulaFremToFfem(truthState, "CRP"),
    gaussian = copulaFremToFfem(standardState, "CRP"),
    flexible = copulaFremToFfem(flexibleState, "CRP"))
  probability <- seq(.01, .999, length.out = 160L)
  crp <- qgamma(probability, truth$gammaShape, scale = truth$gammaScale)
  newdata <- matrix(crp, ncol = 1L, dimnames = list(NULL, "CRP"))
  typicalCL <- c(truth = truth$CL,
    gaussian = standardFit@results@fixed.effects[2L],
    flexible = flexibleFit@results@fixed.effects[2L])
  curve <- do.call(rbind, lapply(names(ffem), function(arm) data.frame(
    replicate = replicate, arm = arm, probability = probability, CRP = crp,
    medianCL = typicalCL[arm] * exp(
      copulaFfemLocation(ffem[[arm]], newdata)[, "CL"]))))

  set.seed(881337L)
  tailProbability <- runif(nValidation, .95, .9995)
  tailCRP <- qgamma(tailProbability, truth$gammaShape,
    scale = truth$gammaScale)
  etaRandom <- matrix(rnorm(nValidation * 2L), nValidation, 2L)
  residualRandom <- matrix(rnorm(nValidation * length(times)),
    nValidation, length(times))
  truthEta <- study_conditional_draw(tailCRP, truthVine, truthMargins,
    etaRandom)
  standardEta <- study_conditional_draw(tailCRP, standardState$vine,
    standardState$margins, etaRandom)
  flexibleEta <- study_conditional_draw(tailCRP, flexibleState$vine,
    flexibleState$margins, etaRandom)
  truthQ <- study_response_quantiles(truthEta, c(truth$V, truth$CL),
    truth$residual, times, residualRandom)
  standardQ <- study_response_quantiles(standardEta,
    as.numeric(standardFit@results@fixed.effects)[1:2],
    study_residual(standardFit), times, residualRandom)
  flexibleQ <- study_response_quantiles(flexibleEta,
    as.numeric(flexibleFit@results@fixed.effects)[1:2],
    study_residual(flexibleFit), times, residualRandom)
  qnames <- c("q10", "q50", "q90")
  vpc <- do.call(rbind, lapply(c("truth", "gaussian", "flexible"), function(arm) {
    value <- switch(arm, truth = truthQ, gaussian = standardQ, flexible = flexibleQ)
    data.frame(replicate = replicate, arm = arm,
      quantile = rep(qnames, each = length(times)),
      time = rep(times, 3L), value = as.vector(t(value)))
  }))
  metrics <- data.frame(replicate = replicate,
    arm = c("gaussian", "flexible"),
    log_rmse = c(study_metric(standardQ, truthQ),
      study_metric(flexibleQ, truthQ)))
  estimates <- data.frame(replicate = replicate,
    arm = c("gaussian", "flexible"),
    V = c(standardFit@results@fixed.effects[1L],
      flexibleFit@results@fixed.effects[1L]),
    CL = c(standardFit@results@fixed.effects[2L],
      flexibleFit@results@fixed.effects[2L]),
    residual = c(study_residual(standardFit),
      study_residual(flexibleFit)),
    runtime_seconds = c(standardTime, flexibleTime))
  result <- list(schema = 1L, replicate = replicate, seed = seed,
    nSubjects = nSubjects, nValidation = nValidation,
    truth = truth, curve = curve, vpc = vpc, metrics = metrics,
    estimates = estimates, standard = study_state_summary(standardState),
    flexible = study_state_summary(flexibleState),
    initializationCorrelation = R0, stockRuntimeSeconds = stockTime,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
  dir.create(outputRoot, recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(file, ".tmp-", Sys.getpid())
  saveRDS(result, temporary)
  if (!file.rename(temporary, file)) stop("could not atomically write ", file)
  copulaClear()
  file
}
