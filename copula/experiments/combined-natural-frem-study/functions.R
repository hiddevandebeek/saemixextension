combined_truth <- function() list(V = 20, CL = 3.5, sdlogV = .22,
  gammaCL = 2.5, gammaCRP = 2, scaleCRP = 20, residual = .12,
  rhoVCL = .20, rhoVCRP = .15, rhoCLCRP = .60, dose = 100,
  times = c(.1, .5, 1, 2, 4, 8, 16, 24))

combined_correlation <- function(truth = combined_truth()) matrix(c(
  1, truth$rhoVCL, truth$rhoVCRP,
  truth$rhoVCL, 1, truth$rhoCLCRP,
  truth$rhoVCRP, truth$rhoCLCRP, 1), 3L, 3L, byrow = TRUE)

combined_structure <- function() rvinecopulib::cvine_structure(c(3, 1, 2))

combined_pk <- function(psi, id, xidep) {
  V <- psi[id, 1L]; CL <- psi[id, 2L]
  xidep[, 1L] / V * exp(-(CL / V) * xidep[, 2L])
}

combined_model <- function(fixed = c(V = 19, CL = 3.2), residual = .12,
                           omega = diag(c(.25^2, .40^2))) saemixModel(
  model = combined_pk, modeltype = "structural",
  description = "Combined natural-parameter and covariate-margin study",
  psi0 = matrix(as.numeric(fixed), 1L,
    dimnames = list(NULL, c("V", "CL"))), transform.par = c(1, 1),
  covariance.model = matrix(1, 2L, 2L), omega.init = omega,
  error.model = "proportional", error.init = c(0, residual), verbose = FALSE)

combined_data <- function(data) saemixData(name.data = data, header = TRUE,
  name.group = "id", name.predictors = c("dose", "time"),
  name.response = "y", verbose = FALSE)

combined_truth_margins <- function(truth = combined_truth()) list(
  copulaNaturalMarginLognormal(truth$sdlogV),
  copulaNaturalMarginGamma(truth$gammaCL),
  copulaMarginCovariateGamma(truth$gammaCRP, truth$scaleCRP))

combined_simulate <- function(seed, n = 250L) {
  set.seed(seed); truth <- combined_truth(); R <- combined_correlation(truth)
  z <- matrix(rnorm(n * 3L), n, 3L) %*% chol(R)
  margins <- combined_truth_margins(truth)
  typical <- matrix(rep(c(truth$V, truth$CL), each = n), nrow = n)
  psi <- copulaNaturalMarginsQuantile(pnorm(z[, 1:2, drop = FALSE]),
    typical, margins[1:2])
  crp <- margins[[3L]]$quantile(pnorm(z[, 3L]), margins[[3L]]$parameters)
  data <- data.frame(id = rep(seq_len(n), each = length(truth$times)),
    dose = truth$dose, time = rep(truth$times, n))
  prediction <- combined_pk(psi, data$id, cbind(data$dose, data$time))
  data$y <- pmax(prediction * (1 + truth$residual * rnorm(nrow(data))), 1e-8)
  list(data = data, psi = psi, crp = crp, truth = truth, margins = margins,
    correlation = R,
    vine = copulaGaussianRvineFromCor(R, combined_structure()))
}

combined_control <- function(seed, iterations = 1500L) list(seed = seed,
  save = FALSE, save.graphs = FALSE, print = FALSE, displayProgress = FALSE,
  warnings = FALSE, nbiter.saemix = c(iterations - 500L, 500L),
  nbiter.mcmc = c(2, 2, 2, 0), ll.is = FALSE, fim = FALSE, map = FALSE)

combined_residual <- function(fit) {
  value <- as.numeric(fit@results@respar[fit@results@indx.res])
  if (length(value) != 1L || !is.finite(value) || value <= 0)
    stop("could not identify fitted proportional residual SD")
  value
}

combined_fit_start <- function(fit, fallback = diag(c(.25^2, .40^2))) {
  fixed <- as.numeric(fit@results@fixed.effects)[1:2]; names(fixed) <- c("V", "CL")
  omega <- as.matrix(fit@results@omega[1:2, 1:2, drop = FALSE])
  if (any(!is.finite(omega)) || inherits(try(chol(omega), silent = TRUE),
      "try-error")) omega <- fallback
  list(fixed = fixed, omega = omega, residual = combined_residual(fit))
}

combined_initial_vine <- function() {
  R <- matrix(c(1, .05, .05, .05, 1, .10, .05, .10, 1), 3L, 3L)
  copulaGaussianRvineFromCor(R, combined_structure())
}

combined_population <- function(vine, margins, crp, scale) copulaPopulation(
  vine, margins = margins, scale = scale,
  conditioning = list(values = matrix(crp, ncol = 1L,
    dimnames = list(NULL, "CRP"))), populationAlgorithm = "score-sa",
  scoreScale = "auto", scoreBurn = 100L, scoreGainScale = .18,
  scoreGainPower = .80, scoreGainOffset = 30, scoreFiniteDifference = 1e-4,
  scoreProjection = 24)

combined_standard_margins <- function(crp) list(copulaMarginNormal(.25),
  copulaMarginNormal(.40), copulaFitCovariateMargin(crp, "normal"))

combined_flexible_covariate_margin <- function(crp)
  copulaFitCovariateMargin(crp, c("normal", "lognormal", "gamma", "weibull"))

combined_natural_lognormal_population <- function(state, crp) {
  scales <- vapply(state$margins[1:2], function(m)
    unname(m$parameters[["sd"]]), numeric(1L))
  margins <- list(copulaNaturalMarginLognormal(scales[1L]),
    copulaNaturalMarginLognormal(scales[2L]), state$margins[[3L]])
  combined_population(state$vine, margins, crp, "parameter")
}

combined_fit_pair <- function(simulation, seedBase, iterations = 1500L,
                              screenDraws = 700L, screenRetry = 2000L,
                              likelihoodDraws = 2500L) {
  sxdata <- combined_data(simulation$data)
  standardPopulation <- combined_population(combined_initial_vine(),
    combined_standard_margins(simulation$crp), simulation$crp,
    "transformed-additive")
  standardTime <- system.time(standard <- saemix(combined_model(), sxdata,
    combined_control(seedBase + 11L, iterations),
    population = standardPopulation))["elapsed"]
  start <- combined_fit_start(standard)

  flexibleCovariate <- combined_flexible_covariate_margin(simulation$crp)
  standardState <- copulaGet(standard)
  incumbentMargins <- c(standardState$margins[1:2], list(flexibleCovariate))
  incumbentPopulation <- combined_population(standardState$vine,
    incumbentMargins, simulation$crp, "transformed-additive")
  incumbentTime <- system.time(incumbent <- saemix(
    combined_model(start$fixed, start$residual, start$omega), sxdata,
    combined_control(seedBase + 21L, iterations),
    population = incumbentPopulation))["elapsed"]

  select <- function(draws, maxIter, seed) suppressWarnings(
    copulaSelectParameterMargins(incumbent,
      supports = c("positive", "positive"),
      candidates = list("lognormal", c("lognormal", "gamma", "weibull")),
      posteriorDraws = draws, max.iter = maxIter, seed = seed,
      optimizerMaxit = 200L, minimumEssFraction = .005,
      minimumPosteriorEssFraction = .01, maximumMcse = 1))
  selectionTime <- system.time(selection <- select(screenDraws, 450L,
    seedBase + 31L))["elapsed"]
  if (!isTRUE(selection$diagnostics$selectionResolved)) {
    retry <- system.time(selection <- select(screenRetry, 700L,
      seedBase + 32L))["elapsed"]
    selectionTime <- selectionTime + retry
  }
  retainStandard <- !isTRUE(selection$diagnostics$selectionResolved)
  selectedFamilies <- if (retainStandard) c("lognormal", "lognormal") else
    selection$families
  finalPopulation <- if (retainStandard)
    combined_natural_lognormal_population(copulaGet(incumbent), simulation$crp) else
    selection$population
  finalStart <- combined_fit_start(incumbent)
  flexibleTime <- system.time(flexible <- saemix(
    combined_model(finalStart$fixed, finalStart$residual, finalStart$omega),
    sxdata, combined_control(seedBase + 41L, iterations),
    population = finalPopulation))["elapsed"]

  for (arm in c("standard", "flexible")) {
    fit <- get(arm); fit@options$nmc.is <- likelihoodDraws
    fit <- suppressWarnings(llisCopula.saemix(fit, defensive = 0,
      batch = 100L, seed = seedBase + if (arm == "standard") 51L else 52L))
    assign(arm, fit)
  }
  selection$retainedStandard <- retainStandard
  list(standard = standard, incumbent = incumbent, flexible = flexible,
    selection = selection, selectedFamilies = selectedFamilies,
    covariateFamily = flexibleCovariate$name,
    elapsed = c(standard = unname(standardTime),
      incumbent = unname(incumbentTime), selection = unname(selectionTime),
      flexible = unname(flexibleTime)))
}

combined_state_ok <- function(state) {
  last <- state$lastJoint
  isTRUE(last$scoreTheory$runtimeConditionsObserved) &&
    identical(last$postFreezeProjectionCount, 0L) &&
    identical(last$postFreezeBacktrackCount, 0L) &&
    identical(last$postFreezeNoMoveCount, 0L) &&
    is.finite(last$scoreAverageMax)
}

combined_parameter_draws <- function(state, fixed, n, seed,
                                     crp = NULL) {
  set.seed(seed); dEta <- state$dEta
  if (is.null(crp)) {
    z <- rvinecopulib::rvinecop(n, state$vine)[, seq_len(dEta), drop = FALSE]
  } else {
    crp <- rep_len(crp, n); R <- copulaGaussianRvineCor(state$vine, state$d)
    margin <- state$margins[[dEta + 1L]]
    zc <- qnorm(margin$cdf(crp, margin$parameters))
    regression <- R[seq_len(dEta), dEta + 1L, drop = FALSE]
    covariance <- R[seq_len(dEta), seq_len(dEta), drop = FALSE] -
      tcrossprod(regression)
    z <- outer(zc, as.numeric(regression)) +
      matrix(rnorm(n * dEta), n, dEta) %*% chol(covariance)
    z <- pnorm(z)
  }
  if (identical(state$populationScale, "parameter")) {
    typical <- matrix(rep(fixed, each = n), nrow = n)
    copulaNaturalMarginsQuantile(z, typical, state$margins[seq_len(dEta)])
  } else {
    eta <- copulaMarginsQuantile(z, state$margins[seq_len(dEta)])
    sweep(exp(eta), 2L, fixed, "*")
  }
}

combined_response_vpc <- function(psi, residual, truth, arm, seed,
                                  probabilities = c(.01, .1, .5, .9, .99)) {
  set.seed(seed); times <- truth$times
  id <- rep(seq_len(nrow(psi)), each = length(times)); time <- rep(times, nrow(psi))
  prediction <- combined_pk(psi, id, cbind(truth$dose, time))
  y <- pmax(prediction * (1 + residual * rnorm(length(prediction))), 1e-8)
  value <- aggregate(y ~ time, data.frame(time, y), function(x)
    quantile(x, probabilities, names = FALSE))
  matrixValue <- if (is.matrix(value$y)) value$y else do.call(rbind, value$y)
  data.frame(time = value$time,
    probability = rep(probabilities, each = nrow(value)),
    value = as.vector(matrixValue), arm = arm)
}

combined_margin_density <- function(state, fixed, grid, coordinate = 2L) {
  margin <- state$margins[[coordinate]]
  if (identical(state$populationScale, "parameter"))
    return(exp(margin$log_density(grid, rep(fixed[coordinate], length(grid)),
      margin$parameters)))
  eta <- log(grid / fixed[coordinate])
  exp(margin$log_density(eta, margin$parameters)) / grid
}

combined_conditional_median <- function(state, fixed, crp) {
  dEta <- state$dEta; R <- copulaGaussianRvineCor(state$vine, state$d)
  marginC <- state$margins[[dEta + 1L]]
  zc <- qnorm(marginC$cdf(crp, marginC$parameters))
  zMedian <- outer(zc, as.numeric(R[seq_len(dEta), dEta + 1L]))
  if (identical(state$populationScale, "parameter")) {
    typical <- matrix(rep(fixed, each = length(crp)), nrow = length(crp))
    return(copulaNaturalMarginsQuantile(pnorm(zMedian), typical,
      state$margins[seq_len(dEta)])[, 2L])
  }
  eta <- copulaMarginsQuantile(pnorm(zMedian), state$margins[seq_len(dEta)])
  fixed[2L] * exp(eta[, 2L])
}

combined_run_replicate <- function(replicate, resultRoot, force = FALSE,
    nSubjects = 250L, nVpc = 50000L, iterations = 1500L,
    screenDraws = 700L, screenRetry = 2000L, likelihoodDraws = 2500L) {
  replicate <- as.integer(replicate)
  path <- file.path(resultRoot, sprintf("replicate_%03d.rds", replicate))
  if (file.exists(path) && !force) return(path)
  started <- Sys.time(); seedBase <- 1400000L + 1000L * replicate
  simulation <- combined_simulate(seedBase + 1L, nSubjects)
  fits <- combined_fit_pair(simulation, seedBase, iterations, screenDraws,
    screenRetry, likelihoodDraws)
  states <- lapply(fits[c("standard", "flexible")], copulaGet)
  if (!all(vapply(states, combined_state_ok, logical(1))))
    stop("post-freeze score-SA runtime conditions failed")
  fixed <- lapply(fits[c("standard", "flexible")], function(f)
    as.numeric(f@results@fixed.effects)[1:2])

  truthState <- list(vine = simulation$vine, margins = simulation$margins,
    d = 3L, dEta = 2L, dConditioning = 1L, populationScale = "parameter")
  set.seed(seedBase + 60L)
  truthPsi <- combined_parameter_draws(truthState,
    c(simulation$truth$V, simulation$truth$CL), nVpc, seedBase + 61L)
  standardPsi <- combined_parameter_draws(states$standard, fixed$standard,
    nVpc, seedBase + 62L)
  flexiblePsi <- combined_parameter_draws(states$flexible, fixed$flexible,
    nVpc, seedBase + 63L)
  populationVpc <- rbind(
    combined_response_vpc(truthPsi, simulation$truth$residual,
      simulation$truth, "Generating model", seedBase + 71L),
    combined_response_vpc(standardPsi, combined_residual(fits$standard),
      simulation$truth, "Gaussian FREM", seedBase + 72L),
    combined_response_vpc(flexiblePsi, combined_residual(fits$flexible),
      simulation$truth, "Flexible FREM", seedBase + 73L))

  set.seed(seedBase + 80L)
  tailProbability <- runif(nVpc, .95, .9995)
  tailCrp <- qgamma(tailProbability, simulation$truth$gammaCRP,
    scale = simulation$truth$scaleCRP)
  truthConditional <- combined_parameter_draws(truthState,
    c(simulation$truth$V, simulation$truth$CL), nVpc, seedBase + 81L, tailCrp)
  standardConditional <- combined_parameter_draws(states$standard,
    fixed$standard, nVpc, seedBase + 82L, tailCrp)
  flexibleConditional <- combined_parameter_draws(states$flexible,
    fixed$flexible, nVpc, seedBase + 83L, tailCrp)
  conditionalVpc <- rbind(
    combined_response_vpc(truthConditional, simulation$truth$residual,
      simulation$truth, "Generating model", seedBase + 91L, c(.1, .5, .9)),
    combined_response_vpc(standardConditional, combined_residual(fits$standard),
      simulation$truth, "Gaussian FREM", seedBase + 92L, c(.1, .5, .9)),
    combined_response_vpc(flexibleConditional, combined_residual(fits$flexible),
      simulation$truth, "Flexible FREM", seedBase + 93L, c(.1, .5, .9)))

  grid <- seq(.2, 14, length.out = 300L)
  density <- rbind(
    data.frame(CL = grid, density = combined_margin_density(truthState,
      c(simulation$truth$V, simulation$truth$CL), grid), arm = "Generating model"),
    data.frame(CL = grid, density = combined_margin_density(states$standard,
      fixed$standard, grid), arm = "Gaussian FREM"),
    data.frame(CL = grid, density = combined_margin_density(states$flexible,
      fixed$flexible, grid), arm = "Flexible FREM"))
  probability <- seq(.01, .999, length.out = 180L)
  crpGrid <- qgamma(probability, simulation$truth$gammaCRP,
    scale = simulation$truth$scaleCRP)
  relationship <- rbind(
    data.frame(probability, CRP = crpGrid,
      medianCL = combined_conditional_median(truthState,
        c(simulation$truth$V, simulation$truth$CL), crpGrid),
      arm = "Generating model"),
    data.frame(probability, CRP = crpGrid,
      medianCL = combined_conditional_median(states$standard,
        fixed$standard, crpGrid), arm = "Gaussian FREM"),
    data.frame(probability, CRP = crpGrid,
      medianCL = combined_conditional_median(states$flexible,
        fixed$flexible, crpGrid), arm = "Flexible FREM"))

  likelihood <- do.call(rbind, lapply(c("standard", "flexible"), function(arm) {
    fit <- fits[[arm]]; state <- states[[arm]]
    covll <- sum(copulaGaussianFremConditioningLogDensity(
      matrix(simulation$crp, ncol = 1L), state$vine, state$margins, 2L))
    data.frame(arm = if (arm == "standard") "Gaussian FREM" else "Flexible FREM",
      total = fit@results@ll.is, covariate = covll,
      conditional_response = fit@results@ll.is - covll,
      mcse = attr(fit, "saemix.copula.likelihood")$se_loglik_total)
  }))
  summary <- data.frame(arm = c("Gaussian FREM", "Flexible FREM"),
    V = c(fixed$standard[1L], fixed$flexible[1L]),
    CL = c(fixed$standard[2L], fixed$flexible[2L]),
    residual = c(combined_residual(fits$standard), combined_residual(fits$flexible)),
    runtime = c(fits$elapsed[["standard"]], fits$elapsed[["flexible"]]),
    familyCL = c("lognormal", fits$selectedFamilies[2L]),
    familyCRP = c("normal", fits$covariateFamily),
    retainedStandard = c(FALSE, isTRUE(fits$selection$retainedStandard)),
    scoreAverageMax = vapply(states, function(s) s$lastJoint$scoreAverageMax,
      numeric(1L)))
  summary$selectionRuntime <- c(0, fits$elapsed[["selection"]])
  summary$incumbentRuntime <- c(0, fits$elapsed[["incumbent"]])
  summary$replicate <- replicate
  summary <- merge(summary, likelihood, by = "arm", sort = FALSE)
  result <- list(schema = 1L, replicate = replicate, dataSeed = seedBase + 1L,
    truth = simulation$truth, summary = summary, density = density,
    relationship = relationship, populationVpc = populationVpc,
    conditionalVpc = conditionalVpc,
    selectionTable = fits$selection$table,
    selectionDiagnostics = fits$selection$diagnostics,
    elapsedSeconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
  dir.create(resultRoot, recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid()); saveRDS(result, temporary)
  if (!file.rename(temporary, path)) stop("could not atomically write ", path)
  copulaClear(); path
}
