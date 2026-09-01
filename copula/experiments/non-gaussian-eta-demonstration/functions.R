ng_pk <- function(psi, id, xidep) {
  dose <- xidep[, 1L]; time <- xidep[, 2L]
  V <- psi[id, 1L]; CL <- psi[id, 2L]
  dose / V * exp(-CL / V * time)
}

ng_model <- function(fixed = c(V = 19, CL = 3.2), residual = .12,
                     omega = diag(c(.25^2, .40^2))) saemixModel(
  model = ng_pk, modeltype = "structural",
  description = "Non-Gaussian clearance eta demonstration",
  psi0 = matrix(as.numeric(fixed), nrow = 1L,
    dimnames = list(NULL, c("V", "CL"))),
  transform.par = c(1, 1), covariance.model = matrix(1, 2L, 2L),
  omega.init = omega, error.model = "proportional",
  error.init = c(0, residual), verbose = FALSE)

ng_data <- function(data) saemixData(
  name.data = data, header = TRUE, name.group = "id",
  name.predictors = c("dose", "time"), name.response = "y",
  verbose = FALSE)

ng_simulate <- function(seed, n = 180L) {
  set.seed(seed)
  truth <- list(V = 20, CL = 3.5, etaV = .22, etaCL = .45,
    etaCLShape = 2.5, rho = .35, residual = .12, dose = 100,
    times = c(.1, .5, 1, 2, 4, 8, 16, 24))
  margins <- list(copulaNaturalMarginLognormal(truth$etaV),
    copulaNaturalMarginGamma(truth$etaCLShape))
  z1 <- rnorm(n)
  z2 <- truth$rho * z1 + sqrt(1 - truth$rho^2) * rnorm(n)
  u <- cbind(pnorm(z1), pnorm(z2))
  typical <- matrix(rep(c(truth$V, truth$CL), each = n), nrow = n)
  psi <- copulaNaturalMarginsQuantile(u, typical, margins)
  colnames(psi) <- c("V", "CL")
  data <- data.frame(id = rep(seq_len(n), each = length(truth$times)),
    dose = truth$dose, time = rep(truth$times, n))
  prediction <- ng_pk(psi, data$id, cbind(data$dose, data$time))
  data$y <- pmax(prediction * (1 + truth$residual * rnorm(nrow(data))), 1e-8)
  list(data = data, psi = psi, truth = truth, margins = margins)
}

ng_control <- function(seed) list(seed = seed, save = FALSE,
  save.graphs = FALSE, print = FALSE, displayProgress = FALSE,
  warnings = FALSE, nbiter.saemix = c(1000, 500),
  nbiter.mcmc = c(2, 2, 2, 0), ll.is = FALSE, fim = FALSE, map = FALSE)

ng_fit_pair <- function(simulation, seedBase) {
  gaussianPopulation <- gaussianCopulaFrem(
    etaSd = c(.25, .40),
    correlation = matrix(c(1, .10, .10, 1), 2L, 2L),
    scoreBurn = 75L, gainScale = .18, gainPower = .80)
  gaussianTime <- system.time(gaussian <- saemix(
    ng_model(), ng_data(simulation$data), ng_control(seedBase + 11L),
    population = gaussianPopulation))["elapsed"]
  selectMargins <- function(draws, maxIter, seed) suppressWarnings(
    copulaSelectParameterMargins(
    gaussian, supports = c("positive", "positive"),
    candidates = list("lognormal", c("lognormal", "gamma", "weibull")),
    posteriorDraws = draws, max.iter = maxIter, seed = seed,
    optimizerMaxit = 200L, minimumEssFraction = .005,
    minimumPosteriorEssFraction = .01, maximumMcse = 1))
  selectionTime <- system.time(selection <- selectMargins(
    700L, 450L, seedBase + 21L))["elapsed"]
  if (!isTRUE(selection$diagnostics$selectionResolved)) {
    retryTime <- system.time(selection <- selectMargins(
      2000L, 700L, seedBase + 22L))["elapsed"]
    selectionTime <- selectionTime + retryTime
  }
  retainStandard <- !isTRUE(selection$diagnostics$selectionResolved)
  selection$retainedStandard <- retainStandard
  selectedFamilies <- if (retainStandard) c("lognormal", "lognormal") else
    selection$families
  startFixed <- as.numeric(gaussian@results@fixed.effects)[1:2]
  names(startFixed) <- c("V", "CL")
  startResidual <- ng_residual(gaussian)
  startOmega <- as.matrix(gaussian@results@omega[1:2, 1:2, drop = FALSE])
  if (any(!is.finite(startOmega)) ||
      inherits(try(chol(startOmega), silent = TRUE), "try-error"))
    startOmega <- diag(c(.25^2, .40^2))
  if (!retainStandard) {
    flexiblePopulation <- selection$population
    flexibleTime <- system.time(flexible <- saemix(
      ng_model(startFixed, startResidual, startOmega),
      ng_data(simulation$data), ng_control(seedBase + 31L),
      population = flexiblePopulation))["elapsed"]
  } else { flexible <- gaussian; flexibleTime <- 0 }
  gaussian@options$nmc.is <- 2500L
  gaussian <- suppressWarnings(llisCopula.saemix(gaussian,
    defensive = .20, batch = 100L, seed = seedBase + 41L))
  if (retainStandard) flexible <- gaussian else {
    flexible@options$nmc.is <- 2500L
    flexible <- suppressWarnings(llisCopula.saemix(flexible,
      defensive = .20, batch = 100L, seed = seedBase + 42L))
  }
  list(gaussian = gaussian, flexible = flexible, selection = selection,
    selectedFamilies = selectedFamilies,
    elapsed = c(gaussian = unname(gaussianTime),
      selection = unname(selectionTime), flexible = unname(flexibleTime)))
}

ng_residual <- function(fit) {
  value <- as.numeric(fit@results@respar[fit@results@indx.res])
  if (length(value) != 1L || !is.finite(value) || value <= 0)
    stop("could not identify fitted proportional residual SD")
  value
}

ng_parameter_draws <- function(fit, n, seed) {
  set.seed(seed); state <- copulaGet(fit)
  u <- rvinecopulib::rvinecop(n, state$vine)
  fixed <- as.numeric(fit@results@fixed.effects)[seq_len(state$dEta)]
  if (identical(state$populationScale, "parameter")) {
    typical <- matrix(rep(fixed, each = n), nrow = n)
    return(copulaNaturalMarginsQuantile(u, typical, state$margins))
  }
  eta <- copulaMarginsQuantile(u, state$margins)
  sweep(exp(eta), 2L, fixed, "*")
}

ng_vpc <- function(psi, residual, truth, arm, seed) {
  set.seed(seed); times <- truth$times
  id <- rep(seq_len(nrow(psi)), each = length(times))
  time <- rep(times, nrow(psi))
  prediction <- ng_pk(psi, id, cbind(truth$dose, time))
  y <- pmax(prediction * (1 + residual * rnorm(length(prediction))), 1e-8)
  value <- aggregate(y ~ time, data.frame(time, y), function(x)
    quantile(x, c(.01, .10, .50, .90, .99), names = FALSE))
  matrixValue <- if (is.matrix(value$y)) value$y else do.call(rbind, value$y)
  data.frame(time = value$time,
    probability = rep(c(.01, .10, .50, .90, .99), each = nrow(value)),
    value = as.vector(matrixValue), arm = arm)
}

ng_run_replicate <- function(replicate, resultRoot, force = FALSE,
                             nVpc = 100000L) {
  replicate <- as.integer(replicate)
  path <- file.path(resultRoot, sprintf("replicate_%03d.rds", replicate))
  if (file.exists(path) && !force) return(readRDS(path))
  seedBase <- 1060000L + 1000L * replicate
  simulation <- ng_simulate(seedBase + 1L)
  fits <- ng_fit_pair(simulation, seedBase)
  states <- lapply(fits[c("gaussian", "flexible")], copulaGet)
  stopifnot((isTRUE(fits$selection$retainedStandard) ||
      identical(states$flexible$lastJoint$scoreMethod,
        "hybrid-fixed-reference-path-score")),
    isTRUE(states$gaussian$lastJoint$scoreTheory$runtimeConditionsObserved),
    isTRUE(states$flexible$lastJoint$scoreTheory$runtimeConditionsObserved))

  set.seed(seedBase + 50L)
  z1 <- rnorm(nVpc); z2 <- simulation$truth$rho * z1 +
    sqrt(1 - simulation$truth$rho^2) * rnorm(nVpc)
  truthTypical <- matrix(rep(c(simulation$truth$V, simulation$truth$CL),
    each = nVpc), nrow = nVpc)
  truthPsi <- copulaNaturalMarginsQuantile(cbind(pnorm(z1), pnorm(z2)),
    truthTypical, simulation$margins)
  gaussianPsi <- ng_parameter_draws(fits$gaussian, nVpc, seedBase + 51L)
  flexiblePsi <- ng_parameter_draws(fits$flexible, nVpc, seedBase + 52L)
  fixedGaussian <- as.numeric(fits$gaussian@results@fixed.effects)[1:2]
  fixedFlexible <- as.numeric(fits$flexible@results@fixed.effects)[1:2]
  vpc <- rbind(
    ng_vpc(truthPsi, simulation$truth$residual, simulation$truth, "Generating model",
      seedBase + 61L),
    ng_vpc(gaussianPsi, ng_residual(fits$gaussian),
      simulation$truth, "Lognormal fit", seedBase + 62L),
    ng_vpc(flexiblePsi, ng_residual(fits$flexible),
      simulation$truth, "Free-margin fit", seedBase + 63L))

  grid <- seq(.2, 14, length.out = 400L)
  cutoff <- seq(3.5, 15, length.out = 250L)
  density <- do.call(rbind, lapply(c("Generating model", "Lognormal fit",
      "Free-margin fit"), function(arm) {
    margin <- switch(arm, "Generating model" = simulation$margins[[2L]],
      "Lognormal fit" = states$gaussian$margins[[2L]],
      states$flexible$margins[[2L]])
    typical <- switch(arm, "Generating model" = simulation$truth$CL,
      "Lognormal fit" = fixedGaussian[2L], fixedFlexible[2L])
    value <- if (inherits(margin, "saemix_natural_parameter_margin"))
      exp(margin$log_density(grid, rep(typical, length(grid)),
        margin$parameters)) else {
      etaGrid <- log(grid / typical)
      exp(margin$log_density(etaGrid, margin$parameters)) / grid
    }
    data.frame(parameter = grid, value = value, arm = arm)
  }))
  tail <- do.call(rbind, lapply(c("Generating model", "Lognormal fit",
      "Free-margin fit"), function(arm) {
    margin <- switch(arm, "Generating model" = simulation$margins[[2L]],
      "Lognormal fit" = states$gaussian$margins[[2L]],
      states$flexible$margins[[2L]])
    typical <- switch(arm, "Generating model" = simulation$truth$CL,
      "Lognormal fit" = fixedGaussian[2L], fixedFlexible[2L])
    value <- if (inherits(margin, "saemix_natural_parameter_margin"))
      1 - margin$cdf(cutoff, rep(typical, length(cutoff)), margin$parameters) else
      1 - margin$cdf(log(cutoff / typical), margin$parameters)
    data.frame(cutoff = cutoff, value = value, arm = arm)
  }))
  summary <- rbind(
    data.frame(arm = "Lognormal fit", family = "lognormal", shape = NA_real_,
      sd = stats::sd(gaussianPsi[, 2L]),
      tail_probability = mean(gaussianPsi[, 2L] > 8),
      V = fixedGaussian[1L], CL = fixedGaussian[2L],
      residual = ng_residual(fits$gaussian),
      log_likelihood = fits$gaussian@results@ll.is,
      likelihood_mcse = attr(fits$gaussian,
        "saemix.copula.likelihood")$se_loglik_total,
      runtime_seconds = fits$elapsed[["gaussian"]]),
    data.frame(arm = "Free-margin fit", family = fits$selectedFamilies[2L],
      shape = if ("shape" %in% names(states$flexible$margins[[2L]]$parameters))
        states$flexible$margins[[2L]]$parameters[["shape"]] else NA_real_,
      sd = stats::sd(flexiblePsi[, 2L]),
      tail_probability = mean(flexiblePsi[, 2L] > 8),
      V = fixedFlexible[1L], CL = fixedFlexible[2L],
      residual = ng_residual(fits$flexible),
      log_likelihood = fits$flexible@results@ll.is,
      likelihood_mcse = attr(fits$flexible,
        "saemix.copula.likelihood")$se_loglik_total,
      runtime_seconds = fits$elapsed[["flexible"]]))
  summary$selection_seconds <- c(0, fits$elapsed[["selection"]])
  summary$selected_families <- c("lognormal/lognormal",
    paste(fits$selectedFamilies, collapse = "/"))
  summary$retained_standard <- c(FALSE,
    isTRUE(fits$selection$retainedStandard))
  summary$replicate <- replicate
  result <- list(schema = 1L, replicate = replicate,
    data_seed = seedBase + 1L, summary = summary, density = density,
    tail = tail, vpc = vpc)
  dir.create(resultRoot, recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, path)
  result
}
