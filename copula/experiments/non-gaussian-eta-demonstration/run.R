## Demonstration of why a non-Gaussian eta margin matters.
##
## Data are generated from a one-compartment IV model with a Normal eta on V
## and a centred Gamma eta on log(CL).  We compare a fixed Gaussian-copula
## model with Normal eta margins against a freshly fitted model with the Gamma
## clearance margin.  Both fits use the same response model and Gaussian
## copula dependence.

suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
  library(ggplot2)
  library(gridExtra)
})

root <- "C:/package/saemix-copula/copula/experiments/non-gaussian-eta-demonstration"
out <- file.path(root, "out")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
force <- "--force" %in% commandArgs(trailingOnly = TRUE)

pk <- function(psi, id, xidep) {
  dose <- xidep[, 1L]
  time <- xidep[, 2L]
  V <- psi[id, 1L]
  CL <- psi[id, 2L]
  dose / V * exp(-CL / V * time)
}

make_model <- function() saemixModel(
  model = pk, modeltype = "structural",
  description = "Non-Gaussian clearance eta demonstration",
  psi0 = matrix(c(V = 19, CL = 3.2), nrow = 1L,
    dimnames = list(NULL, c("V", "CL"))),
  transform.par = c(1, 1), covariance.model = matrix(1, 2L, 2L),
  omega.init = diag(c(.25^2, .40^2)), error.model = "proportional",
  error.init = c(0, .12), verbose = FALSE)

make_data <- function(data) saemixData(
  name.data = data, header = TRUE, name.group = "id",
  name.predictors = c("dose", "time"), name.response = "y",
  verbose = FALSE)

simulate_data <- function(seed = 940101L, n = 180L) {
  set.seed(seed)
  truth <- list(V = 20, CL = 3.5, etaV = .22, etaCL = .45,
    etaCLShape = 2.5, rho = .35, residual = .12, dose = 100,
    times = c(.1, .5, 1, 2, 4, 8, 16, 24))
  etaCLMargin <- copulaMarginCenteredGamma(truth$etaCLShape, truth$etaCL)
  z1 <- rnorm(n)
  z2 <- truth$rho * z1 + sqrt(1 - truth$rho^2) * rnorm(n)
  eta <- cbind(eta_V = truth$etaV * z1,
    eta_CL = etaCLMargin$quantile(pnorm(z2), etaCLMargin$parameters))
  psi <- cbind(V = truth$V * exp(eta[, 1L]),
    CL = truth$CL * exp(eta[, 2L]))
  data <- data.frame(id = rep(seq_len(n), each = length(truth$times)),
    dose = truth$dose, time = rep(truth$times, n))
  prediction <- pk(psi, data$id, cbind(data$dose, data$time))
  data$y <- pmax(prediction * (1 + truth$residual * rnorm(nrow(data))), 1e-8)
  list(data = data, eta = eta, psi = psi, truth = truth,
    etaCLMargin = etaCLMargin)
}

control <- function(seed) list(seed = seed, save = FALSE,
  save.graphs = FALSE, print = FALSE, displayProgress = FALSE,
  warnings = FALSE, nbiter.saemix = c(1000, 500),
  nbiter.mcmc = c(2, 2, 2, 0), ll.is = FALSE, fim = FALSE, map = FALSE)

fit_models <- function(simulation) {
  cache <- file.path(out, "fits.rds")
  if (file.exists(cache) && !force) return(readRDS(cache))
  gaussian_population <- gaussianCopulaFrem(
    etaMargins = list(copulaMarginNormal(.25), copulaMarginNormal(.40)),
    correlation = matrix(c(1, .10, .10, 1), 2L, 2L),
    scoreBurn = 75L, gainScale = .18, gainPower = .80)
  gaussian_time <- system.time(gaussian <- saemix(
    make_model(), make_data(simulation$data), control(940111L),
    population = gaussian_population))["elapsed"]

  flexible_population <- gaussianCopulaFrem(
    etaMargins = list(copulaMarginNormal(.25),
      copulaMarginCenteredGamma(shape = 4, sd = .40)),
    correlation = matrix(c(1, .10, .10, 1), 2L, 2L),
    scoreBurn = 75L, gainScale = .18, gainPower = .80)
  flexible_time <- system.time(flexible <- saemix(
    make_model(), make_data(simulation$data), control(940131L),
    population = flexible_population))["elapsed"]

  gaussian@options$nmc.is <- 2500L
  flexible@options$nmc.is <- 2500L
  gaussian <- suppressWarnings(llisCopula.saemix(gaussian,
    defensive = .20, batch = 100L, seed = 940141L))
  flexible <- suppressWarnings(llisCopula.saemix(flexible,
    defensive = .20, batch = 100L, seed = 940142L))
  result <- list(gaussian = gaussian, flexible = flexible,
    elapsed = c(gaussian = gaussian_time, flexible = flexible_time))
  saveRDS(result, cache)
  result
}

margin_density <- function(fit, coordinate, grid) {
  margin <- copulaGet(fit)$margins[[coordinate]]
  exp(margin$log_density(grid, margin$parameters))
}

population_draws <- function(fit, n, seed) {
  set.seed(seed)
  state <- copulaGet(fit)
  copulaMarginsQuantile(rvinecopulib::rvinecop(n, state$vine),
    state$margins)[, seq_len(state$dEta), drop = FALSE]
}

fitted_residual <- function(fit) {
  value <- as.numeric(fit@results@respar[fit@results@indx.res])
  if (length(value) != 1L || !is.finite(value) || value <= 0)
    stop("could not identify fitted proportional residual SD")
  value
}

vpc_summary <- function(eta, fixed, residual, truth, arm, seed) {
  set.seed(seed)
  times <- truth$times
  psi <- sweep(exp(eta), 2L, fixed, "*")
  id <- rep(seq_len(nrow(eta)), each = length(times))
  time <- rep(times, nrow(eta))
  prediction <- pk(psi, id, cbind(truth$dose, time))
  y <- pmax(prediction * (1 + residual * rnorm(length(prediction))), 1e-8)
  value <- aggregate(y ~ time, data.frame(time, y), function(x)
    quantile(x, c(.01, .10, .50, .90, .99), names = FALSE))
  matrix_value <- if (is.matrix(value$y)) value$y else do.call(rbind, value$y)
  data.frame(time = value$time, probability = rep(c(.01, .10, .50, .90, .99),
    each = nrow(value)), value = as.vector(matrix_value), arm = arm)
}

simulation <- simulate_data()
fits <- fit_models(simulation)
states <- lapply(fits[c("gaussian", "flexible")], copulaGet)
stopifnot(identical(states$flexible$lastJoint$scoreMethod,
    "hybrid-fixed-reference-path-score"),
  isTRUE(states$flexible$lastJoint$scoreTheory$runtimeConditionsObserved))

grid <- seq(-1.65, 1.65, length.out = 600L)
density <- rbind(
  data.frame(eta = grid,
    density = exp(simulation$etaCLMargin$log_density(grid,
      simulation$etaCLMargin$parameters)), arm = "Generating model"),
  data.frame(eta = grid, density = margin_density(fits$gaussian, 2L, grid),
    arm = "Gaussian fit"),
  data.frame(eta = grid, density = margin_density(fits$flexible, 2L, grid),
    arm = "Flexible-margin fit"))

cutoff <- seq(.02, 1.5, length.out = 300L)
tail_curve <- function(margin, arm) data.frame(cutoff = cutoff,
  probability = margin$cdf(-cutoff, margin$parameters) +
    1 - margin$cdf(cutoff, margin$parameters), arm = arm)
tailProbability <- rbind(
  data.frame(cutoff = cutoff,
    probability = simulation$etaCLMargin$cdf(-cutoff,
      simulation$etaCLMargin$parameters) + 1 -
      simulation$etaCLMargin$cdf(cutoff,
        simulation$etaCLMargin$parameters),
    arm = "Generating model"),
  tail_curve(states$gaussian$margins[[2L]], "Gaussian fit"),
  tail_curve(states$flexible$margins[[2L]], "Flexible-margin fit"))

n_vpc <- 150000L
set.seed(940160L)
z1 <- rnorm(n_vpc)
z2 <- simulation$truth$rho * z1 + sqrt(1 - simulation$truth$rho^2) *
  rnorm(n_vpc)
truth_eta <- cbind(simulation$truth$etaV * z1,
  simulation$etaCLMargin$quantile(pnorm(z2),
    simulation$etaCLMargin$parameters))
gaussian_eta <- population_draws(fits$gaussian, n_vpc, 940161L)
flexible_eta <- population_draws(fits$flexible, n_vpc, 940162L)
fixed_gaussian <- as.numeric(fits$gaussian@results@fixed.effects)[1:2]
fixed_flexible <- as.numeric(fits$flexible@results@fixed.effects)[1:2]
vpc <- rbind(
  vpc_summary(truth_eta, c(simulation$truth$V, simulation$truth$CL),
    simulation$truth$residual, simulation$truth, "Generating model", 940171L),
  vpc_summary(gaussian_eta, fixed_gaussian,
    fitted_residual(fits$gaussian), simulation$truth,
    "Gaussian fit", 940172L),
  vpc_summary(flexible_eta, fixed_flexible,
    fitted_residual(fits$flexible), simulation$truth,
    "Flexible-margin fit", 940173L))

vpcTruth <- subset(vpc, arm == "Generating model",
  select = c("time", "probability", "value"))
names(vpcTruth)[3L] <- "truth"
vpcMetrics <- do.call(rbind, lapply(c("Gaussian fit", "Flexible-margin fit"),
  function(armName) {
    fitted <- merge(vpc[vpc$arm == armName, ], vpcTruth,
      by = c("time", "probability"))
    data.frame(arm = armName,
      all_quantile_log_rmse = sqrt(mean((log(fitted$value) -
        log(fitted$truth))^2)),
      extreme_quantile_log_rmse = sqrt(mean((log(fitted$value[
        fitted$probability %in% c(.01, .99)]) - log(fitted$truth[
          fitted$probability %in% c(.01, .99)]))^2)))
  }))

tail_probability <- function(draws, cutoff = 1) mean(abs(draws) > cutoff)
summary <- data.frame(
  arm = c("Generating model", "Gaussian fit", "Flexible-margin fit"),
  eta_CL_family = c("gamma-centered", states$gaussian$margins[[2L]]$name,
    states$flexible$margins[[2L]]$name),
  eta_CL_shape = c(simulation$truth$etaCLShape, NA_real_,
    states$flexible$margins[[2L]]$parameters[["shape"]]),
  eta_CL_sd = c(simulation$truth$etaCL,
    copulaMarginScales(states$gaussian$margins)[2L],
    copulaMarginScales(states$flexible$margins)[2L]),
  eta_CL_abs_gt_1 = c(tail_probability(truth_eta[, 2L]),
    tail_probability(gaussian_eta[, 2L]),
    tail_probability(flexible_eta[, 2L])),
  V = c(simulation$truth$V, fixed_gaussian[1L], fixed_flexible[1L]),
  CL = c(simulation$truth$CL, fixed_gaussian[2L], fixed_flexible[2L]),
  residual = c(simulation$truth$residual,
    fitted_residual(fits$gaussian), fitted_residual(fits$flexible)),
  log_likelihood = c(NA_real_, fits$gaussian@results@ll.is,
    fits$flexible@results@ll.is),
  likelihood_mcse = c(NA_real_,
    attr(fits$gaussian, "saemix.copula.likelihood")$se_loglik_total,
    attr(fits$flexible, "saemix.copula.likelihood")$se_loglik_total),
  score_method = c(NA_character_, states$gaussian$lastJoint$scoreMethod,
    states$flexible$lastJoint$scoreMethod),
  runtime_seconds = c(NA_real_, unname(fits$elapsed[1L]),
    unname(fits$elapsed[2L])))

colours <- c("Generating model" = "#222222", "Gaussian fit" = "#2878B5",
  "Flexible-margin fit" = "#E07A24")
p_density <- ggplot(density, aes(eta, density, colour = arm,
    linetype = arm)) +
  geom_line(linewidth = .95) +
  scale_colour_manual(values = colours) +
  scale_linetype_manual(values = c("Generating model" = "solid",
    "Gaussian fit" = "dashed", "Flexible-margin fit" = "dotdash")) +
  coord_cartesian(xlim = c(-1.5, 1.5)) +
  labs(title = "A  Population clearance variability",
    subtitle = "The centred Gamma eta is skewed but retains E(eta) = 0",
    x = expression(eta[CL]), y = "Density", colour = NULL, linetype = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom",
    panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

p_tail <- ggplot(tailProbability, aes(cutoff, probability, colour = arm,
    linetype = arm)) +
  geom_line(linewidth = .95) +
  scale_y_log10() +
  scale_colour_manual(values = colours) +
  scale_linetype_manual(values = c("Generating model" = "solid",
    "Gaussian fit" = "dashed", "Flexible-margin fit" = "dotdash")) +
  labs(title = "B  Probability of an extreme clearance effect",
    subtitle = expression(P(abs(eta[CL]) > x)~"reveals the tail deficit"),
    x = expression("Threshold "*x), y = "Exceedance probability (log scale)",
    colour = NULL, linetype = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom",
    panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

vpc$probability <- factor(vpc$probability,
  levels = c(.01, .10, .50, .90, .99),
  labels = c("1st", "10th", "Median", "90th", "99th"))
p_vpc <- ggplot(vpc, aes(time, value, colour = arm, group =
    interaction(arm, probability), linetype = probability)) +
  geom_line(linewidth = .9) +
  scale_y_log10() +
  scale_colour_manual(values = colours) +
  scale_linetype_manual(values = c("1st" = "dotted", "10th" = "dashed",
    "Median" = "solid", "90th" = "dashed", "99th" = "dotted")) +
  labs(title = "C  Population predictive distribution",
    subtitle = "Extreme concentration quantiles expose the eta-tail assumption",
    x = "Time", y = "Concentration (log scale)", colour = NULL,
    linetype = "Quantile") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom",
    panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

png(file.path(out, "non_gaussian_eta_demonstration.png"), width = 2600,
  height = 1850, res = 260)
grid.arrange(p_density, p_tail, p_vpc, ncol = 2L,
  layout_matrix = rbind(c(1, 2), c(3, 3)), heights = c(1, 1.18))
dev.off()

write.csv(summary, file.path(out, "fit_summary.csv"), row.names = FALSE)
write.csv(vpc, file.path(out, "vpc_quantiles.csv"), row.names = FALSE)
write.csv(vpcMetrics, file.path(out, "vpc_metrics.csv"), row.names = FALSE)
write.csv(tailProbability, file.path(out, "eta_tail_probabilities.csv"),
  row.names = FALSE)
saveRDS(list(simulation = simulation, fits = fits, density = density,
  tailProbability = tailProbability, vpc = vpc, vpcMetrics = vpcMetrics,
  summary = summary), file.path(out, "analysis.rds"))

print(summary)
print(vpcMetrics)
