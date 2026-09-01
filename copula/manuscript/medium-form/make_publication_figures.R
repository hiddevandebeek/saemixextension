suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
  library(rvinecopulib)
  library(ggplot2)
  library(gridExtra)
})

repo <- "C:/package/saemix-copula"
summaryRoot <- file.path(repo,
  "copula/study_frem_validation_large/out/summary")
figureRoot <- file.path(repo, "copula/manuscript/medium-form/figures")
dir.create(figureRoot, recursive = TRUE, showWarnings = FALSE)

scenarioLabels <- c(normal_normal = "Normal / Normal",
  lognormal_normal = "Lognormal / Normal",
  gamma_normal = "Gamma / Normal",
  weibull_lognormal = "Weibull / Lognormal")
armLabels <- c(standard = "Gaussian FREM", flexible = "Flexible-margin FREM")
colours <- c("Gaussian FREM" = "#2878B5",
  "Flexible-margin FREM" = "#E07A24")

structural <- read.csv(file.path(summaryRoot,
  "structural_recovery_summary.csv"), stringsAsFactors = FALSE)
structural$scenario_label <- factor(scenarioLabels[structural$scenario],
  levels = scenarioLabels)
structural$arm_label <- factor(armLabels[structural$arm], levels = armLabels)
structural$parameter <- factor(structural$parameter,
  levels = c("V", "CL", "residual_sd"),
  labels = c("V", "CL", "Residual SD"))
structural$lower <- structural$relative_bias_percent -
  1.96 * 100 * structural$mcse_bias / structural$truth
structural$upper <- structural$relative_bias_percent +
  1.96 * 100 * structural$mcse_bias / structural$truth
pRecovery <- ggplot(structural, aes(parameter, relative_bias_percent,
    colour = arm_label, shape = arm_label)) +
  geom_hline(yintercept = 0, colour = "grey55") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = .08,
    position = position_dodge(.35), linewidth = .45) +
  geom_point(position = position_dodge(.35), size = 2) +
  facet_wrap(~scenario_label, ncol = 2) +
  scale_colour_manual(values = colours, name = NULL) +
  scale_shape_manual(values = c(16, 17), name = NULL) +
  labs(x = NULL, y = "Mean relative bias (%, 95% Monte Carlo CI)",
    title = "A  Structural and residual-error recovery") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 10))

decomposition <- read.csv(file.path(summaryRoot,
  "likelihood_decomposition.csv"), stringsAsFactors = FALSE)
longDecomposition <- rbind(
  data.frame(cell = decomposition$cell, scenario = decomposition$scenario,
    component = "Covariate density", value = decomposition$delta_covariate),
  data.frame(cell = decomposition$cell, scenario = decomposition$scenario,
    component = "Conditional response",
    value = decomposition$delta_conditional_response))
means <- aggregate(value ~ scenario + component, longDecomposition, mean)
means$scenario_label <- factor(scenarioLabels[means$scenario],
  levels = scenarioLabels)
means$component <- factor(means$component,
  levels = c("Covariate density", "Conditional response"))
jointIntervals <- do.call(rbind, lapply(split(decomposition,
  decomposition$scenario), function(x) data.frame(scenario = x$scenario[1L],
    lower = quantile(x$delta_joint, .1), median = median(x$delta_joint),
    upper = quantile(x$delta_joint, .9))))
jointIntervals$scenario_label <- factor(
  scenarioLabels[jointIntervals$scenario], levels = scenarioLabels)
pLikelihood <- ggplot(means, aes(scenario_label, value, fill = component)) +
  geom_hline(yintercept = 0, colour = "grey45") +
  geom_col(width = .62) +
  geom_errorbar(data = jointIntervals,
    aes(scenario_label, ymin = lower, ymax = upper), inherit.aes = FALSE,
    width = .16, colour = "#222222", linewidth = .55) +
  geom_point(data = jointIntervals, aes(scenario_label, median),
    inherit.aes = FALSE, shape = 21, fill = "white", colour = "#222222",
    size = 2.1, stroke = .6) +
  scale_fill_manual(values = c("Covariate density" = "#5B8DB8",
    "Conditional response" = "#E9A15B"), name = NULL) +
  labs(x = NULL, y = expression(Delta*" joint log likelihood"),
    title = "B  Source of the likelihood difference") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
    axis.text.x = element_text(angle = 24, hjust = 1),
    plot.title = element_text(face = "bold", size = 10))

png(file.path(figureRoot, "figure1_reestimation.png"), width = 2200,
  height = 2100, res = 240)
grid.arrange(pRecovery, pLikelihood, ncol = 1, heights = c(1.18, .82))
dev.off()

mapping <- read.csv(file.path(repo,
  "copula/examples/gamma-tail-vpc/out/ffem_translation_mapping.csv"),
  stringsAsFactors = FALSE)

## Recreate the simulated individuals used in the fitted example. The saved
## response data do not contain individual clearance because it is latent in an
## NLME analysis; these values are available only because this is a simulation.
example <- readRDS(file.path(repo,
  "copula/examples/gamma-tail-vpc/out/example_fit.rds"))
truth <- example$truth
truthMargins <- list(copulaMarginNormal(truth$etaV),
  copulaMarginNormal(truth$etaCL),
  copulaMarginCovariateGamma(truth$gammaShape, truth$gammaScale),
  copulaMarginCovariateNormal(truth$albuminMean, truth$albuminSd))
truthVine <- copulaGaussianRvineFromCor(example$Rtruth,
  cvine_structure(c(3, 1, 4, 2)))
set.seed(471103L)
joint <- copulaMarginsQuantile(rvinecop(300L, truthVine), truthMargins)
stopifnot(max(abs(joint[, 3:4] - example$conditioning)) < 1e-12)
individuals <- data.frame(CRP = joint[, 3L],
  CL = truth$CL * exp(joint[, 2L]))

## Independent algebraic verification of the plotted Gaussian FREM curve.
standardState <- example$standardState
Rfit <- copulaGaussianRvineCor(standardState$vine, standardState$d)
crpMargin <- standardState$margins[[3L]]
etaScale <- copulaMarginScales(standardState$margins)[2L]
gaussianRows <- mapping$arm == "Gaussian FREM / FFEM"
zDirect <- (mapping$CRP[gaussianRows] - crpMargin$parameters[["mean"]]) /
  crpMargin$parameters[["sd"]]
scoreError <- max(abs(mapping$score[gaussianRows] - zDirect))
medianDirect <- example$standardFit@results@fixed.effects[2L] *
  exp(etaScale * Rfit[2L, 3L] * zDirect)
medianError <- max(abs(mapping$medianCL[gaussianRows] - medianDirect))
write.csv(data.frame(check = c("Gaussian score", "Gaussian conditional median"),
  maximum_absolute_difference = c(scoreError, medianError)),
  file.path(figureRoot, "figure2_algebraic_verification.csv"), row.names = FALSE)
stopifnot(scoreError < 1e-10, medianError < 1e-10)

mapping$arm <- factor(mapping$arm,
  levels = c("Generating model", "Gaussian FREM / FFEM", "Translated FFEM"),
  labels = c("Generating model", "Gaussian FREM", "Flexible-margin FREM / FFEM"))
effectColours <- c("Generating model" = "#222222",
  "Gaussian FREM" = "#2878B5",
  "Flexible-margin FREM / FFEM" = "#E07A24")
effectTypes <- c("Generating model" = "solid", "Gaussian FREM" = "dashed",
  "Flexible-margin FREM / FFEM" = "solid")
pScore <- ggplot(mapping, aes(CRP, score, colour = arm, linetype = arm)) +
  geom_line(linewidth = .9) +
  geom_rug(data = individuals, aes(x = CRP), inherit.aes = FALSE,
    sides = "b", alpha = .16, colour = "#555555") +
  scale_colour_manual(values = effectColours, name = NULL) +
  scale_linetype_manual(values = effectTypes, name = NULL) +
  labs(x = "CRP", y = "Gaussian covariate score",
    title = "A  Marginal transformation") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "none",
    plot.title = element_text(face = "bold", size = 10))
pEffect <- ggplot(mapping,
    aes(CRP, medianCL, colour = arm, linetype = arm)) +
  geom_point(data = individuals, aes(CRP, CL), inherit.aes = FALSE,
    colour = "#777777", alpha = .18, size = .65) +
  geom_line(linewidth = .9) +
  scale_colour_manual(values = effectColours, name = NULL) +
  scale_linetype_manual(values = effectTypes, name = NULL) +
  labs(x = "CRP", y = "Conditional median clearance",
    title = "B  Conditional FFEM relationship") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 10))
png(file.path(figureRoot, "figure2_conditional_effect.png"), width = 2300,
  height = 1100, res = 240)
grid.arrange(pScore, pEffect, ncol = 2, widths = c(1, 1.12))
dev.off()

pEffectMain <- pEffect +
  labs(title = NULL, x = "CRP", y = "Conditional median clearance") +
  theme(plot.title = element_blank())
ggsave(file.path(figureRoot, "figure2_conditional_relationship.png"),
  pEffectMain, width = 7.4, height = 4.8, dpi = 300)
