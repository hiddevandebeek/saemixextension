args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1L] else "pilot"
if (!mode %in% c("pilot", "final")) stop("mode must be pilot or final")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })
root <- "C:/package/saemix-copula/copula/experiments/combined-natural-frem-study"
input <- file.path(root, if (mode == "pilot") "pilot" else "replicates")
output <- file.path(root, if (mode == "pilot") "pilot-out" else "out")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
files <- sort(list.files(input, "^replicate_[0-9]+\\.rds$", full.names = TRUE))
if (!length(files)) stop("no replicate files")
results <- lapply(files, readRDS)
replicates <- vapply(results, `[[`, integer(1), "replicate")
if (anyDuplicated(replicates)) stop("duplicate replicate identifiers")

bind <- function(name) bind_rows(lapply(results, `[[`, name))
summary <- bind("summary"); density <- bind("density")
relationship <- bind("relationship"); populationVpc <- bind("populationVpc")
conditionalVpc <- bind("conditionalVpc")
for (object in c("density", "relationship", "populationVpc", "conditionalVpc")) {
  value <- get(object)
  value$replicate <- rep(replicates, vapply(lapply(results, `[[`, object),
    nrow, integer(1)))
  assign(object, value)
}

quantile_frame <- function(data, groups, value = "value") data |>
  group_by(across(all_of(groups))) |>
  summarise(median = median(.data[[value]]), lower = quantile(.data[[value]], .025),
    upper = quantile(.data[[value]], .975), .groups = "drop")

metric <- function(data, truthArm = "Generating model") {
  truth <- data |> filter(arm == truthArm) |>
    select(replicate, time, probability, truth = value)
  data |> filter(arm != truthArm) |> left_join(truth,
    by = c("replicate", "time", "probability")) |>
    group_by(replicate, arm) |>
    summarise(log_rmse = sqrt(mean((log(value) - log(truth))^2)), .groups = "drop")
}
populationMetric <- metric(populationVpc)
conditionalMetric <- metric(conditionalVpc)
relationshipTruth <- relationship |> filter(arm == "Generating model") |>
  select(replicate, probability, truth = medianCL)
relationshipMetric <- relationship |> filter(arm != "Generating model") |>
  left_join(relationshipTruth, by = c("replicate", "probability")) |>
  group_by(replicate, arm) |>
  summarise(log_rmse = sqrt(mean((log(medianCL) - log(truth))^2)), .groups = "drop")

wideSummary <- summary |> select(replicate, arm, V, CL, residual, total,
  covariate, conditional_response, mcse, familyCL, familyCRP,
  retainedStandard, scoreAverageMax) |>
  pivot_wider(names_from = arm, values_from = -replicate)
quality <- data.frame(
  replicates = length(results),
  unique_data_seeds = length(unique(vapply(results, `[[`, integer(1), "dataSeed"))),
  gamma_parameter_selections = sum(summary$arm == "Flexible FREM" & summary$familyCL == "gamma"),
  gamma_covariate_selections = sum(summary$arm == "Flexible FREM" & summary$familyCRP == "gamma"),
  retained_standard = sum(summary$arm == "Flexible FREM" & summary$retainedStandard),
  likelihood_better = sum(wideSummary$`total_Flexible FREM` > wideSummary$`total_Gaussian FREM`),
  conditional_likelihood_better = sum(wideSummary$`conditional_response_Flexible FREM` >
    wideSummary$`conditional_response_Gaussian FREM`),
  likelihood_over_2mcse = sum((wideSummary$`total_Flexible FREM` -
    wideSummary$`total_Gaussian FREM`) > 2 * sqrt(
      wideSummary$`mcse_Flexible FREM`^2 + wideSummary$`mcse_Gaussian FREM`^2)),
  population_vpc_better = sum((populationMetric |> filter(arm == "Flexible FREM"))$log_rmse <
    (populationMetric |> filter(arm == "Gaussian FREM"))$log_rmse),
  conditional_vpc_better = sum((conditionalMetric |> filter(arm == "Flexible FREM"))$log_rmse <
    (conditionalMetric |> filter(arm == "Gaussian FREM"))$log_rmse),
  relationship_better = sum((relationshipMetric |> filter(arm == "Flexible FREM"))$log_rmse <
    (relationshipMetric |> filter(arm == "Gaussian FREM"))$log_rmse),
  median_total_gain = median(wideSummary$`total_Flexible FREM` - wideSummary$`total_Gaussian FREM`),
  median_conditional_gain = median(wideSummary$`conditional_response_Flexible FREM` -
    wideSummary$`conditional_response_Gaussian FREM`),
  max_score_average = max(summary$scoreAverageMax))
write.csv(summary, file.path(output, "fit_summary.csv"), row.names = FALSE)
write.csv(populationMetric, file.path(output, "population_vpc_metrics.csv"), row.names = FALSE)
write.csv(conditionalMetric, file.path(output, "conditional_vpc_metrics.csv"), row.names = FALSE)
write.csv(relationshipMetric, file.path(output, "relationship_metrics.csv"), row.names = FALSE)
write.csv(quality, file.path(output, "quality_summary.csv"), row.names = FALSE)
metricSummary <- bind_rows(
  populationMetric |> group_by(arm) |> summarise(endpoint = "population VPC",
    median_log_rmse = median(log_rmse), lower = quantile(log_rmse, .025),
    upper = quantile(log_rmse, .975), .groups = "drop"),
  conditionalMetric |> group_by(arm) |> summarise(endpoint = "conditional VPC",
    median_log_rmse = median(log_rmse), lower = quantile(log_rmse, .025),
    upper = quantile(log_rmse, .975), .groups = "drop"),
  relationshipMetric |> group_by(arm) |> summarise(endpoint = "conditional CL relationship",
    median_log_rmse = median(log_rmse), lower = quantile(log_rmse, .025),
    upper = quantile(log_rmse, .975), .groups = "drop"))
familySummary <- summary |> filter(arm == "Flexible FREM") |>
  count(familyCL, familyCRP, retainedStandard, name = "datasets")
estimateSummary <- summary |> group_by(arm) |>
  summarise(across(c(V, CL, residual), list(mean = mean, median = median,
    lower = ~quantile(.x, .025), upper = ~quantile(.x, .975))), .groups = "drop")
gainSummary <- data.frame(endpoint = c("joint likelihood", "conditional response likelihood"),
  median = c(median(wideSummary$`total_Flexible FREM` - wideSummary$`total_Gaussian FREM`),
    median(wideSummary$`conditional_response_Flexible FREM` -
      wideSummary$`conditional_response_Gaussian FREM`)),
  lower = c(quantile(wideSummary$`total_Flexible FREM` - wideSummary$`total_Gaussian FREM`, .025),
    quantile(wideSummary$`conditional_response_Flexible FREM` -
      wideSummary$`conditional_response_Gaussian FREM`, .025)),
  upper = c(quantile(wideSummary$`total_Flexible FREM` - wideSummary$`total_Gaussian FREM`, .975),
    quantile(wideSummary$`conditional_response_Flexible FREM` -
      wideSummary$`conditional_response_Gaussian FREM`, .975)))
write.csv(metricSummary, file.path(output, "metric_summary.csv"), row.names = FALSE)
write.csv(familySummary, file.path(output, "family_summary.csv"), row.names = FALSE)
write.csv(estimateSummary, file.path(output, "estimate_summary.csv"), row.names = FALSE)
write.csv(gainSummary, file.path(output, "likelihood_gain_summary.csv"), row.names = FALSE)

truthColour <- "#202020"; flexibleColour <- "#E76F00"; gaussianColour <- "#2878B5"
fittedColours <- c("Flexible FREM" = flexibleColour, "Gaussian FREM" = gaussianColour)
themePaper <- theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank(),
  legend.position = "bottom", plot.title = element_text(face = "bold", size = 13),
  plot.subtitle = element_text(size = 10), strip.background = element_rect(fill = "grey94"))

densitySummary <- quantile_frame(density |> filter(arm != "Generating model"),
  c("arm", "CL"), "density")
truthDensity <- density |> filter(arm == "Generating model") |>
  group_by(CL) |> summarise(density = median(density), .groups = "drop")
pA <- ggplot(densitySummary, aes(CL, median, colour = arm, fill = arm)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .14, colour = NA) +
  geom_line(linewidth = .8) + geom_line(data = truthDensity,
    aes(CL, density), inherit.aes = FALSE, colour = truthColour, linewidth = .9) +
  scale_colour_manual(values = fittedColours) + scale_fill_manual(values = fittedColours) +
  labs(title = "A  Natural clearance distribution", subtitle = "Median and 95% range across fitted datasets",
    x = "Clearance", y = "Density", colour = NULL, fill = NULL) + themePaper

relationshipSummary <- quantile_frame(relationship |> filter(arm != "Generating model"),
  c("arm", "probability", "CRP"), "medianCL")
truthRelationship <- relationship |> filter(arm == "Generating model") |>
  group_by(probability, CRP) |> summarise(medianCL = median(medianCL), .groups = "drop")
pB <- ggplot(relationshipSummary, aes(probability, median,
    colour = arm, fill = arm)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .14, colour = NA) +
  geom_line(linewidth = .8) + geom_line(data = truthRelationship,
    aes(probability, medianCL), inherit.aes = FALSE,
    colour = truthColour, linewidth = .9) +
  scale_colour_manual(values = fittedColours) + scale_fill_manual(values = fittedColours) +
  scale_x_continuous(trans = scales::logit_trans(),
    breaks = c(.01, .25, .5, .75, .95, .999),
    labels = c("1", "25", "50", "75", "95", "99.9")) +
  labs(title = "B  Conditional clearance relationship",
    subtitle = "Median CL given CRP percentile", x = "CRP percentile", y = "Conditional median CL",
    colour = NULL, fill = NULL) + themePaper

vpc_panel <- function(data, title, subtitle, probabilities) {
  fit <- quantile_frame(data |> filter(arm != "Generating model",
    probability %in% probabilities), c("arm", "time", "probability"))
  truth <- data |> filter(arm == "Generating model", probability %in% probabilities) |>
    group_by(time, probability) |> summarise(value = median(value), .groups = "drop")
  shapes <- setNames(c(17, 16, 15)[seq_along(probabilities)], probabilities)
  ggplot(fit, aes(time, median, colour = arm, fill = arm,
      group = interaction(arm, probability))) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .11, colour = NA) +
    geom_line(linewidth = .75, linetype = "dashed") +
    geom_line(data = truth, aes(time, value, group = probability),
      inherit.aes = FALSE, colour = truthColour, linewidth = .65) +
    geom_point(data = truth, aes(time, value, shape = factor(probability)),
      inherit.aes = FALSE, colour = truthColour, size = 2) +
    scale_colour_manual(values = fittedColours) + scale_fill_manual(values = fittedColours) +
    scale_shape_manual(values = shapes) + scale_y_log10() +
    guides(shape = "none") +
    labs(title = title, subtitle = subtitle, x = "Time", y = "Concentration (log scale)",
      colour = NULL, fill = NULL, shape = "Generating quantile") + themePaper
}
pC <- vpc_panel(populationVpc, "C  Population VPC",
  "1st, median and 99th percentiles", c(.01, .5, .99))
pD <- vpc_panel(conditionalVpc, "D  Conditional VPC",
  "Upper 5% of the CRP distribution", c(.1, .5, .9))

if (!requireNamespace("patchwork", quietly = TRUE)) stop("patchwork is required")
figure <- (pA + pB) / (pC + pD) + patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
ggsave(file.path(output, "combined_frem_figure.png"), figure,
  width = 13, height = 9.6, dpi = 220, bg = "white")
print(quality)
