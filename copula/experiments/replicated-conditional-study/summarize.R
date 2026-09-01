suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})
root <- "C:/package/saemix-copula/copula/experiments/replicated-conditional-study"
files <- sort(list.files(file.path(root, "results"), pattern = "[.]rds$",
  full.names = TRUE))
if (!length(files)) stop("no replicate results")
results <- lapply(files, readRDS)
stopifnot(all(vapply(results, `[[`, integer(1), "schema") == 1L))
replicateId <- vapply(results, `[[`, integer(1), "replicate")
stopifnot(!anyDuplicated(replicateId), all(is.finite(replicateId)))
diagnostics <- do.call(rbind, lapply(results, function(x) data.frame(
  replicate = x$replicate,
  standard_projection = x$standard$projectionCount,
  standard_backtrack = x$standard$postFreezeBacktrackCount,
  standard_no_move = x$standard$postFreezeNoMoveCount,
  standard_runtime_conditions = x$standard$runtimeConditionsObserved,
  flexible_projection = x$flexible$projectionCount,
  flexible_backtrack = x$flexible$postFreezeBacktrackCount,
  flexible_no_move = x$flexible$postFreezeNoMoveCount,
  flexible_runtime_conditions = x$flexible$runtimeConditionsObserved)))
stopifnot(all(diagnostics$standard_projection == 0L),
  all(diagnostics$standard_backtrack == 0L),
  all(diagnostics$standard_no_move == 0L),
  all(diagnostics$standard_runtime_conditions),
  all(diagnostics$flexible_projection == 0L),
  all(diagnostics$flexible_backtrack == 0L),
  all(diagnostics$flexible_no_move == 0L),
  all(diagnostics$flexible_runtime_conditions))
curves <- do.call(rbind, lapply(results, `[[`, "curve"))
vpc <- do.call(rbind, lapply(results, `[[`, "vpc"))
metrics <- do.call(rbind, lapply(results, `[[`, "metrics"))
estimates <- do.call(rbind, lapply(results, `[[`, "estimates"))

interval <- function(x) c(lower = unname(quantile(x, .025)),
  median = median(x), upper = unname(quantile(x, .975)))
summarize_curve <- function(value) do.call(rbind, lapply(split(value,
  interaction(value$arm, value$probability, drop = TRUE)), function(x) {
  limits <- interval(x$medianCL)
  data.frame(arm = x$arm[1L], probability = x$probability[1L],
    CRP = x$CRP[1L], t(limits))
}))
curveSummary <- summarize_curve(curves)
curve100 <- summarize_curve(subset(curves, replicate <= 100L))
curveStability <- merge(curve100, curveSummary,
  by = c("arm", "probability", "CRP"), suffixes = c("_100", "_all"))
for (name in c("lower", "median", "upper"))
  curveStability[[paste0("change_", name)]] <-
    curveStability[[paste0(name, "_all")]] -
    curveStability[[paste0(name, "_100")]]
stabilitySummary <- do.call(rbind, lapply(split(curveStability,
  curveStability$arm), function(x) data.frame(arm = x$arm[1L],
    max_abs_median_change = max(abs(x$change_median)),
    max_abs_lower_change = max(abs(x$change_lower)),
    max_abs_upper_change = max(abs(x$change_upper)))))
vpcSummary <- do.call(rbind, lapply(split(vpc,
  interaction(vpc$arm, vpc$quantile, vpc$time, drop = TRUE)), function(x) {
    value <- interval(x$value)
    data.frame(arm = x$arm[1L], quantile = x$quantile[1L],
      time = x$time[1L], t(value))
  }))
metricSummary <- do.call(rbind, lapply(split(metrics, metrics$arm), function(x) {
  value <- interval(x$log_rmse)
  data.frame(arm = x$arm[1L], n = nrow(x), t(value),
    mean = mean(x$log_rmse))
}))
paired <- reshape(metrics, idvar = "replicate", timevar = "arm",
  direction = "wide")
paired$reduction_percent <- 100 * (1 - paired$log_rmse.flexible /
  paired$log_rmse.gaussian)
paired$difference <- paired$log_rmse.gaussian - paired$log_rmse.flexible
halfWidth <- qt(.975, nrow(paired) - 1L) * sd(paired$difference) /
  sqrt(nrow(paired))
pairedSummary <- data.frame(n = nrow(paired),
  flexible_better = sum(paired$log_rmse.flexible < paired$log_rmse.gaussian),
  mean_paired_difference = mean(paired$difference),
  lower_mean_difference = mean(paired$difference) - halfWidth,
  upper_mean_difference = mean(paired$difference) + halfWidth,
  median_reduction_percent = median(paired$reduction_percent),
  lower_reduction_percent = unname(quantile(paired$reduction_percent, .025)),
  upper_reduction_percent = unname(quantile(paired$reduction_percent, .975)))

write.csv(curves, file.path(root, "curves_all.csv"), row.names = FALSE)
write.csv(curveSummary, file.path(root, "curves_summary.csv"), row.names = FALSE)
write.csv(curveStability, file.path(root, "curve_stability_100_vs_all.csv"),
  row.names = FALSE)
write.csv(stabilitySummary, file.path(root,
  "curve_stability_summary_100_vs_all.csv"), row.names = FALSE)
write.csv(vpc, file.path(root, "vpc_all.csv"), row.names = FALSE)
write.csv(vpcSummary, file.path(root, "vpc_summary.csv"), row.names = FALSE)
write.csv(metrics, file.path(root, "metrics_all.csv"), row.names = FALSE)
write.csv(metricSummary, file.path(root, "metrics_summary.csv"), row.names = FALSE)
write.csv(paired, file.path(root, "metrics_paired.csv"), row.names = FALSE)
write.csv(pairedSummary, file.path(root, "paired_summary.csv"), row.names = FALSE)
write.csv(estimates, file.path(root, "estimates_all.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(root, "diagnostics.csv"), row.names = FALSE)

labels <- c(truth = "Generating model", gaussian = "Gaussian FREM",
  flexible = "Flexible-margin FREM")
colours <- c(truth = "#222222", gaussian = "#2878B5", flexible = "#E07A24")
pCurve <- ggplot(subset(curveSummary, arm != "truth"), aes(CRP, median,
    colour = arm, fill = arm)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .16, colour = NA) +
  geom_line(linewidth = .85) +
  geom_line(data = subset(curveSummary, arm == "truth"),
    aes(CRP, median), inherit.aes = FALSE, colour = "#222222", linewidth = .8) +
  scale_colour_manual(values = colours, labels = labels, name = NULL) +
  scale_fill_manual(values = colours, labels = labels, name = NULL) +
  labs(x = "CRP", y = "Conditional median clearance",
    title = "A  Conditional FFEM relationship") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 10))

metrics$arm <- factor(metrics$arm, levels = c("gaussian", "flexible"),
  labels = c("Gaussian FREM", "Flexible-margin FREM"))
pMetric <- ggplot(metrics, aes(arm, log_rmse, fill = arm)) +
  geom_boxplot(width = .5, outlier.shape = NA, alpha = .65) +
  geom_line(aes(group = replicate), colour = "grey45", alpha = .28,
    linewidth = .3) +
  geom_point(aes(group = replicate), position = position_jitter(width = .07),
    alpha = .45, size = 1) +
  scale_fill_manual(values = c("Gaussian FREM" = "#2878B5",
    "Flexible-margin FREM" = "#E07A24"), guide = "none") +
  labs(x = NULL, y = "Conditional VPC log-RMSE",
    title = "B  Predictive error across fits") +
  theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 10))

png(file.path(root, "replicated_conditional_figure.png"), width = 2400,
  height = 1050, res = 300)
gridExtra::grid.arrange(pCurve, pMetric, nrow = 1L, widths = c(1.45, 1))
dev.off()
print(metricSummary)
print(pairedSummary)
print(stabilitySummary)
