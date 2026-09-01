suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
})

root <- "C:/package/saemix-copula/copula/experiments/non-gaussian-eta-demonstration"
args <- commandArgs(trailingOnly = TRUE)
expected <- if (length(args)) as.integer(args[1L]) else NA_integer_
files <- sort(list.files(file.path(root, "replicates"),
  pattern = "^replicate_[0-9]{3}[.]rds$", full.names = TRUE))
if (!length(files) || (!is.na(expected) && length(files) != expected))
  stop("replicate result count does not match the requested study size")
results <- lapply(files, readRDS)
stopifnot(all(vapply(results, `[[`, integer(1), "schema") == 1L),
  !anyDuplicated(vapply(results, `[[`, integer(1), "replicate")),
  !anyDuplicated(vapply(results, `[[`, integer(1), "data_seed")))
nFits <- length(results)

bind <- function(field) {
  values <- lapply(results, function(result) {
    value <- result[[field]]; value$replicate <- result$replicate; value
  })
  columns <- unique(unlist(lapply(values, names)))
  values <- lapply(values, function(value) {
    for (name in setdiff(columns, names(value))) value[[name]] <- NA
    value[, columns, drop = FALSE]
  })
  do.call(rbind, values)
}
estimates <- bind("summary"); density <- bind("density")
tail <- bind("tail"); vpc <- bind("vpc")
renameArm <- function(value) {
  value[value == "Gaussian fit"] <- "Lognormal fit"
  value[value %in% c("Flexible-margin fit", "Gamma fit")] <- "Free-margin fit"
  value
}
estimates$arm <- renameArm(estimates$arm)
density$arm <- renameArm(density$arm)
tail$arm <- renameArm(tail$arm)
vpc$arm <- renameArm(vpc$arm)

interval <- function(x) c(lower = unname(quantile(x, .025)),
  median = median(x), upper = unname(quantile(x, .975)))
summarize_grid <- function(data, keys) do.call(rbind, lapply(split(data,
  interaction(data[, keys, drop = FALSE], drop = TRUE)), function(x)
    cbind(x[1L, keys, drop = FALSE], as.data.frame(t(interval(x$value))))))
densitySummary <- summarize_grid(density, c("arm", "parameter"))
tailSummary <- summarize_grid(tail, c("arm", "cutoff"))
vpcSummary <- summarize_grid(vpc, c("arm", "probability", "time"))

paired <- merge(subset(estimates, arm == "Lognormal fit"),
  subset(estimates, arm == "Free-margin fit"), by = "replicate",
  suffixes = c("_gaussian", "_flexible"))
paired$likelihood_gain <- paired$log_likelihood_flexible -
  paired$log_likelihood_gaussian
paired$likelihood_difference_mcse <- sqrt(paired$likelihood_mcse_gaussian^2 +
  paired$likelihood_mcse_flexible^2)
paired$runtime_ratio <- paired$runtime_seconds_flexible /
  paired$runtime_seconds_gaussian

truthVpc <- subset(vpc, arm == "Generating model",
  select = c("replicate", "time", "probability", "value"))
names(truthVpc)[4L] <- "truth"
vpcMetrics <- do.call(rbind, lapply(c("Lognormal fit", "Free-margin fit"),
  function(armName) do.call(rbind, lapply(split(vpc[vpc$arm == armName, ],
    vpc[vpc$arm == armName, "replicate"]), function(fitted) {
      joined <- merge(fitted, truthVpc,
        by = c("replicate", "time", "probability"))
      data.frame(replicate = fitted$replicate[1L], arm = armName,
        all_quantile_log_rmse = sqrt(mean((log(joined$value) -
          log(joined$truth))^2)),
        extreme_quantile_log_rmse = sqrt(mean((log(joined$value[
          joined$probability %in% c(.01, .99)]) - log(joined$truth[
            joined$probability %in% c(.01, .99)]))^2)))
    }))))
metricSummary <- do.call(rbind, lapply(split(vpcMetrics, vpcMetrics$arm),
  function(x) data.frame(arm = x$arm[1L],
    all_median = median(x$all_quantile_log_rmse),
    all_lower = quantile(x$all_quantile_log_rmse, .025),
    all_upper = quantile(x$all_quantile_log_rmse, .975),
    extreme_median = median(x$extreme_quantile_log_rmse),
    extreme_lower = quantile(x$extreme_quantile_log_rmse, .025),
    extreme_upper = quantile(x$extreme_quantile_log_rmse, .975))))
pairedMetrics <- merge(subset(vpcMetrics, arm == "Lognormal fit"),
  subset(vpcMetrics, arm == "Free-margin fit"), by = "replicate",
  suffixes = c("_standard", "_free"))
pairedSummary <- data.frame(n = nrow(paired),
  flexible_likelihood_better = sum(paired$likelihood_gain > 0),
  likelihood_gain_median = median(paired$likelihood_gain),
  likelihood_gain_lower = quantile(paired$likelihood_gain, .025),
  likelihood_gain_upper = quantile(paired$likelihood_gain, .975),
  likelihood_gain_over_2mcse = sum(paired$likelihood_gain >
    2 * paired$likelihood_difference_mcse),
  retained_standard_count = sum(as.logical(paired$retained_standard_flexible),
    na.rm = TRUE),
  gamma_selected_count = sum(grepl("/gamma$", paired$selected_families_flexible)),
  weibull_selected_count = sum(grepl("/weibull$",
    paired$selected_families_flexible)),
  lognormal_selected_count = sum(grepl("/lognormal$",
    paired$selected_families_flexible)),
  all_vpc_better_count = sum(pairedMetrics$all_quantile_log_rmse_free <
    pairedMetrics$all_quantile_log_rmse_standard),
  extreme_vpc_better_count = sum(
    pairedMetrics$extreme_quantile_log_rmse_free <
      pairedMetrics$extreme_quantile_log_rmse_standard),
  gamma_shape_median = median(paired$shape_flexible, na.rm = TRUE),
  gamma_sd_median = median(paired$sd_flexible),
  gaussian_tail_probability_median = median(paired$tail_probability_gaussian),
  flexible_tail_probability_median = median(paired$tail_probability_flexible),
  runtime_ratio_median = median(paired$runtime_ratio),
  selection_seconds_median = median(paired$selection_seconds_flexible))

colours <- c("Generating model" = "#222222", "Lognormal fit" = "#2878B5",
  "Free-margin fit" = "#E07A24")
fitDensity <- subset(densitySummary, arm != "Generating model")
truthDensity <- subset(densitySummary, arm == "Generating model")
pDensity <- ggplot(fitDensity, aes(parameter, median, colour = arm, fill = arm)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .16, colour = NA) +
  geom_line(linewidth = .9) +
  geom_line(data = truthDensity, aes(parameter, median), inherit.aes = FALSE,
    colour = colours[["Generating model"]], linewidth = .85) +
  scale_colour_manual(values = colours) + scale_fill_manual(values = colours) +
  labs(title = "A  Recovered individual clearance distribution",
    subtitle = paste("Median and 95% range across", nFits,
      "independently fitted datasets"),
    x = "Clearance", y = "Density", colour = NULL, fill = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom",
    panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

fitTail <- subset(tailSummary, arm != "Generating model")
truthTail <- subset(tailSummary, arm == "Generating model")
pTail <- ggplot(fitTail, aes(cutoff, median, colour = arm, fill = arm)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .16, colour = NA) +
  geom_line(linewidth = .9) +
  geom_line(data = truthTail, aes(cutoff, median), inherit.aes = FALSE,
    colour = colours[["Generating model"]], linewidth = .85) +
  scale_y_log10() + scale_colour_manual(values = colours) +
  scale_fill_manual(values = colours) +
  labs(title = "B  Upper tail of individual clearance",
    subtitle = bquote(P(CL > x)~"across"~.(nFits)~"fitted datasets"),
    x = expression("Threshold "*x), y = "Exceedance probability (log scale)",
    colour = NULL, fill = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom",
    panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

displayProbabilities <- c(.01, .50, .99)
fitVpc <- subset(vpcSummary, arm != "Generating model" &
  probability %in% displayProbabilities)
truthVpcSummary <- subset(vpcSummary, arm == "Generating model" &
  probability %in% displayProbabilities)
fitVpc$arm <- factor(fitVpc$arm,
  levels = c("Free-margin fit", "Lognormal fit"))
truthVpcFacet <- do.call(rbind, lapply(levels(fitVpc$arm), function(label) {
  value <- truthVpcSummary; value$arm <- factor(label, levels = levels(fitVpc$arm))
  value
}))
shape <- c(`0.01` = 17, `0.5` = 16, `0.99` = 15)
pVpc <- ggplot(fitVpc, aes(time, median, group = probability,
    colour = arm, fill = arm)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = .16, colour = NA) +
  geom_line(linetype = "dashed", linewidth = .9) +
  geom_line(data = truthVpcFacet, aes(time, median, group = probability),
    inherit.aes = FALSE, colour = "#222222", linewidth = .65) +
  geom_point(data = truthVpcFacet, aes(time, median,
    shape = factor(probability)), inherit.aes = FALSE, colour = "#222222",
    size = 2.1) +
  facet_wrap(~arm, nrow = 1L) + scale_y_log10() +
  scale_colour_manual(values = colours) + scale_fill_manual(values = colours) +
  scale_shape_manual(values = shape, labels = c("1st", "Median", "99th")) +
  labs(title = "C  Replicated population VPC",
    subtitle = "Fitted medians and 95% ranges; black symbols show generating quantiles",
    x = "Time", y = "Concentration (log scale)", colour = NULL, fill = NULL,
    shape = "Generating") +
  theme_bw(base_size = 10) + theme(legend.position = "bottom",
    panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))

output <- file.path(root, "out")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
figurePath <- file.path(output, paste0("non_gaussian_eta_", nFits,
  "_dataset.png"))
png(figurePath, width = 2600,
  height = 1900, res = 260)
grid.arrange(pDensity, pTail, pVpc, ncol = 2L,
  layout_matrix = rbind(c(1, 2), c(3, 3)), heights = c(1, 1.2))
dev.off()
write.csv(estimates, file.path(output, "replicated_fit_summary.csv"),
  row.names = FALSE)
write.csv(paired, file.path(output, "replicated_paired_summary.csv"),
  row.names = FALSE)
write.csv(pairedSummary, file.path(output, "replicated_overall_summary.csv"),
  row.names = FALSE)
write.csv(vpcMetrics, file.path(output, "replicated_vpc_metrics.csv"),
  row.names = FALSE)
write.csv(metricSummary, file.path(output, "replicated_vpc_metrics_summary.csv"),
  row.names = FALSE)
write.csv(vpcSummary, file.path(output, "replicated_vpc_summary.csv"),
  row.names = FALSE)
print(pairedSummary)
print(metricSummary)
