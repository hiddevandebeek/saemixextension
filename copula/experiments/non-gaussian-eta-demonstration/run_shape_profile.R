suppressPackageStartupMessages(library(parallel))

repo <- "C:/package/saemix-copula"
root <- file.path(repo, "copula", "experiments",
  "non-gaussian-eta-demonstration")
output <- file.path(root, "shape_profile")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
shapeGrid <- c(1.5, 2, 2.5, 3, 4, 6, 10)
workers <- 4L

cluster <- makeCluster(workers, outfile = "")
on.exit(stopCluster(cluster), add = TRUE)
clusterExport(cluster, c("repo", "root", "output", "shapeGrid"))
clusterEvalQ(cluster, {
  suppressPackageStartupMessages({
    library(devtools)
    load_all(repo, quiet = TRUE)
  })
  source(file.path(root, "functions.R"))
  NULL
})
results <- parLapplyLB(cluster, seq_along(shapeGrid), function(index) {
  shape <- shapeGrid[index]
  path <- file.path(output, sprintf("shape_%04.1f.rds", shape))
  simulation <- ng_simulate(940101L)
  gamma <- copulaNaturalMarginGamma(shape)
  gamma$free[names(gamma$parameters) == "shape"] <- FALSE
  population <- gaussianCopulaFrem(
    parameterMargins = list(copulaNaturalMarginLognormal(.25), gamma),
    correlation = matrix(c(1, .10, .10, 1), 2L, 2L),
    scoreScale = 1.0, scoreBurn = 75L, gainScale = .18, gainPower = .80)
  elapsed <- system.time(fit <- saemix(ng_model(), ng_data(simulation$data),
    ng_control(975000L + index), population = population))["elapsed"]
  fit@options$nmc.is <- 4000L
  fit <- suppressWarnings(llisCopula.saemix(fit, defensive = .20,
    batch = 100L, seed = 975100L))
  state <- copulaGet(fit)
  stopifnot(isTRUE(state$lastJoint$scoreTheory$runtimeConditionsObserved),
    length(state$margins[[2L]]$free) == 1L,
    !state$margins[[2L]]$free[1L])
  answer <- data.frame(shape = shape,
    fitted_sd = sqrt(shape) * as.numeric(fit@results@fixed.effects)[2L] *
      exp(-digamma(shape)),
    V = as.numeric(fit@results@fixed.effects)[1L],
    CL = as.numeric(fit@results@fixed.effects)[2L],
    residual = ng_residual(fit), log_likelihood = fit@results@ll.is,
    likelihood_mcse = attr(fit,
      "saemix.copula.likelihood")$se_loglik_total,
    runtime_seconds = unname(elapsed))
  saveRDS(answer, path)
  answer
})
profile <- do.call(rbind, results)
profile$relative_log_likelihood <- profile$log_likelihood -
  max(profile$log_likelihood)
write.csv(profile, file.path(root, "out", "gamma_shape_profile.csv"),
  row.names = FALSE)

suppressPackageStartupMessages(library(ggplot2))
p <- ggplot(profile, aes(shape, relative_log_likelihood)) +
  geom_hline(yintercept = -1.92, linetype = "dashed", colour = "grey45") +
  geom_errorbar(aes(ymin = relative_log_likelihood - 1.96 * likelihood_mcse,
    ymax = relative_log_likelihood + 1.96 * likelihood_mcse), width = .12,
    colour = "#E07A24") +
  geom_line(colour = "#E07A24", linewidth = .9) +
  geom_point(colour = "#E07A24", size = 2.4) +
  geom_vline(xintercept = 2.5, colour = "#222222", linewidth = .7) +
  scale_x_log10(breaks = shapeGrid) +
  labs(title = "Profile likelihood for centred-Gamma eta shape",
    subtitle = "Other parameters re-estimated for 1,500 iterations at each fixed shape",
    x = "Fixed Gamma shape (log scale)",
    y = "Log likelihood relative to the maximum") +
  theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"))
ggsave(file.path(root, "out", "gamma_shape_profile.png"), p,
  width = 8.2, height = 5.2, dpi = 260)
print(profile)
