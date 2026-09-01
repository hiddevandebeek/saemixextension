root <- "C:/package/saemix-copula/copula/experiments/combined-natural-frem-study"
quality <- read.csv(file.path(root, "pilot-out", "quality_summary.csv"))
fits <- read.csv(file.path(root, "pilot-out", "fit_summary.csv"))
checks <- c(
  ten_complete_unique = quality$replicates == 10L && quality$unique_data_seeds == 10L,
  non_gaussian_parameter_selected = quality$gamma_parameter_selections >= 6L,
  non_gaussian_covariate_selected = sum(fits$arm == "Flexible FREM" &
    fits$familyCRP != "normal") >= 8L,
  finite_stable_estimates = all(is.finite(fits$V)) && all(is.finite(fits$CL)) &&
    all(is.finite(fits$residual)) && all(fits$V > 10 & fits$V < 30) &&
    all(fits$CL > 1.5 & fits$CL < 6) && all(fits$residual > .05 & fits$residual < .25),
  joint_likelihood_signal = quality$likelihood_better >= 8L &&
    quality$likelihood_over_2mcse >= 8L,
  conditional_likelihood_signal = quality$conditional_likelihood_better >= 8L,
  population_prediction_signal = quality$population_vpc_better >= 6L,
  conditional_prediction_signal = quality$conditional_vpc_better >= 8L &&
    quality$relationship_better >= 8L)
print(data.frame(check = names(checks), passed = unname(checks)))
if (!all(checks)) stop("combined-study pilot did not pass inspection")
cat("Combined-study pilot passed all prespecified inspection gates.\n")
