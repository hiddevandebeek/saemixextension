## Internal test of eta-family screening from full incumbent-posterior draws.
## The observation model is y = eta + epsilon, so direct marginal likelihoods
## are available for independent verification.

log_add_exp <- function(a, b) {
  maximum <- pmax(a, b)
  maximum + log(exp(a - maximum) + exp(b - maximum))
}

row_log_mean_exp <- function(value) {
  maximum <- apply(value, 1L, max)
  maximum + log(rowMeans(exp(value - maximum)))
}

log_laplace <- function(x, sd) {
  scale <- sd / sqrt(2)
  -log(2 * scale) - abs(x) / scale
}

log_normal_laplace <- function(y, eta_sd, residual_sd) {
  scale <- eta_sd / sqrt(2)
  common <- -log(2 * scale) + residual_sd^2 / (2 * scale^2)
  left <- -y / scale + stats::pnorm(
    y / residual_sd - residual_sd / scale, log.p = TRUE)
  right <- y / scale + stats::pnorm(
    -y / residual_sd - residual_sd / scale, log.p = TRUE)
  common + log_add_exp(left, right)
}

direct_fit <- function(y, residual_sd, family) {
  objective <- function(log_sd) {
    sd <- exp(log_sd)
    value <- if (family == "normal")
      stats::dnorm(y, 0, sqrt(sd^2 + residual_sd^2), log = TRUE) else
      log_normal_laplace(y, sd, residual_sd)
    -sum(value)
  }
  fit <- stats::optimize(function(x) objective(x), log(c(.08, 1.2)))
  c(sd = exp(fit$minimum), loglik = -fit$objective)
}

posterior_screen <- function(y, residual_sd, draws, incumbent_sd, family) {
  incumbent_loglik <- stats::dnorm(y, 0,
    sqrt(incumbent_sd^2 + residual_sd^2), log = TRUE)
  incumbent_logprior <- stats::dnorm(draws, 0, incumbent_sd, log = TRUE)
  evaluate <- function(log_sd, diagnostics = FALSE) {
    sd <- exp(log_sd)
    candidate <- if (family == "normal")
      stats::dnorm(draws, 0, sd, log = TRUE) else log_laplace(draws, sd)
    log_weight <- candidate - incumbent_logprior
    ratio <- row_log_mean_exp(log_weight)
    loglik <- sum(incumbent_loglik + ratio)
    if (!diagnostics) return(-loglik)
    shifted <- exp(log_weight - apply(log_weight, 1L, max))
    ess <- rowSums(shifted)^2 / rowSums(shifted^2)
    relative_variance <- apply(shifted, 1L, stats::var) /
      (ncol(shifted) * rowMeans(shifted)^2)
    list(sd = sd, loglik = loglik, ess = ess,
      mcse = sqrt(sum(relative_variance)))
  }
  fit <- stats::optimize(function(x) evaluate(x), log(c(.08, 1.2)))
  evaluate(fit$minimum, diagnostics = TRUE)
}

run_replicate <- function(generating_family, replicate, n = 100L,
                          posterior_draws = 2000L, eta_sd = .45,
                          residual_sd = .25) {
  set.seed(910000L + 1000L * match(generating_family,
    c("normal", "laplace")) + replicate)
  eta <- if (generating_family == "normal")
    stats::rnorm(n, 0, eta_sd) else {
    scale <- eta_sd / sqrt(2)
    sign(stats::runif(n) - .5) * stats::rexp(n, 1 / scale)
  }
  y <- eta + stats::rnorm(n, 0, residual_sd)
  incumbent_sd <- sqrt(max(mean(y^2) - residual_sd^2, .02^2))
  posterior_variance <- 1 / (1 / incumbent_sd^2 + 1 / residual_sd^2)
  posterior_mean <- posterior_variance * y / residual_sd^2
  draws <- matrix(stats::rnorm(n * posterior_draws,
    mean = rep(posterior_mean, posterior_draws),
    sd = sqrt(posterior_variance)), nrow = n, ncol = posterior_draws)

  rows <- lapply(c("normal", "laplace"), function(candidate) {
    direct <- direct_fit(y, residual_sd, candidate)
    posterior <- posterior_screen(y, residual_sd, draws, incumbent_sd,
      candidate)
    data.frame(generating_family = generating_family, replicate = replicate,
      candidate = candidate, incumbent_sd = incumbent_sd,
      direct_sd = direct[["sd"]], posterior_sd = posterior$sd,
      direct_loglik = direct[["loglik"]],
      posterior_loglik = posterior$loglik,
      posterior_minus_direct = posterior$loglik - direct[["loglik"]],
      posterior_mcse = posterior$mcse,
      median_ess_fraction = stats::median(posterior$ess) / posterior_draws,
      minimum_ess_fraction = min(posterior$ess) / posterior_draws)
  })
  do.call(rbind, rows)
}

results <- do.call(rbind, lapply(c("normal", "laplace"), function(family)
  do.call(rbind, lapply(seq_len(12L), function(replicate)
    run_replicate(family, replicate)))))
selection <- do.call(rbind, lapply(split(results,
  interaction(results$generating_family, results$replicate, drop = TRUE)),
  function(x) data.frame(generating_family = x$generating_family[1L],
    replicate = x$replicate[1L],
    selected_direct = x$candidate[which.max(x$direct_loglik)],
    selected_posterior = x$candidate[which.max(x$posterior_loglik)])))

summary <- do.call(rbind, lapply(split(results,
  interaction(results$generating_family, results$candidate, drop = TRUE)),
  function(x) data.frame(generating_family = x$generating_family[1L],
    candidate = x$candidate[1L],
    mean_loglik_error = mean(x$posterior_minus_direct),
    rmse_loglik_error = sqrt(mean(x$posterior_minus_direct^2)),
    mean_reported_mcse = mean(x$posterior_mcse),
    median_ess_fraction = stats::median(x$median_ess_fraction),
    minimum_ess_fraction = min(x$minimum_ess_fraction))))

selection_summary <- data.frame(
  generating_family = c("normal", "laplace"),
  direct_correct = vapply(c("normal", "laplace"), function(family)
    mean(selection$selected_direct[selection$generating_family == family] ==
      family), numeric(1)),
  posterior_correct = vapply(c("normal", "laplace"), function(family)
    mean(selection$selected_posterior[
      selection$generating_family == family] == family), numeric(1)),
  posterior_matches_direct = vapply(c("normal", "laplace"), function(family)
    mean(selection$selected_posterior[
      selection$generating_family == family] == selection$selected_direct[
        selection$generating_family == family]), numeric(1)))

output <- "copula/experiments/posterior-margin-selection/results"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(output, "informative_results.csv"),
  row.names = FALSE)
write.csv(selection, file.path(output, "informative_selection.csv"),
  row.names = FALSE)
write.csv(summary, file.path(output, "informative_summary.csv"),
  row.names = FALSE)
write.csv(selection_summary,
  file.path(output, "informative_selection_summary.csv"),
  row.names = FALSE)
print(summary)
print(selection_summary)

weak <- do.call(rbind, lapply(c("normal", "laplace"), function(family)
  do.call(rbind, lapply(seq_len(12L), function(replicate)
    run_replicate(family, replicate, residual_sd = .60)))))
weakSelection <- do.call(rbind, lapply(split(weak,
  interaction(weak$generating_family, weak$replicate, drop = TRUE)),
  function(x) data.frame(generating_family = x$generating_family[1L],
    replicate = x$replicate[1L],
    selected_direct = x$candidate[which.max(x$direct_loglik)],
    selected_posterior = x$candidate[which.max(x$posterior_loglik)])))
weakSelectionSummary <- data.frame(
  generating_family = c("normal", "laplace"),
  direct_correct = vapply(c("normal", "laplace"), function(family)
    mean(weakSelection$selected_direct[
      weakSelection$generating_family == family] == family), numeric(1)),
  posterior_correct = vapply(c("normal", "laplace"), function(family)
    mean(weakSelection$selected_posterior[
      weakSelection$generating_family == family] == family), numeric(1)),
  posterior_matches_direct = vapply(c("normal", "laplace"), function(family)
    mean(weakSelection$selected_posterior[
      weakSelection$generating_family == family] ==
      weakSelection$selected_direct[
        weakSelection$generating_family == family]), numeric(1)))
mismatch <- subset(weakSelection,
  selected_direct != selected_posterior & generating_family == "laplace")
weakRerun <- if (nrow(mismatch)) do.call(rbind, lapply(mismatch$replicate,
  function(replicate) run_replicate("laplace", replicate,
    posterior_draws = 20000L, residual_sd = .60))) else data.frame()
write.csv(weak, file.path(output, "weak_information_results.csv"),
  row.names = FALSE)
write.csv(weakSelection, file.path(output, "weak_information_selection.csv"),
  row.names = FALSE)
write.csv(weakSelectionSummary,
  file.path(output, "weak_information_selection_summary.csv"),
  row.names = FALSE)
write.csv(weakRerun,
  file.path(output, "weak_information_rerun_20000.csv"), row.names = FALSE)
print(weakSelectionSummary)
if (nrow(weakRerun)) print(weakRerun[, c("replicate", "candidate",
  "direct_loglik", "posterior_loglik", "posterior_mcse",
  "minimum_ess_fraction")])
