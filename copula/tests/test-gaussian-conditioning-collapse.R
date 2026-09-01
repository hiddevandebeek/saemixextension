## Exact analytical collapse of a Gaussian R-vine with arbitrary continuous
## conditioning margins to an ordinary Gaussian random-effects model.
source("copula/tests/helper-load.R")

n_fail <- 0L
ok <- function(label, pass, detail = "") {
  if (!isTRUE(pass)) n_fail <<- n_fail + 1L
  cat(sprintf("%-72s %s %s\n", label,
              if (isTRUE(pass)) "PASS" else "**FAIL**", detail))
}

dmvn <- function(x, mean, covariance) {
  x <- as.matrix(x)
  centred <- if (is.matrix(mean) && identical(dim(mean), dim(x))) x - mean else
    sweep(x, 2, as.numeric(mean), "-")
  U <- chol(covariance)
  standardized <- forwardsolve(t(U), t(centred))
  -.5 * (ncol(x) * log(2 * pi) + 2 * sum(log(diag(U))) +
           colSums(standardized^2))
}

set.seed(82441)
d_eta <- 2L
d_cov <- 2L
d <- d_eta + d_cov
A0 <- matrix(rnorm(d * d), d, d)
R <- cov2cor(crossprod(A0) + diag(d))
vine <- copulaGaussianRvineFromCor(
  R, rvinecopulib::cvine_structure(c(3, 1, 4, 2)))
margins <- list(
  copulaMarginNormal(.32), copulaMarginNormal(.47),
  copulaMarginCovariateLognormal(log(72), .22),
  copulaMarginCovariateWeibull(3.1, 96))

conditioning <- cbind(
  WT = rlnorm(40, log(72), .22),
  eGFR = rweibull(40, 3.1, 96))
u_cov <- copulaMarginsEvaluate(conditioning, margins[3:4],
  numericalPolicy = "exact")$u
z_cov <- qnorm(u_cov)

ee <- seq_len(d_eta)
cc <- d_eta + seq_len(d_cov)
D_eta <- diag(copulaMarginScales(margins[ee]), d_eta)
regression <- D_eta %*% R[ee, cc, drop = FALSE] %*%
  solve(R[cc, cc, drop = FALSE])
omega <- D_eta %*% (R[ee, ee, drop = FALSE] -
  R[ee, cc, drop = FALSE] %*% solve(R[cc, cc, drop = FALSE]) %*%
    R[cc, ee, drop = FALSE]) %*% D_eta

subject <- rep(seq_len(nrow(conditioning)), each = 20L)
eta <- matrix(rnorm(length(subject) * d_eta), ncol = d_eta) %*% chol(omega)
eta <- eta + z_cov[subject, , drop = FALSE] %*% t(regression)
joint <- cbind(eta, conditioning[subject, , drop = FALSE])
log_joint <- copulaLogPrior(joint, vine, margins = margins,
  numericalPolicy = "exact")

log_cov_margins <- copulaMarginsEvaluate(
  conditioning[subject, , drop = FALSE], margins[cc],
  numericalPolicy = "exact")$log_margin
log_cov <- dmvn(z_cov[subject, , drop = FALSE], c(0, 0),
                R[cc, cc, drop = FALSE]) -
  rowSums(dnorm(z_cov[subject, , drop = FALSE], log = TRUE)) +
  log_cov_margins
log_conditional <- vapply(seq_along(subject), function(i) {
  dmvn(eta[i, , drop = FALSE], regression %*% z_cov[subject[i], ], omega)
}, numeric(1))

collapse_error <- max(abs((log_joint - log_cov) - log_conditional))
ok("GC1 joint Gaussian vine minus covariate density equals collapsed prior",
   collapse_error < 2e-9, sprintf("maxabs=%.3g", collapse_error))

set.seed(82442)
one <- 7L
draws <- matrix(rnorm(200000L * d_eta), ncol = d_eta) %*% chol(omega)
draws <- sweep(draws, 2, regression %*% z_cov[one, ], "+")
mean_error <- max(abs(colMeans(draws) - regression %*% z_cov[one, ]))
cov_error <- max(abs(cov(draws) - omega))
ok("GC2 collapsed conditional sampler has the analytic mean",
   mean_error < 3e-3, sprintf("maxabs=%.3g", mean_error))
ok("GC3 collapsed conditional sampler has the Schur covariance",
   cov_error < 1.5e-3, sprintf("maxabs=%.3g", cov_error))

## The deeper collapse does not require Gaussian eta margins.  On latent
## normal-score coordinates every continuous marginal is just a monotone
## observation map eta_j = F_j^{-1}{Phi(z_j)}.  The Gaussian copula therefore
## remains an ordinary conditional Gaussian latent-variable model.
latent_margins <- list(
  copulaMarginStudent(.32, 5), copulaMarginLaplace(.47),
  margins[[3]], margins[[4]])
latent_cov <- R[ee, ee, drop = FALSE] -
  R[ee, cc, drop = FALSE] %*% solve(R[cc, cc, drop = FALSE]) %*%
    R[cc, ee, drop = FALSE]
latent_regression <- R[ee, cc, drop = FALSE] %*%
  solve(R[cc, cc, drop = FALSE])
set.seed(82443)
z_eta <- matrix(rnorm(length(subject) * d_eta), ncol = d_eta) %*%
  chol(latent_cov)
z_eta <- z_eta + z_cov[subject, , drop = FALSE] %*%
  t(latent_regression)
eta_flexible <- copulaMarginsQuantile(pnorm(z_eta), latent_margins[ee])
joint_flexible <- cbind(eta_flexible,
                        conditioning[subject, , drop = FALSE])
log_joint_flexible <- copulaLogPrior(joint_flexible, vine,
  margins = latent_margins, numericalPolicy = "exact")
eta_eval <- copulaMarginsEvaluate(eta_flexible, latent_margins[ee],
  numericalPolicy = "exact")
log_conditional_flexible <- dmvn(
  z_eta,
  z_cov[subject, , drop = FALSE] %*% t(latent_regression),
  latent_cov) - rowSums(dnorm(z_eta, log = TRUE)) +
    eta_eval$log_margin
flexible_error <- max(abs((log_joint_flexible - log_cov) -
                          log_conditional_flexible))
ok("GC4 latent-normal collapse is exact for non-Gaussian eta margins",
   flexible_error < 2e-9, sprintf("maxabs=%.3g", flexible_error))

set.seed(82444)
z_draw <- matrix(rnorm(200000L * d_eta), ncol = d_eta) %*%
  chol(latent_cov)
z_draw <- sweep(z_draw, 2, latent_regression %*% z_cov[one, ], "+")
eta_draw <- copulaMarginsQuantile(pnorm(z_draw), latent_margins[ee])
roundtrip <- copulaMarginsEvaluate(eta_draw, latent_margins[ee],
  numericalPolicy = "exact")$u
transport_error <- max(abs(qnorm(roundtrip) - z_draw))
ok("GC5 flexible-margin sampler is the exact monotone transport",
   transport_error < 2e-8, sprintf("maxabs=%.3g", transport_error))

## Subject-specific missing covariates are integrated rather than imputed with
## a plug-in value. Compare each observed-data density with the corresponding
## independently constructed marginal MVN density on normal-score coordinates.
conditioning_missing <- conditioning[1:4, , drop = FALSE]
conditioning_missing[2, 2] <- NA_real_
conditioning_missing[3, 1] <- NA_real_
conditioning_missing[4, ] <- NA_real_
eta_missing <- eta_flexible[1:4, , drop = FALSE]
joint_missing <- cbind(eta_missing, conditioning_missing)
collapsed_missing <- copulaGaussianFremLogPrior(
  joint_missing, vine, latent_margins, d_eta, likelihoodTarget = "joint")
independent_missing <- vapply(seq_len(4L), function(i) {
  observed_cov <- which(!is.na(conditioning_missing[i, ]))
  index <- c(ee, d_eta + observed_cov)
  values <- joint_missing[i, index, drop = FALSE]
  evaluated <- copulaMarginsEvaluate(values, latent_margins[index],
                                     numericalPolicy = "exact")
  scores <- qnorm(evaluated$u)
  dmvn(scores, rep(0, length(index)), R[index, index, drop = FALSE]) -
    sum(dnorm(scores, log = TRUE)) + sum(evaluated$log_margin)
}, numeric(1))
missing_density_error <- max(abs(collapsed_missing - independent_missing))
ok("GC6 missing covariates are exactly marginalized pattern by pattern",
   missing_density_error < 2e-9,
   sprintf("maxabs=%.3g", missing_density_error))

## The independence proposal is the exact conditional prior for every pattern.
## Check on the latent normal-score scale so arbitrary eta margins do not obscure
## the Gaussian conditional moments.
set.seed(82445)
proposal_conditioning <- conditioning_missing[rep(seq_len(4L), each = 50000L), ]
proposal_eta <- copulaGaussianFremRandEta(
  proposal_conditioning, vine, latent_margins, d_eta)
proposal_score <- copulaMarginsEvaluate(
  proposal_eta, latent_margins[ee], numericalPolicy = "exact")$u
proposal_score <- qnorm(proposal_score)
conditional <- copulaGaussianFremConditional(
  conditioning_missing, vine, latent_margins, d_eta)
proposal_mean_error <- proposal_cov_error <- 0
for (i in seq_len(4L)) {
  rows <- ((i - 1L) * 50000L + 1L):(i * 50000L)
  proposal_mean_error <- max(proposal_mean_error,
    abs(colMeans(proposal_score[rows, , drop = FALSE]) - conditional$mean[i, ]))
  proposal_cov_error <- max(proposal_cov_error,
    abs(cov(proposal_score[rows, , drop = FALSE]) - conditional$covariance[[i]]))
}
ok("GC7 conditional proposals have the correct pattern-specific means",
   proposal_mean_error < 8e-3,
   sprintf("maxabs=%.3g", proposal_mean_error))
ok("GC8 conditional proposals have the correct pattern-specific covariances",
   proposal_cov_error < 8e-3,
   sprintf("maxabs=%.3g", proposal_cov_error))

conditioning_for_fit <- conditioning
conditioning_for_fit[c(2, 5, 11), 1] <- NA_real_
conditioning_for_fit[c(3, 7, 13), 2] <- NA_real_
auto_with_missing <- copulaFitCovariateMargins(conditioning_for_fit)
ok("GC9 automatic continuous-margin fitting uses all observed values",
   length(auto_with_missing) == 2L &&
     all(vapply(auto_with_missing, inherits, logical(1),
                "saemix_copula_margin")), "")

mstep_conditioning <- conditioning
mstep_conditioning[seq(2, nrow(conditioning), by = 5), 1] <- NA_real_
mstep_conditioning[seq(3, nrow(conditioning), by = 7), 2] <- NA_real_

old_state <- as.list(.cop)
copulaSet(vine, margins = latent_margins,
  conditioning = list(values = mstep_conditioning,
                      variableName = c("WT", "eGFR")),
  warmStartOnActivate = FALSE)
set.seed(82447)
api_proposal <- copulaRandEta(nrow(mstep_conditioning), mstep_conditioning)
ok("GC10 public population state accepts NA and exposes exact prior proposals",
   all(dim(api_proposal) == c(nrow(mstep_conditioning), d_eta)) &&
     all(is.finite(api_proposal)), "")
copulaRestoreState(old_state)

## Random-walk MH ratios may use p(eta | c_obs) instead of p(eta,c_obs).
## The omitted covariate density must cancel row by row, including under
## arbitrary eta margins and missing-covariate patterns.
eta_a <- eta_flexible[seq_len(4L), , drop = FALSE]
eta_b <- eta_a + matrix(c(.03, -.02, -.01, .04, .02, .01, -.03, -.02),
  nrow = 4L)
conditional_u <- copulaGaussianFremConditionalUetaEvaluator(
  conditioning_missing, vine, latent_margins, d_eta)
conditional_difference <- conditional_u(eta_b) - conditional_u(eta_a)
joint_difference <- -copulaGaussianFremLogPrior(
  cbind(eta_b, conditioning_missing), vine, latent_margins, d_eta, "joint") +
  copulaGaussianFremLogPrior(
    cbind(eta_a, conditioning_missing), vine, latent_margins, d_eta, "joint")
ok("GC11 cached conditional-prior evaluator preserves arbitrary-margin MH ratios",
   max(abs(conditional_difference - joint_difference)) < 2e-10,
   sprintf("maxabs=%.3g", max(abs(conditional_difference - joint_difference))))

normal_eta_margins <- latent_margins
normal_eta_margins[[1L]] <- copulaMarginNormal(.43)
normal_eta_margins[[2L]] <- copulaMarginNormal(.31)
normal_u <- copulaGaussianFremConditionalUetaEvaluator(
  conditioning_missing, vine, normal_eta_margins, d_eta)
normal_a <- matrix(rnorm(8, sd = .3), 4L, d_eta)
normal_b <- normal_a + matrix(rnorm(8, sd = .04), 4L, d_eta)
normal_difference <- normal_u(normal_b) - normal_u(normal_a)
normal_joint_difference <- -copulaGaussianFremLogPrior(
  cbind(normal_b, conditioning_missing), vine, normal_eta_margins,
  d_eta, "joint") + copulaGaussianFremLogPrior(
    cbind(normal_a, conditioning_missing), vine, normal_eta_margins,
    d_eta, "joint")
ok("GC12 Normal-eta conditional collapse preserves exact MH ratios",
   max(abs(normal_difference - normal_joint_difference)) < 2e-10,
   sprintf("maxabs=%.3g", max(abs(normal_difference - normal_joint_difference))))

rowwise_probe <- rbind(
  c(normal_a[1L, ], conditioning_missing[1L, ]),
  c(normal_a[1L, ], conditioning_missing[1L, ]))
rowwise_probe[2L, 4L] <- -1
rowwise_density <- copulaGaussianFremLogPrior(rowwise_probe,
  vine, normal_eta_margins, d_eta, "joint")
ok("GC13 invalid proposal handling is rowwise",
   is.finite(rowwise_density[1L]) && !is.finite(rowwise_density[2L]), "")

cat(sprintf("\n%d failure(s)\n", n_fail))
if (n_fail > 0L) quit(status = 1L)
