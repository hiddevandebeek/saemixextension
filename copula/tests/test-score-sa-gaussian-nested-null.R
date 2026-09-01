## The paper estimator must reduce to the same likelihood as Gaussian FREM
## when every margin is Normal and the copula is Gaussian.
suppressPackageStartupMessages({
  source("copula/tests/helper-load.R")
  library(mvtnorm)
})

set.seed(290826L)
d <- 4L
R <- matrix(c(1, .25, -.15, .35,
              .25, 1, .30, -.10,
              -.15, .30, 1, .20,
              .35, -.10, .20, 1), d, d, byrow = TRUE)
R <- stats::cov2cor(copulaGaussianFremEnsurePd(R, 1e-8))
sd <- c(.22, .31, 14, 25)
margins <- list(copulaMarginNormal(sd[1]), copulaMarginNormal(sd[2]),
  copulaMarginCovariateNormal(70, sd[3]),
  copulaMarginCovariateNormal(90, sd[4]))
vine <- copulaGaussianRvineFromCor(R,
  rvinecopulib::cvine_structure(c(3, 1, 4, 2)))
x <- sweep(rmvnorm(300, sigma = diag(sd) %*% R %*% diag(sd)), 2,
  c(0, 0, 70, 90), "+")
oracle <- dmvnorm(x, mean = c(0, 0, 70, 90),
  sigma = diag(sd) %*% R %*% diag(sd), log = TRUE)
exact <- copulaGaussianFremLogPrior(x, vine, margins, dEta = 2L,
  likelihoodTarget = "joint")
stopifnot(max(abs(exact - oracle)) < 2e-10)

pairPath <- file.path(copulaTestRepo, "copula", "study_frem_validation",
  "out", "score_pair_theory_final", "score_pair_summary.csv")
if (file.exists(pairPath)) {
  pair <- read.csv(pairPath, stringsAsFactors = FALSE)
  pair <- pair[pair$scenario == "normal_normal", ]
  combined <- sqrt(pair$score_mcse^2 + pair$common_mcse^2)
  stopifnot(nrow(pair) == 10L,
    all(abs(pair$difference) <= 2 * combined),
    abs(mean(pair$difference)) < .1)
}

cat("score-SA Gaussian nested-null checks passed\n")
