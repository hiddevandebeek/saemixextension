## Proper binary/ordinal Gaussian-copula FREM likelihood and conditional kernel.
suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})
nfail <- 0L
ok <- function(label, pass, detail = "") {
  if (!isTRUE(pass)) nfail <<- nfail + 1L
  cat(sprintf("%-76s %s %s\n", label,
    if (isTRUE(pass)) "PASS" else "**FAIL**", detail))
}

rho <- .58; R <- matrix(c(1, rho, rho, 1), 2L)
structure <- rvinecopulib::cvine_structure(c(2, 1))
continuousVine <- copulaGaussianRvineFromCor(R, structure)
binary <- copulaMarginBernoulli(.3)
margins <- list(copulaMarginNormal(.4), binary)
vine <- copulaVineForMargins(continuousVine, margins)
eta <- matrix(seq(-.7, .7, length.out = 31L), ncol = 1L)
category <- rep(c(0, 1), length.out = nrow(eta))
joint <- cbind(eta, category); z <- eta[, 1L] / .4; cut <- qnorm(.7)
p0 <- pnorm((cut - rho * z) / sqrt(1 - rho^2))
oracle <- dnorm(eta[, 1L], sd = .4, log = TRUE) +
  ifelse(category == 0, log(p0), log1p(-p0))
actual <- copulaGaussianFremLogPrior(joint, vine, margins, 1L, "joint")
ok("CAT1 Bernoulli/Normal Gaussian-copula likelihood equals probit rectangle",
  max(abs(actual - oracle)) < 2e-11,
  sprintf("max error=%.3g", max(abs(actual - oracle))))

ordinal <- copulaMarginOrdinal(c(.2, .5, .3), labels = 0:2)
ordinalMargins <- list(copulaMarginNormal(.4), ordinal)
ordinalVine <- copulaVineForMargins(continuousVine, ordinalMargins)
ordinalCategory <- rep(0:2, length.out = nrow(eta))
upper <- qnorm(c(.2, .7, 1))[ordinalCategory + 1L]
lower <- qnorm(c(0, .2, .7))[ordinalCategory + 1L]
ordinalProbability <- pnorm((upper - rho * z) / sqrt(1 - rho^2)) -
  pnorm((lower - rho * z) / sqrt(1 - rho^2))
ordinalOracle <- dnorm(eta[, 1L], sd = .4, log = TRUE) + log(ordinalProbability)
ordinalActual <- copulaGaussianFremLogPrior(cbind(eta, ordinalCategory),
  ordinalVine, ordinalMargins, 1L, "joint")
ok("CAT2 three-level ordinal likelihood equals its latent-normal intervals",
  max(abs(ordinalActual - ordinalOracle)) < 2e-11,
  sprintf("max error=%.3g", max(abs(ordinalActual - ordinalOracle))))

## Missing category integrates to one and leaves the eta marginal unchanged.
missingActual <- copulaGaussianFremLogPrior(cbind(eta, NA_real_),
  vine, margins, 1L, "joint")
ok("CAT3 a missing categorical covariate is marginalized exactly",
  max(abs(missingActual - dnorm(eta[, 1L], sd = .4, log = TRUE))) < 2e-11)

## Conditional and joint prior ratios are identical for categorical events.
conditioning <- matrix(rep(c(0, 1), each = 8L), ncol = 1L)
eta0 <- matrix(seq(-.4, .35, length.out = 16L), ncol = 1L)
eta1 <- eta0 + .07
kernel <- copulaGaussianFremConditionalKernel(conditioning, vine, margins, 1L)
conditionalDifference <- kernel$negative(eta1) - kernel$negative(eta0)
jointDifference <- -copulaGaussianFremLogPrior(cbind(eta1, conditioning),
  vine, margins, 1L, "joint") +
  copulaGaussianFremLogPrior(cbind(eta0, conditioning),
    vine, margins, 1L, "joint")
ok("CAT4 categorical conditional kernel preserves exact MH ratios",
  max(abs(conditionalDifference - jointDifference)) < 2e-11,
  sprintf("max error=%.3g", max(abs(conditionalDifference - jointDifference))))

## Exact truncated-normal conditional sampler moments for each binary category.
set.seed(280831L)
n <- 60000L; cond <- matrix(rep(c(0, 1), each = n), ncol = 1L)
sampler <- copulaGaussianFremConditionalKernel(cond, vine, margins, 1L)
draw <- sampler$random()[, 1L] / .4
truncatedMoments <- function(a, b) {
  pa <- pnorm(a); pb <- pnorm(b); mass <- pb - pa
  mean <- (dnorm(a) - dnorm(b)) / mass
  ad <- if (is.finite(a)) a * dnorm(a) else 0
  bd <- if (is.finite(b)) b * dnorm(b) else 0
  variance <- 1 + (ad - bd) / mass - mean^2
  c(mean = mean, variance = variance)
}
m0 <- truncatedMoments(-Inf, cut); m1 <- truncatedMoments(cut, Inf)
expectedMean <- c(rho * m0["mean"], rho * m1["mean"])
expectedVariance <- c(1 - rho^2 + rho^2 * m0["variance"],
  1 - rho^2 + rho^2 * m1["variance"])
sampleMean <- c(mean(draw[seq_len(n)]), mean(draw[n + seq_len(n)]))
sampleVariance <- c(var(draw[seq_len(n)]), var(draw[n + seq_len(n)]))
ok("CAT5 categorical conditional sampler has selection-normal moments",
  max(abs(sampleMean - expectedMean)) < .012 &&
    max(abs(sampleVariance - expectedVariance)) < .015,
  sprintf("mean error=%.3g variance error=%.3g",
    max(abs(sampleMean - expectedMean)),
    max(abs(sampleVariance - expectedVariance))))

## The public state accepts a discrete conditioning margin and samples eta.
oldState <- as.list(.cop)
copulaSet(vine, margins = margins,
  conditioning = list(values = matrix(c(0, 1, NA, 0), ncol = 1L),
    variableName = "SEX"), warmStartOnActivate = FALSE)
set.seed(280832L); publicDraw <- copulaRandEta(4L, matrix(c(0, 1, NA, 0), ncol = 1L))
ok("CAT6 public Gaussian-copula FREM supports observed/missing binary covariates",
  all(dim(publicDraw) == c(4L, 1L)) && all(is.finite(publicDraw)))
copulaRestoreState(oldState)

## Literal common-Q fitting jointly updates the Bernoulli mass and Gaussian
## dependence; no continuous-density shortcut is used for the category.
set.seed(280833L)
nfit <- 650L
latent <- matrix(rnorm(2L * nfit), ncol = 2L) %*% chol(R)
fitData <- cbind(.4 * latent[, 1L], as.numeric(latent[, 2L] > cut))
fitMargins0 <- list(copulaMarginNormal(.34), copulaMarginBernoulli(.5))
fitVine0 <- copulaVineForMargins(
  copulaGaussianRvineFromCor(matrix(c(1, .08, .08, 1), 2L), structure),
  fitMargins0)
weights <- rep(1 / nfit, nfit)
q0 <- sum(weights * copulaGaussianFremLogPrior(
  fitData, fitVine0, fitMargins0, 1L, "joint"))
fitted <- copulaMaximiseGaussianFrem(fitData, weights, fitMargins0,
  fitVine0, 2L, 1L, maxit = 160L, diagnostics = TRUE)
q1 <- sum(weights * copulaGaussianFremLogPrior(
  sweep(fitData, 2L, fitted$delta, "-"), fitted$vine,
  fitted$margins, 1L, "joint"))
fittedProbability <- copulaCategoricalProbabilities(
  fitted$margins[[2L]]$parameters)[2L]
fittedRho <- copulaGaussianRvineCor(fitted$vine)[1L, 2L]
ok("CAT7 Bernoulli probability and dependence maximize the literal common Q",
  q1 + 1e-10 >= q0 && abs(fittedProbability - .3) < .06 &&
    abs(fittedRho - rho) < .10,
  sprintf("Q %.5f -> %.5f prob=%.3f rho=%.3f",
    q0, q1, fittedProbability, fittedRho))

badNominal <- try(copulaMarginCategorical(c(.2, .5, .3),
  labels = c(10, 20, 30)), silent = TRUE)
nominal <- copulaMarginCategorical(c(.2, .5, .3),
  labels = c(10, 20, 30), latent_order = c(30, 10, 20))
nominalMargins <- list(copulaMarginNormal(.4), nominal)
nominalVine <- copulaVineForMargins(continuousVine, nominalMargins)
nominalMass <- vapply(c(10, 20, 30), function(label)
  exp(copulaGaussianFremLogPrior(cbind(eta, label), nominalVine,
    nominalMargins, 1L, "joint")), numeric(nrow(eta)))
nominalSum <- rowSums(nominalMass)
ok("CAT8 nominal margins require an explicit order and remain normalized",
  inherits(badNominal, "try-error") &&
    max(abs(nominalSum - dnorm(eta[, 1L], sd = .4))) < 2e-11,
  sprintf("max normalization error=%.3g",
    max(abs(nominalSum - dnorm(eta[, 1L], sd = .4)))))

Rmulti <- matrix(c(1, .45, -.25, .45, 1, .30, -.25, .30, 1), 3L)
multiStructure <- rvinecopulib::cvine_structure(c(2, 1, 3))
multiMargins <- list(copulaMarginNormal(.4),
  copulaMarginBernoulli(.3), copulaMarginBernoulli(.6))
multiVine <- copulaVineForMargins(
  copulaGaussianRvineFromCor(Rmulti, multiStructure), multiMargins)
combinations <- expand.grid(c1 = 0:1, c2 = 0:1)
multiMass <- vapply(seq_len(nrow(combinations)), function(k)
  exp(copulaGaussianFremLogPrior(cbind(eta,
    combinations$c1[k], combinations$c2[k]), multiVine,
    multiMargins, 1L, "joint")), numeric(nrow(eta)))
multiError <- max(abs(rowSums(multiMass) - dnorm(eta[, 1L], sd = .4)))
ok("CAT9 multiple categorical covariates use normalized MVN rectangles",
  multiError < 2e-8, sprintf("max normalization error=%.3g", multiError))

## The theorem-backed score uses fixed-support latent uniforms instead of
## differentiating the numerical rectangle probability.
set.seed(280837L)
multiConditioning <- cbind(
  c(0, 1, NA, 0, 1, NA), c(1, 0, 1, NA, 1, 0))
multiEta <- matrix(seq(-.45, .45, length.out = 6L), ncol = 1L)
augmented <- copulaGaussianFremAugmentMixedConditioning(
  multiEta, multiConditioning, multiVine, multiMargins, 1L)
completed <- augmented$conditioning
uniform <- augmented$categoricalUniform
observedPreserved <- all(completed[!is.na(multiConditioning)] ==
  multiConditioning[!is.na(multiConditioning)])
ok("CAT10 multivariate categorical augmentation completes missing values on fixed support",
  observedPreserved && !anyNA(completed) &&
    all(uniform[, 2:3] > 0 & uniform[, 2:3] < 1),
  sprintf("uniform range=[%.4g, %.4g]",
    min(uniform[, 2:3]), max(uniform[, 2:3])))

augmentedRows <- cbind(multiEta, completed)
augmentedEvaluation <- copulaGaussianFremAugmentedEvaluateMargins(
  augmentedRows, multiMargins, uniform)
ok("CAT11 fixed-support augmented mixed density is finite",
  all(augmentedEvaluation$valid) &&
    all(is.finite(augmentedEvaluation$z)) &&
    all(is.finite(augmentedEvaluation$logMargin)))

scoreStep <- copulaGaussianFremPopulationScoreStep(
  augmentedRows, rep(1 / nrow(augmentedRows), nrow(augmentedRows)),
  multiMargins, multiVine, 3L, 1L, gain = .02,
  categoricalUniform = uniform, scoreScale = .005)
ok("CAT12 multiple categorical coordinates use the score recursion without rectangle derivatives",
  identical(scoreStep$scoreMethod, "global-centered-difference") &&
    all(is.finite(scoreStep$score)) &&
    grepl("fixed-support", scoreStep$scoreTheory$categoricalAugmentation),
  sprintf("max score=%.3g", scoreStep$scoreMax))

cat(sprintf("\n%d failure(s)\n", nfail))
if (nfail) quit(status = 1L)
