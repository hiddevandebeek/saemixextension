## Fisher-identity oracle for two dependent categorical population coordinates.
source("copula/tests/helper-load.R")

set.seed(300840L)
R <- matrix(c(1, .35, -.20, .35, 1, .30, -.20, .30, 1), 3L)
structure <- rvinecopulib::cvine_structure(c(2, 1, 3))
margins <- list(copulaMarginNormal(.4),
  copulaMarginBernoulli(.32),
  copulaMarginOrdinal(c(.25, .50, .25), labels = 0:2))
vine <- copulaVineForMargins(
  copulaGaussianRvineFromCor(R, structure), margins)
etaValue <- .17
categoryValue <- matrix(c(1, 2), nrow = 1L)

layout <- copulaGaussianFremScoreLayout(
  margins, vine, 3L, 1L, withMu = FALSE)
internal <- copulaScoreToInternal(layout$native, layout$lower, layout$upper)
observedObjective <- function(v) {
  candidate <- copulaGaussianFremScoreMaterialize(v, layout)
  copulaGaussianFremLogPrior(cbind(etaValue, categoryValue),
    candidate$vine, candidate$margins, 1L, "joint")
}
h <- 1e-5
observedScore <- vapply(seq_along(internal), function(j) {
  plus <- minus <- internal; plus[j] <- plus[j] + h; minus[j] <- minus[j] - h
  (observedObjective(plus) - observedObjective(minus)) / (2 * h)
}, numeric(1))

nBatch <- 24L; batchSize <- 1000L
batchScore <- matrix(NA_real_, nBatch, length(internal))
for (b in seq_len(nBatch)) {
  eta <- matrix(rep(etaValue, batchSize), ncol = 1L)
  conditioning <- categoryValue[rep(1L, batchSize), , drop = FALSE]
  augmented <- copulaGaussianFremAugmentMixedConditioning(
    eta, conditioning, vine, margins, 1L)
  answer <- copulaGaussianFremPopulationScoreStep(
    cbind(eta, augmented$conditioning), rep(1 / batchSize, batchSize),
    margins, vine, 3L, 1L, gain = .01, withMu = FALSE,
    finiteDifference = 1e-4,
    categoricalUniform = augmented$categoricalUniform)
  batchScore[b, ] <- answer$score
}
estimated <- colMeans(batchScore)
mcse <- apply(batchScore, 2L, stats::sd) / sqrt(nBatch)
error <- abs(estimated - observedScore)
stopifnot(all(is.finite(observedScore)), all(is.finite(estimated)),
  max(error / pmax(mcse, 1e-5)) < 4.5,
  max(error) < .025)

cat(sprintf(paste0("multi-categorical Fisher oracle passed; max error %.4g, ",
  "max MCSE %.4g\n"), max(error), max(mcse)))
