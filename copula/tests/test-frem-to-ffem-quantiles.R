source("copula/tests/helper-load.R")

R <- matrix(c(1, .20, .45, .10,
              .20, 1, .25, .55,
              .45, .25, 1, .15,
              .10, .55, .15, 1), 4L, 4L, byrow = TRUE)
vine <- copulaGaussianRvineFromCor(R,
  rvinecopulib::cvine_structure(c(3, 1, 4, 2)))
conditioning <- cbind(WT = seq(35, 100, length.out = 40L),
  eGFR = seq(45, 130, length.out = 40L))

normalMargins <- list(copulaMarginNormal(.2), copulaMarginNormal(.3),
  copulaMarginCovariateGamma(4, 14),
  copulaMarginCovariateNormal(90, 20))
copulaSet(vine, margins = normalMargins,
  conditioning = list(values = conditioning,
    variableName = c("WT", "eGFR")),
  warmStartOnActivate = FALSE, guard = FALSE)
normalState <- copulaSnapshot(copulaGet(),
  variableName = c("eta_V", "eta_CL"))
linear <- copulaFremToFfem(normalState, "WT")
newdata <- matrix(qgamma(.8, 4, scale = 14), ncol = 1L,
  dimnames = list(NULL, "WT"))
expectedMean <- as.numeric(copulaFfemLocation(linear, newdata))
set.seed(4901L)
draws <- copulaFfemSimulate(linear, newdata, n = 100000L)
stopifnot(max(abs(colMeans(draws) - expectedMean)) < .003,
  max(abs(stats::cov(draws) - linear$omega)) < .002)

generalMargins <- normalMargins
generalMargins[[1L]] <- copulaMarginCenteredGamma(4, .2)
copulaSet(vine, margins = generalMargins,
  conditioning = list(values = conditioning,
    variableName = c("WT", "eGFR")),
  warmStartOnActivate = FALSE, guard = FALSE)
generalState <- copulaSnapshot(copulaGet(),
  variableName = c("eta_V", "eta_CL"))
general <- copulaFremToFfem(generalState, "WT")
stopifnot(inherits(try(copulaFremToFfem(generalState, "WT", "linear"),
  silent = TRUE), "try-error"))
probabilities <- c(.05, .5, .95)
analytic <- copulaFfemQuantile(general, newdata, probabilities)
set.seed(4902L)
draws <- copulaFfemSimulate(general, newdata, n = 150000L)
empirical <- do.call(rbind, lapply(seq_len(ncol(draws)), function(j)
  data.frame(eta = colnames(draws)[j], probability = probabilities,
    value = as.numeric(quantile(draws[, j], probabilities, names = FALSE)))))
merged <- merge(analytic, empirical, by = c("eta", "probability"),
  suffixes = c("_analytic", "_empirical"))
stopifnot(max(abs(merged$value_analytic - merged$value_empirical)) < .006)

copulaClear()
cat("FREM-to-FFEM conditional quantile checks passed\n")
