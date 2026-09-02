suppressPackageStartupMessages({
  library(devtools); devtools::load_all("C:/package/saemix-copula", quiet = TRUE)
  library(rvinecopulib)
})
set.seed(902301)
n <- 180L; samples <- 300L; rows <- n * samples
R <- matrix(c(1, .2, .15, .2, 1, .6, .15, .6, 1), 3L, 3L)
vine <- copulaGaussianRvineFromCor(R, cvine_structure(c(3, 1, 2)))
incumbent <- list(vine = vine, d = 3L, dEta = 2L, dConditioning = 1L,
  margins = list(copulaMarginNormal(.22), copulaMarginNormal(.4),
    copulaMarginCovariateGamma(2, 20)))
predictor <- matrix(rep(log(c(20, 3.5)), each = rows), nrow = rows)
eta <- cbind(rnorm(rows, 0, .22), rnorm(rows, 0, .4))
conditioning <- matrix(rgamma(rows, 2, scale = 20), ncol = 1L)
transform <- c(1L, 1L)
data <- list(n = n, samples = samples, eta = eta, predictor = predictor,
  typical = copulaWorkingToNatural(predictor, transform),
  natural = copulaWorkingToNatural(predictor + eta, transform),
  transform = transform, conditioning = conditioning)
candidate <- list(copulaNaturalMarginLognormal(.22),
  copulaNaturalMarginGamma(2.5))
genericMargins <- c(candidate, incumbent$margins[3L])
generic <- copulaNaturalFremLogPrior(eta, conditioning, vine,
  genericMargins, 2L, predictor, transform, "joint")
prepared <- copulaNaturalBridgePrepare(data, incumbent)
cached <- copulaNaturalBridgeCandidateLogPrior(candidate, prepared, incumbent)
stopifnot(max(abs(generic - cached)) < 1e-10)

oldTime <- system.time(for (i in 1:10) copulaNaturalFremLogPrior(eta,
  conditioning, vine, genericMargins, 2L, predictor, transform, "joint"))["elapsed"]
newTime <- system.time(for (i in 1:10)
  copulaNaturalBridgeCandidateLogPrior(candidate, prepared, incumbent))["elapsed"]
stopifnot(newTime < oldTime)
cat("screen cache exactness passed; speedup", round(oldTime / newTime, 2), "x\n")
