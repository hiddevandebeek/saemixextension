suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})

structure <- rvinecopulib::cvine_structure(c(2, 1))
vine <- copulaGaussianRvineFromCor(matrix(c(1, .35, .35, 1), 2L), structure)
conditioning <- matrix(seq(40, 80, length.out = 30L), ncol = 1L)

copulaSet(vine, margins = list(copulaMarginNormal(.25),
    copulaMarginCovariateNormal(mean(conditioning), stats::sd(conditioning))),
  conditioning = conditioning, guard = FALSE, warmStartOnActivate = FALSE)
state <- copulaGet()
stopifnot(identical(state$populationAlgorithmRequested, "score-sa"),
  identical(state$populationAlgorithm, "score-sa"),
  identical(state$populationAlgorithmReason, "paper likelihood-score estimator"))

copulaSet(vine, margins = list(copulaMarginNormal(.25),
    copulaMarginCovariateGamma(4, 14)),
  conditioning = conditioning, guard = FALSE, warmStartOnActivate = FALSE)
state <- copulaGet()
stopifnot(identical(state$populationAlgorithm, "score-sa"))

badGain <- try(copulaSet(vine, margins = list(copulaMarginNormal(.25),
    copulaMarginCovariateGamma(4, 14)),
  conditioning = conditioning, guard = FALSE, warmStartOnActivate = FALSE,
  scoreGainPower = .75), silent = TRUE)
stopifnot(inherits(badGain, "try-error"))

customEtaUndeclared <- copulaMarginDistribution("norm",
  c(mean = 0, sd = .25), c(0, 1e-6), c(0, 1e3),
  free = c(FALSE, TRUE), centered = TRUE,
  scale = function(par) par["sd"], scale_is_sd = TRUE,
  roles = c("location", "scale"))
badSupport <- try(copulaSet(vine,
  margins = list(customEtaUndeclared, copulaMarginCovariateGamma(4, 14)),
  conditioning = conditioning, guard = FALSE, warmStartOnActivate = FALSE),
  silent = TRUE)
stopifnot(inherits(badSupport, "try-error"))

customEtaFixed <- copulaMarginDistribution("norm",
  c(mean = 0, sd = .25), c(0, 1e-6), c(0, 1e3),
  free = c(FALSE, TRUE), centered = TRUE,
  scale = function(par) par["sd"], scale_is_sd = TRUE,
  roles = c("location", "scale"), support_fixed = TRUE)
copulaSet(vine,
  margins = list(customEtaFixed, copulaMarginCovariateGamma(4, 14)),
  conditioning = conditioning, guard = FALSE, warmStartOnActivate = FALSE)
stopifnot(isTRUE(copulaGet()$margins[[1L]]$metadata$parameter_independent_support))

copulaSet(vine, margins = list(copulaMarginNormal(.25),
    copulaMarginCovariateGamma(4, 14)),
  conditioning = conditioning, guard = FALSE, warmStartOnActivate = FALSE,
  populationAlgorithm = "common-q")
state <- copulaGet()
stopifnot(identical(state$populationAlgorithmRequested, "common-q"),
  identical(state$populationAlgorithm, "common-q"),
  identical(state$populationAlgorithmReason,
    "explicit compatibility/reference backend"))

vine3 <- copulaGaussianRvineFromCor(
  matrix(c(1, .2, .3, .2, 1, .1, .3, .1, 1), 3L),
  rvinecopulib::cvine_structure(c(3, 1, 2)))
copulaSet(vine3,
  margins = list(copulaMarginNormal(.25), copulaMarginBernoulli(.3),
    copulaMarginBernoulli(.6)),
  conditioning = cbind(C1 = rep(c(0, 1), 15L),
    C2 = rep(c(0, 0, 1), 10L)),
  guard = FALSE, warmStartOnActivate = FALSE,
  populationAlgorithm = "score-sa")
categoricalState <- copulaGet()
stopifnot(identical(categoricalState$populationAlgorithm, "score-sa"),
  categoricalState$dConditioning == 2L,
  all(vapply(categoricalState$margins[2:3], function(m)
    identical(m$type, "discrete"), logical(1))))

copulaClear()
cat("likelihood-score population-algorithm checks passed\n")
