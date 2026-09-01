source("copula/tests/helper-load.R")

set.seed(290827L)
covariates <- cbind(
  WT = stats::rlnorm(40, log(70), .2),
  eGFR = stats::rnorm(40, 90, 22))
automatic <- gaussianCopulaFrem(
  etaSd = c(.22, .30), covariates = covariates)
stopifnot(length(automatic$arguments$margins) == 4L,
  all(vapply(automatic$arguments$margins[3:4], function(margin)
    identical(margin$metadata$selection$criterion, "AIC"), logical(1))))

population <- gaussianCopulaFrem(
  etaSd = c(.22, .30), covariates = covariates,
  covariateMargins = list(
    copulaMarginCovariateLognormal(log(70), .2),
    copulaMarginCovariateNormal(90, 22)))
stopifnot(inherits(population, "saemixPopulation"),
  identical(population$arguments$populationAlgorithm, "score-sa"))

copulaUsePopulation(population)
state <- copulaGet()
stopifnot(state$dEta == 2L, state$dConditioning == 2L,
  copulaIsFullGaussianVine(state$vine, 4L),
  identical(vapply(state$margins, `[[`, character(1), "name"),
    c("normal", "normal", "lognormal", "normal")))
copulaClear()

categorical <- cbind(SEX = rep(0:1, 20L),
  STAGE = rep(0:2, length.out = 40L))
categorical[c(4L, 17L), 1L] <- NA_real_
categorical[c(8L, 17L), 2L] <- NA_real_
multiCategorical <- gaussianCopulaFrem(
  etaSd = c(.22, .30), covariates = categorical,
  covariateMargins = list(copulaMarginBernoulli(.5),
    copulaMarginOrdinal(c(.3, .4, .3), labels = 0:2)))
copulaUsePopulation(multiCategorical)
multiState <- copulaGet()
stopifnot(multiState$dConditioning == 2L,
  sum(vapply(multiState$margins, function(m)
    identical(m$type, "discrete"), logical(1))) == 2L,
  identical(multiState$populationAlgorithm, "score-sa"))
copulaClear()

cat("minimal Gaussian-copula FREM API checks passed\n")
