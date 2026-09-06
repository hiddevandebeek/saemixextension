## Small end-to-end Gaussian-copula FREM fit with a skewed continuous
## covariate and subject-specific missingness. This exercises the actual SAEM
## E-step, exact conditional independence proposal, empirical-Q pool, and
## collapsed M-step rather than only the standalone density helpers.
suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})
setwd("C:/package/saemix-copula/copula")
source("R/simpleModels.R")

gR <- matrix(c(
  1, .35, .55,
  .35, 1, -.30,
  .55, -.30, 1), 3L, 3L, byrow = TRUE)
stopifnot(all(eigen(gR, symmetric = TRUE, only.values = TRUE)$values > 0))
vine <- copulaGaussianRvineFromCor(gR,
  rvinecopulib::cvine_structure(c(3, 1, 2)))
margins <- list(copulaMarginNormal(.28), copulaMarginNormal(.34),
                copulaMarginCovariateLognormal(log(72), .20))

set.seed(82601)
N <- 24L
u <- rvinecopulib::rvinecop(N, vine)
joint <- copulaMarginsQuantile(u, margins)
eta <- joint[, 1:2, drop = FALSE]
weight <- joint[, 3]
model_definition <- MODELS$iv1
psi <- sweep(exp(eta), 2L, model_definition$true, "*")
nt <- length(model_definition$times)
dataset <- data.frame(id = rep(seq_len(N), each = nt),
                      dose = model_definition$dose,
                      time = rep(model_definition$times, N))
prediction <- model_definition$f(
  psi, dataset$id, cbind(dataset$dose, dataset$time))
dataset$y <- pmax(prediction * (1 + .08 * rnorm(nrow(dataset))), 1e-6)

conditioning <- matrix(weight, ncol = 1L, dimnames = list(NULL, "WT"))
conditioning[seq(4L, N, by = 4L), 1L] <- NA_real_
copulaSet(vine, margins = margins,
  conditioning = list(values = conditioning, variableName = "WT"),
  warmStartOnActivate = FALSE, familySet = NULL,
  refitEvery = 2L, jointMaxit = 8L, jointFinalMaxit = 15L,
  guard = FALSE)

control <- list(seed = 82602, save = FALSE, save.graphs = FALSE,
  print = FALSE, displayProgress = FALSE,
  nbiter.saemix = c(20, 10), nbiter.mcmc = c(1, 1, 1, 0),
  warnings = FALSE, ll.is = FALSE, fim = FALSE, map = FALSE)
fit <- saemix(saemixModelFor(model_definition),
             saemixDataFor(dataset), control)
state <- copulaGet(fit)
stopifnot(inherits(fit, "SaemixObject"),
          inherits(state, "saemixCopulaSnapshot"),
          identical(state$lastJoint$backend,
                    "gaussian-copula-frem-mvnormal-em"),
          all(is.finite(fit@results@fixed.effects)),
          all(is.finite(state$sd)),
          anyNA(state$conditioning))

## The affine/conjugate likelihood shortcut must fail closed for this nonlinear
## exponential-parameter PK model with proportional residual error.
fit@options$nmc.is <- 200L
likelihoodFit <- suppressWarnings(llisCopula.saemix(fit, defensive = 0,
  batch = 50L, seed = 82603L))
likelihoodDiagnostic <- attr(likelihoodFit, "saemix.copula.likelihood")
stopifnot(is.finite(likelihoodFit@results@ll.is),
  !isTRUE(likelihoodDiagnostic$exact))

## The affine/constant-error submodel must take the conjugate observed-
## likelihood path and match an independently assembled Gaussian integral.
fixedLin <- c(2.1, -.16)
psiLin <- sweep(eta, 2L, fixedLin, "+")
datasetLin <- dataset
datasetLin$y <- psiLin[datasetLin$id, 1L] +
  psiLin[datasetLin$id, 2L] * datasetLin$time +
  .12 * rnorm(nrow(datasetLin))
linear <- function(psi, id, xidep)
  psi[id, 1L] + psi[id, 2L] * xidep[, 2L]
linearModel <- saemixModel(model = linear, modeltype = "structural",
  psi0 = matrix(fixedLin, nrow = 1L,
    dimnames = list(NULL, c("BASE", "SLOPE"))),
  transform.par = c(0, 0), covariance.model = matrix(1, 2L, 2L),
  omega.init = diag(c(.3^2, .15^2)), error.model = "constant",
  error.init = c(.12, 0), verbose = FALSE)
linearData <- saemixData(name.data = datasetLin, header = TRUE,
  name.group = "id", name.predictors = c("dose", "time"),
  name.response = "y", verbose = FALSE)
copulaSet(vine, margins = margins,
  conditioning = list(values = conditioning, variableName = "WT"),
  warmStartOnActivate = FALSE, familySet = NULL, refitEvery = 2L,
  jointMaxit = 8L, jointFinalMaxit = 15L, guard = FALSE)
linearFit <- saemix(linearModel, linearData, control)
linearFit <- llisCopula.saemix(linearFit, defensive = 0, seed = 82604L)
linearDiagnostic <- attr(linearFit, "saemix.copula.likelihood")
stopifnot(isTRUE(linearDiagnostic$exact),
  identical(linearDiagnostic$method, "exact-linear-gaussian"),
  linearDiagnostic$se_loglik_total == 0)

linearState <- copulaGet(linearFit)
conditionalLin <- copulaGaussianFremConditional(conditioning,
  linearState$vine, linearState$margins, 2L)
etaScale <- diag(vapply(linearState$margins[1:2], function(m)
  m$scale(m$parameters), numeric(1L)))
predictorLin <- linearFit@results@mean.phi[, 1:2, drop = FALSE]
responseReference <- numeric(N)
for (i in seq_len(N)) {
  priorMean <- predictorLin[i, ] +
    as.numeric(etaScale %*% conditionalLin$mean[i, ])
  priorCov <- etaScale %*% conditionalLin$covariance[[i]] %*% etaScale
  rows <- datasetLin$id == i
  design <- cbind(1, datasetLin$time[rows])
  responseCov <- design %*% priorCov %*% t(design) +
    diag(linearFit@results@respar[1L]^2, sum(rows))
  responseReference[i] <- copulaGaussianLogDensity(matrix(
    datasetLin$y[rows] - as.numeric(design %*% priorMean), nrow = 1L),
    responseCov)
}
jointReference <- responseReference +
  copulaGaussianFremConditioningLogDensity(conditioning,
    linearState$vine, linearState$margins, 2L)
stopifnot(max(abs(linearDiagnostic$per_subject_loglik - jointReference)) < 1e-10,
  abs(linearFit@results@ll.is - sum(jointReference)) < 1e-10)

cat("Gaussian-copula FREM missing-covariate end-to-end check passed\n")
