## Public constructor for the fixed Gaussian-copula score model.

gaussianCopulaFrem <- function(etaSd = NULL, etaMargins = NULL,
                               covariates = NULL,
                               covariateMargins = "auto",
                               correlation = NULL, structure = NULL,
                               covariateNames = NULL,
                               scoreScale = 0.05, scoreBurn = 50L,
                               gainScale = 0.2, gainPower = 0.8,
                               finiteDifference = 1e-4,
                               projection = 24) {
  if (is.null(etaMargins)) {
    if (is.null(etaSd) || !length(etaSd) ||
        any(!is.finite(etaSd)) || any(etaSd <= 0))
      stop("supply positive etaSd or one etaMargins object per random effect")
    etaMargins <- lapply(as.numeric(etaSd), copulaMarginNormal)
  } else {
    if (!is.null(etaSd))
      stop("supply etaSd or etaMargins, not both")
    if (!is.list(etaMargins) || !length(etaMargins))
      stop("etaMargins must be a non-empty list")
    lapply(etaMargins, copulaMarginValidate)
  }

  conditioning <- NULL
  covMargins <- list()
  if (!is.null(covariates)) {
    covariates <- as.matrix(covariates)
    storage.mode(covariates) <- "double"
    if (!nrow(covariates) || !ncol(covariates) || any(is.infinite(covariates)))
      stop("covariates must be a non-empty numeric matrix; NA is allowed")
    if (is.character(covariateMargins)) {
      if (length(covariateMargins) != 1L ||
          !identical(tolower(covariateMargins), "auto"))
        stop("character covariateMargins must be 'auto'")
      covMargins <- copulaFitCovariateMargins(covariates)
    } else {
      if (!is.list(covariateMargins) ||
          length(covariateMargins) != ncol(covariates))
        stop("covariateMargins must contain one margin per covariate")
      lapply(covariateMargins, copulaMarginValidate)
      covMargins <- covariateMargins
    }
    if (is.null(covariateNames)) covariateNames <- colnames(covariates)
    if (is.null(covariateNames) || length(covariateNames) != ncol(covariates))
      covariateNames <- paste0("COV", seq_len(ncol(covariates)))
    conditioning <- list(values = covariates,
      variableName = as.character(covariateNames))
  } else if (!identical(covariateMargins, "auto")) {
    stop("covariateMargins requires covariates")
  }

  margins <- c(etaMargins, covMargins)
  d <- length(margins)
  if (d < 2L)
    stop("a copula population requires at least two joint coordinates")
  if (is.null(correlation)) correlation <- diag(d)
  correlation <- copulaValidateCorrelation(correlation)
  if (nrow(correlation) != d)
    stop("correlation dimension must equal eta plus covariate coordinates")
  if (is.null(structure))
    structure <- rvinecopulib::cvine_structure(seq_len(d))
  vine <- copulaGaussianRvineFromCor(correlation, structure)

  copulaPopulation(vine, margins = margins,
    scale = "transformed-additive", conditioning = conditioning,
    populationAlgorithm = "score-sa",
    scoreScale = scoreScale, scoreBurn = scoreBurn,
    scoreGainScale = gainScale, scoreGainPower = gainPower,
    scoreFiniteDifference = finiteDifference,
    scoreProjection = projection)
}
