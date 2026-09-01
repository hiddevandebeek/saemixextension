## Support-based registry for automatic natural-parameter marginal screening.

.naturalMarginRegistry <- new.env(parent = emptyenv())

copulaRegisterNaturalMargin <- function(name, support, start,
                                        replace = FALSE) {
  name <- tolower(as.character(name)); support <- match.arg(support,
    c("real", "positive", "unit"))
  if (length(name) != 1L || !nzchar(name) || !is.function(start))
    stop("invalid natural-margin registry entry")
  if (exists(name, .naturalMarginRegistry, inherits = FALSE) && !replace)
    stop("natural-margin family is already registered: ", name)
  assign(name, list(name = name, support = support, start = start),
    .naturalMarginRegistry)
  invisible(name)
}

copulaNaturalMarginRegistry <- function() {
  if (!exists(".builtins", .naturalMarginRegistry, inherits = FALSE)) {
    assign(".builtins", TRUE, .naturalMarginRegistry)
    logDeviationSd <- function(values, typical)
      stats::sd(log(values / typical))
    copulaRegisterNaturalMargin("normal", "real", function(values, typical)
      copulaNaturalMarginNormal(max(stats::sd(values - typical), 1e-3)), TRUE)
    copulaRegisterNaturalMargin("student", "real", function(values, typical)
      copulaNaturalMarginStudent(max(stats::sd(values - typical), 1e-3), 6), TRUE)
    copulaRegisterNaturalMargin("laplace", "real", function(values, typical)
      copulaNaturalMarginLaplace(max(stats::sd(values - typical), 1e-3)), TRUE)
    copulaRegisterNaturalMargin("lognormal", "positive", function(values, typical)
      copulaNaturalMarginLognormal(max(logDeviationSd(values, typical), 1e-3)), TRUE)
    copulaRegisterNaturalMargin("gamma", "positive", function(values, typical) {
      target <- max(logDeviationSd(values, typical)^2, 1e-6)
      shape <- stats::uniroot(function(x) trigamma(x) - target,
        c(.05, 1e4))$root
      copulaNaturalMarginGamma(shape)
    }, TRUE)
    copulaRegisterNaturalMargin("weibull", "positive", function(values, typical) {
      shape <- pi / (sqrt(6) * max(logDeviationSd(values, typical), 1e-3))
      copulaNaturalMarginWeibull(shape)
    }, TRUE)
  }
  entries <- as.list(.naturalMarginRegistry)
  entries[setdiff(names(entries), ".builtins")]
}

copulaNaturalMarginFamilies <- function(support = c("real", "positive", "unit")) {
  support <- match.arg(support)
  registry <- copulaNaturalMarginRegistry()
  names(Filter(function(entry) identical(entry$support, support), registry))
}

copulaNaturalMarginStart <- function(family, values, typical, support) {
  registry <- copulaNaturalMarginRegistry(); family <- tolower(family)
  entry <- registry[[family]]
  if (is.null(entry)) stop("unsupported natural parameter family: ", family)
  if (!identical(entry$support, support))
    stop(family, " is incompatible with declared support ", support)
  margin <- entry$start(as.numeric(values), as.numeric(typical))
  copulaNaturalMarginValidate(margin)
  margin
}
