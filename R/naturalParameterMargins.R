## Natural-scale individual-parameter margins. The distribution is independent
## of the saemix working-coordinate transform; mapping and Jacobians are applied
## separately by copulaNaturalToWorking()/copulaWorkingToNatural().

copulaNaturalMargin <- function(name, support, anchor, parameters, lower, upper,
    log_density, cdf, quantile, random = NULL,
    free = rep(TRUE, length(parameters)), roles = NULL, metadata = list()) {
  support <- match.arg(support, c("real", "positive", "unit"))
  anchor <- match.arg(anchor, c("mean", "median", "geometric-mean", "custom"))
  parameters <- as.numeric(parameters) |>
    stats::setNames(names(parameters) %||% paste0("par", seq_along(parameters)))
  if (is.null(roles)) roles <- rep("shape", length(parameters))
  margin <- structure(list(name = as.character(name), support = support,
    anchor = anchor, type = "continuous", parameters = parameters,
    lower = as.numeric(lower), upper = as.numeric(upper),
    free = as.logical(free), roles = as.character(roles),
    log_density = log_density, cdf = cdf, quantile = quantile,
    random = random, metadata = metadata),
    class = "saemix_natural_parameter_margin")
  copulaNaturalMarginValidate(margin)
  margin
}

copulaNaturalMarginValidate <- function(margin, probe = TRUE) {
  if (!inherits(margin, "saemix_natural_parameter_margin"))
    stop("margin must inherit from saemix_natural_parameter_margin")
  p <- margin$parameters
  if (length(p) != length(margin$lower) || length(p) != length(margin$upper) ||
      length(p) != length(margin$free) || length(p) != length(margin$roles) ||
      any(!is.finite(p)) || any(p < margin$lower | p > margin$upper) ||
      !all(vapply(margin[c("log_density", "cdf", "quantile")],
        is.function, logical(1)))) stop("invalid natural parameter margin")
  if (isTRUE(probe)) {
    typical <- rep(switch(margin$support, real = .4, positive = 2, unit = .4), 3L)
    u <- c(.2, .5, .8)
    q <- margin$quantile(u, typical, p)
    if (any(!is.finite(q)) || any(!is.finite(margin$log_density(q, typical, p))) ||
        max(abs(margin$cdf(q, typical, p) - u)) > 1e-7)
      stop("natural parameter margin failed density/CDF/quantile probe")
  }
  invisible(TRUE)
}

copulaNaturalMarginWithParameters <- function(margin, parameters) {
  parameters <- stats::setNames(as.numeric(parameters), names(margin$parameters))
  if (length(parameters) != length(margin$parameters) ||
      any(!is.finite(parameters)) || any(parameters < margin$lower) ||
      any(parameters > margin$upper)) stop("natural margin parameters violate bounds")
  margin$parameters <- parameters
  copulaNaturalMarginValidate(margin, probe = FALSE)
  margin
}

copulaNaturalMarginNormal <- function(sd = 1) {
  copulaNaturalMargin("normal", "real", "mean", c(sd = sd), 1e-8, 1e4,
    roles = "scale",
    log_density = function(x, typical, par)
      stats::dnorm(x, typical, par["sd"], log = TRUE),
    cdf = function(x, typical, par) stats::pnorm(x, typical, par["sd"]),
    quantile = function(u, typical, par) stats::qnorm(u, typical, par["sd"]),
    random = function(n, typical, par) stats::rnorm(n, typical, par["sd"]))
}

copulaNaturalMarginStudent <- function(sd = 1, df = 6) {
  rawScale <- function(par) par["sd"] * sqrt((par["df"] - 2) / par["df"])
  copulaNaturalMargin("student", "real", "mean", c(sd = sd, df = df),
    c(1e-8, 2 + 1e-6), c(1e4, 200), roles = c("scale", "shape"),
    log_density = function(x, typical, par) {
      s <- rawScale(par)
      stats::dt((x - typical) / s, par["df"], log = TRUE) - log(s)
    },
    cdf = function(x, typical, par)
      stats::pt((x - typical) / rawScale(par), par["df"]),
    quantile = function(u, typical, par)
      typical + rawScale(par) * stats::qt(u, par["df"]),
    random = function(n, typical, par)
      typical + rawScale(par) * stats::rt(n, par["df"]))
}

copulaNaturalMarginLaplace <- function(sd = 1) {
  b <- function(par) par["sd"] / sqrt(2)
  copulaNaturalMargin("laplace", "real", "mean", c(sd = sd), 1e-8, 1e4,
    roles = "scale",
    log_density = function(x, typical, par)
      -log(2 * b(par)) - abs(x - typical) / b(par),
    cdf = function(x, typical, par) {
      z <- (x - typical) / b(par)
      ifelse(z < 0, .5 * exp(z), 1 - .5 * exp(-z))
    },
    quantile = function(u, typical, par)
      typical + ifelse(u < .5, b(par) * log(2 * u),
        -b(par) * log(2 * (1 - u))),
    random = function(n, typical, par) {
      u <- stats::runif(n)
      typical + ifelse(u < .5, b(par) * log(2 * u),
        -b(par) * log(2 * (1 - u)))
    })
}

copulaNaturalMarginLognormal <- function(sdlog = .3,
    anchor = c("median", "mean", "geometric-mean")) {
  anchor <- match.arg(anchor)
  meanlog <- function(typical, par) log(typical) -
    if (anchor == "mean") .5 * par["sdlog"]^2 else 0
  copulaNaturalMargin("lognormal", "positive", anchor, c(sdlog = sdlog),
    1e-8, 5, roles = "scale",
    log_density = function(x, typical, par)
      stats::dlnorm(x, meanlog(typical, par), par["sdlog"], log = TRUE),
    cdf = function(x, typical, par)
      stats::plnorm(x, meanlog(typical, par), par["sdlog"]),
    quantile = function(u, typical, par)
      stats::qlnorm(u, meanlog(typical, par), par["sdlog"]),
    random = function(n, typical, par)
      stats::rlnorm(n, meanlog(typical, par), par["sdlog"]))
}

copulaNaturalMarginGamma <- function(shape = 4,
    anchor = c("geometric-mean", "median", "mean")) {
  anchor <- match.arg(anchor)
  scaleValue <- function(typical, par) switch(anchor,
    mean = typical / par["shape"],
    median = typical / stats::qgamma(.5, par["shape"], scale = 1),
    `geometric-mean` = typical * exp(-digamma(par["shape"])))
  copulaNaturalMargin("gamma", "positive", anchor, c(shape = shape),
    .05, 1e4, roles = "shape",
    log_density = function(x, typical, par)
      stats::dgamma(x, par["shape"], scale = scaleValue(typical, par),
        log = TRUE),
    cdf = function(x, typical, par)
      stats::pgamma(x, par["shape"], scale = scaleValue(typical, par)),
    quantile = function(u, typical, par)
      stats::qgamma(u, par["shape"], scale = scaleValue(typical, par)),
    random = function(n, typical, par)
      stats::rgamma(n, par["shape"], scale = scaleValue(typical, par)))
}

copulaNaturalMarginWeibull <- function(shape = 3,
    anchor = c("geometric-mean", "median", "mean")) {
  anchor <- match.arg(anchor); eulerGamma <- -digamma(1)
  scaleValue <- function(typical, par) switch(anchor,
    mean = typical / gamma(1 + 1 / par["shape"]),
    median = typical / log(2)^(1 / par["shape"]),
    `geometric-mean` = typical * exp(eulerGamma / par["shape"]))
  copulaNaturalMargin("weibull", "positive", anchor, c(shape = shape),
    .05, 1e3, roles = "shape",
    log_density = function(x, typical, par)
      stats::dweibull(x, par["shape"], scaleValue(typical, par), log = TRUE),
    cdf = function(x, typical, par)
      stats::pweibull(x, par["shape"], scaleValue(typical, par)),
    quantile = function(u, typical, par)
      stats::qweibull(u, par["shape"], scaleValue(typical, par)),
    random = function(n, typical, par)
      stats::rweibull(n, par["shape"], scaleValue(typical, par)))
}

copulaWorkingToNatural <- function(phi, transform) {
  transform <- rep_len(as.integer(transform), ncol(as.matrix(phi)))
  transphi(as.matrix(phi), transform)
}

copulaNaturalToWorking <- function(psi, transform) {
  transform <- rep_len(as.integer(transform), ncol(as.matrix(psi)))
  transpsi(as.matrix(psi), transform)
}

copulaWorkingLogJacobian <- function(phi, transform) {
  phi <- as.matrix(phi); transform <- rep_len(as.integer(transform), ncol(phi))
  answer <- matrix(0, nrow(phi), ncol(phi))
  for (j in seq_len(ncol(phi))) answer[, j] <- switch(
    as.character(transform[j]),
    `0` = 0,
    `1` = phi[, j],
    `2` = stats::dnorm(phi[, j], log = TRUE),
    `3` = {
      p <- stats::plogis(phi[, j]); log(p) + log1p(-p)
    }, stop("unsupported saemix parameter transform"))
  answer
}

copulaNaturalMarginsEvaluate <- function(psi, typical, margins) {
  psi <- as.matrix(psi); typical <- as.matrix(typical)
  if (any(dim(psi) != dim(typical)) || ncol(psi) != length(margins))
    stop("natural values, typical values, and margins do not align")
  z <- logMargin <- matrix(NA_real_, nrow(psi), ncol(psi))
  valid <- rep(TRUE, nrow(psi))
  for (j in seq_along(margins)) {
    margin <- margins[[j]]
    density <- margin$log_density(psi[, j], typical[, j], margin$parameters)
    probability <- margin$cdf(psi[, j], typical[, j], margin$parameters)
    ok <- is.finite(density) & is.finite(probability) &
      probability > 0 & probability < 1
    valid <- valid & ok
    z[ok, j] <- stats::qnorm(probability[ok]); logMargin[ok, j] <- density[ok]
  }
  list(z = z, logMargin = logMargin, valid = valid)
}

copulaNaturalMarginsQuantile <- function(u, typical, margins) {
  u <- as.matrix(u); typical <- as.matrix(typical)
  if (any(dim(u) != dim(typical)) || ncol(u) != length(margins))
    stop("uniforms, typical values, and margins do not align")
  answer <- u
  for (j in seq_along(margins)) answer[, j] <- margins[[j]]$quantile(
    u[, j], typical[, j], margins[[j]]$parameters)
  answer
}

copulaNaturalWorkingPriorKernel <- function(vine, margins, predictor,
                                            transform) {
  predictor <- as.matrix(predictor); d <- ncol(predictor)
  transform <- rep_len(as.integer(transform), d)
  if (length(margins) != d || !copulaIsFullGaussianVine(vine, d))
    stop("invalid natural-parameter prior kernel")
  typical <- copulaWorkingToNatural(predictor, transform)
  correlation <- copulaGaussianRvineCor(vine, d)
  U <- chol(correlation); logDet <- 2 * sum(log(diag(U)))
  negative <- function(eta) {
    eta <- as.matrix(eta)
    if (any(dim(eta) != dim(predictor)))
      stop("natural-parameter eta batch does not align with its predictor")
    phi <- predictor + eta
    psi <- copulaWorkingToNatural(phi, transform)
    evaluated <- copulaNaturalMarginsEvaluate(psi, typical, margins)
    answer <- rep(Inf, nrow(eta)); rows <- which(evaluated$valid)
    if (length(rows)) {
      z <- evaluated$z[rows, , drop = FALSE]
      standardized <- forwardsolve(t(U), t(z))
      logGaussian <- -.5 * (d * log(2 * pi) + logDet +
        colSums(standardized^2))
      logDensity <- logGaussian - rowSums(stats::dnorm(z, log = TRUE)) +
        rowSums(evaluated$logMargin[rows, , drop = FALSE]) +
        rowSums(copulaWorkingLogJacobian(phi[rows, , drop = FALSE], transform))
      answer[rows] <- -logDensity
    }
    answer
  }
  random <- function() {
    u <- rvinecopulib::rvinecop(nrow(predictor), vine)
    psi <- copulaNaturalMarginsQuantile(u, typical, margins)
    copulaNaturalToWorking(psi, transform) - predictor
  }
  list(negative = negative, random = random, predictor = predictor,
    typical = typical, transform = transform,
    method = "natural-parameter-working-coordinate-prior")
}
