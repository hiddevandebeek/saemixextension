## Marginal-distribution interface for copula random effects.
##
## A margin is deliberately independent of the vine.  It owns a proper
## univariate law and a finite native parameter vector.  The population
## location remains X beta, as in ordinary saemix, so built-in eta margins are
## centred at zero to avoid duplicating the intercept.

copulaMargin <- function(name, parameters, lower, upper,
                         free = rep(TRUE, length(parameters)),
                         log_density, cdf, quantile,
                         random = NULL, type = c("continuous", "discrete"),
                         cdf_left = NULL, scale = NULL,
                         set_scale = NULL,
                         centered = TRUE, roles = NULL, scale_is_sd = FALSE,
                         metadata = list(), support_fixed = NULL) {
  type <- match.arg(type)
  parameter_names <- names(parameters)
  parameters <- as.numeric(parameters)
  names(parameters) <- parameter_names %||% paste0("par", seq_along(parameters))
  lower <- as.numeric(lower); upper <- as.numeric(upper)
  free <- as.logical(free)
  if (is.null(roles)) {
    roles <- rep("shape", length(parameters))
    roles[tolower(names(parameters)) %in% c("mean", "mu", "location", "loc")] <-
      "location"
  }
  roles <- as.character(roles)
  if (length(parameters) != length(lower) || length(lower) != length(upper))
    stop("copulaMargin: parameters, lower, and upper must have equal length")
  if (length(free) != length(parameters) || anyNA(free))
    stop("copulaMargin: free must be one logical value per parameter")
  if (length(roles) != length(parameters) ||
      any(!roles %in% c("location", "scale", "shape")))
    stop("copulaMargin: roles must be location, scale, or shape")
  if (any(!is.finite(parameters)) || any(parameters < lower) || any(parameters > upper))
    stop("copulaMargin: initial parameters must be finite and within bounds")
  if (!is.function(log_density) || !is.function(cdf) || !is.function(quantile))
    stop("copulaMargin: log_density, cdf, and quantile must be functions")
  if (type == "discrete" && !is.function(cdf_left))
    stop("copulaMargin: a discrete margin requires cdf_left")
  if (!is.null(random) && !is.function(random))
    stop("copulaMargin: random must be NULL or a function")
  if (is.null(scale)) scale <- function(par) NA_real_
  if (!is.function(scale)) stop("copulaMargin: scale must be a function")
  if (!is.null(set_scale) && !is.function(set_scale))
    stop("copulaMargin: set_scale must be NULL or a function")
  if (!is.null(support_fixed)) {
    if (length(support_fixed) != 1L || is.na(support_fixed))
      stop("copulaMargin: support_fixed must be TRUE, FALSE, or NULL")
    metadata$parameter_independent_support <- isTRUE(support_fixed)
    if (!isTRUE(support_fixed) && identical(type, "continuous"))
      metadata$fixed_reference_quantile <- TRUE
  }
  structure(list(
    name = as.character(name), type = type, parameters = parameters, free = free,
    lower = lower, upper = upper, log_density = log_density, cdf = cdf,
    cdf_left = cdf_left, quantile = quantile, random = random, scale = scale,
    set_scale = set_scale,
    centered = isTRUE(centered), roles = roles,
    scale_is_sd = isTRUE(scale_is_sd), metadata = metadata),
    class = "saemix_copula_margin")
}

## Public adapter for a user-defined continuous eta distribution. The supplied
## callbacks describe a raw distribution X; the returned margin describes
## eta = X - E(X), keeping the saemix population location in C %*% mu. `scale`
## must return the actual marginal standard deviation for every parameter
## value. The caller must explicitly declare whether the centered native
## support is parameter independent.
copulaMarginCenteredCustom <- function(
    name, parameters, lower, upper, log_density, cdf, quantile,
    center, scale, support_fixed,
    random = NULL, set_scale = NULL,
    free = rep(TRUE, length(parameters)), roles = NULL,
    metadata = list()) {
  callbacks <- list(log_density = log_density, cdf = cdf,
    quantile = quantile, center = center, scale = scale)
  if (any(!vapply(callbacks, is.function, logical(1))))
    stop("copulaMarginCenteredCustom: density, CDF, quantile, center, and scale callbacks are required")
  if (length(support_fixed) != 1L || is.na(support_fixed))
    stop("copulaMarginCenteredCustom: support_fixed must be TRUE or FALSE")
  if (!is.null(random) && !is.function(random))
    stop("copulaMarginCenteredCustom: random must be NULL or a function")
  if (is.null(roles)) {
    roles <- rep("shape", length(parameters))
    roles[tolower(names(parameters)) %in% c("sd", "sigma", "scale")] <- "scale"
  }
  centered_log_density <- function(eta, par)
    log_density(eta + center(par), par)
  centered_cdf <- function(eta, par)
    cdf(eta + center(par), par)
  centered_quantile <- function(u, par)
    quantile(u, par) - center(par)
  centered_random <- if (is.null(random)) NULL else function(n, par)
    random(n, par) - center(par)
  initialCenter <- center(parameters)
  initialScale <- scale(parameters)
  if (length(initialCenter) != 1L || !is.finite(initialCenter) ||
      length(initialScale) != 1L || !is.finite(initialScale) || initialScale <= 0)
    stop("copulaMarginCenteredCustom: center and scale must return one finite value; scale must be positive")
  metadata$adapter <- "centered-custom-continuous"
  metadata$centering_definition <- "eta = raw variate - center(parameters)"
  metadata$scale_definition <- "marginal standard deviation"
  margin <- copulaMargin(name, parameters, lower, upper, free = free,
    log_density = centered_log_density, cdf = centered_cdf,
    quantile = centered_quantile, random = centered_random,
    type = "continuous", scale = scale, set_scale = set_scale,
    centered = TRUE, roles = roles, scale_is_sd = TRUE,
    metadata = metadata, support_fixed = support_fixed)
  copulaMarginValidate(margin)
  margin
}

`%||%` <- function(x, y) if (is.null(x)) y else x

copulaMarginValidate <- function(margin, probe = TRUE) {
  if (!inherits(margin, "saemix_copula_margin"))
    stop("margin must inherit from saemix_copula_margin")
  p <- margin$parameters
  if (length(p) != length(margin$lower) || length(p) != length(margin$upper) ||
      length(p) != length(margin$free) ||
      length(p) != length(margin$roles) ||
      any(!is.finite(p)) || any(p < margin$lower) || any(p > margin$upper))
    stop("invalid marginal parameter vector or bounds")
  if (isTRUE(probe)) {
    q <- margin$quantile(c(.2, .5, .8), p)
    if (length(q) != 3L || any(!is.finite(q)))
      stop("margin quantile probe returned invalid values")
    Fq <- margin$cdf(q, p)
    ld <- margin$log_density(q, p)
    if (length(Fq) != 3L || any(!is.finite(Fq)) || any(Fq < 0 | Fq > 1) ||
        length(ld) != 3L || any(!is.finite(ld)))
      stop("margin density/CDF probe returned invalid values")
  }
  invisible(TRUE)
}

copulaMarginCenteredGamma <- function(shape = 4, sd = 1) {
  scaleFrom <- function(par) par["sd"] / sqrt(par["shape"])
  shift <- function(par) par["shape"] * scaleFrom(par)
  copulaMargin("gamma-centered", c(shape = unname(shape), sd = unname(sd)),
    c(1e-4, 1e-6), c(1e6, 1e6),
    log_density = function(x, par)
      stats::dgamma(x + shift(par), shape = par["shape"],
        scale = scaleFrom(par), log = TRUE),
    cdf = function(x, par)
      stats::pgamma(x + shift(par), shape = par["shape"],
        scale = scaleFrom(par)),
    quantile = function(u, par)
      stats::qgamma(u, shape = par["shape"], scale = scaleFrom(par)) -
        shift(par),
    random = function(n, par)
      stats::rgamma(n, shape = par["shape"], scale = scaleFrom(par)) -
        shift(par),
    scale = function(par) par["sd"], centered = TRUE,
    set_scale = function(par, value) { par["sd"] <- value; par },
    roles = c("shape", "scale"), scale_is_sd = TRUE,
    metadata = list(scale_definition = "standard deviation"),
    support_fixed = FALSE)
}

copulaMarginWithParameters <- function(margin, parameters) {
  parameters <- stats::setNames(as.numeric(parameters), names(margin$parameters))
  if (length(parameters) != length(margin$parameters) ||
      any(!is.finite(parameters)) || any(parameters < margin$lower) ||
      any(parameters > margin$upper))
    stop("marginal parameters are outside their native bounds")
  margin$parameters <- parameters
  copulaMarginValidate(margin, probe = FALSE)
  margin
}

copulaMarginNormal <- function(sd = 1) {
  copulaMargin("normal", c(sd = unname(sd)), 1e-6, 1e3,
    log_density = function(x, par) stats::dnorm(x, sd = par[1], log = TRUE),
    cdf = function(x, par) stats::pnorm(x, sd = par[1]),
    quantile = function(u, par) stats::qnorm(u, sd = par[1]),
    random = function(n, par) stats::rnorm(n, sd = par[1]),
    scale = function(par) par[1], centered = TRUE,
    set_scale = function(par, value) { par[1] <- value; par },
    roles = "scale", scale_is_sd = TRUE,
    metadata = list(scale_definition = "standard deviation"),
    support_fixed = TRUE)
}

## Centred Student margin parameterised by its finite standard deviation.
copulaMarginStudent <- function(sd = 1, df = 6) {
  raw_scale <- function(par) par[1] * sqrt((par[2] - 2) / par[2])
  copulaMargin("student", c(sd = unname(sd), df = unname(df)),
    c(1e-6, 2 + 1e-6), c(1e3, 200),
    log_density = function(x, par) {
      s <- raw_scale(par); stats::dt(x / s, df = par[2], log = TRUE) - log(s)
    },
    cdf = function(x, par) stats::pt(x / raw_scale(par), df = par[2]),
    quantile = function(u, par) raw_scale(par) * stats::qt(u, df = par[2]),
    random = function(n, par) raw_scale(par) * stats::rt(n, df = par[2]),
    scale = function(par) par[1], centered = TRUE,
    set_scale = function(par, value) { par[1] <- value; par },
    roles = c("scale", "shape"), scale_is_sd = TRUE,
    metadata = list(scale_definition = "standard deviation"),
    support_fixed = TRUE)
}

copulaMarginLaplace <- function(sd = 1) {
  b <- function(par) par[1] / sqrt(2)
  copulaMargin("laplace", c(sd = unname(sd)), 1e-6, 1e3,
    log_density = function(x, par) -log(2 * b(par)) - abs(x) / b(par),
    cdf = function(x, par) ifelse(x < 0, .5 * exp(x / b(par)),
                                  1 - .5 * exp(-x / b(par))),
    quantile = function(u, par) ifelse(u < .5, b(par) * log(2 * u),
                                       -b(par) * log(2 * (1 - u))),
    random = function(n, par) {
      u <- stats::runif(n); ifelse(u < .5, b(par) * log(2 * u),
                                   -b(par) * log(2 * (1 - u)))
    },
    scale = function(par) par[1], centered = TRUE,
    set_scale = function(par, value) { par[1] <- value; par },
    roles = "scale", scale_is_sd = TRUE,
    metadata = list(scale_definition = "standard deviation"),
    support_fixed = TRUE)
}

copulaCategoricalProbabilities <- function(logits) {
  logits <- as.numeric(logits)
  if (!length(logits)) return(1)
  stick <- stats::plogis(logits); remaining <- 1
  probability <- numeric(length(logits) + 1L)
  for (j in seq_along(logits)) {
    probability[j] <- remaining * stick[j]
    remaining <- remaining * (1 - stick[j])
  }
  probability[length(probability)] <- remaining
  probability / sum(probability)
}

## Proper ordered categorical margin. Native optimizer coordinates are
## stick-breaking logits, which enforce positive probabilities summing to one
## without constrained simplex optimization. The category labels may be
## numeric or character; their supplied order defines the ordinal order.
copulaMarginOrdinal <- function(probabilities, labels = seq_along(probabilities),
                                free = TRUE, name = "ordinal") {
  probabilities <- as.numeric(probabilities); labels <- as.character(labels)
  if (length(probabilities) < 2L || length(labels) != length(probabilities) ||
      any(!is.finite(probabilities)) || any(probabilities <= 0) ||
      abs(sum(probabilities) - 1) > 1e-8 || anyDuplicated(labels))
    stop("ordinal margin requires distinct labels and positive probabilities summing to one")
  remaining <- 1; logits <- numeric(length(probabilities) - 1L)
  for (j in seq_along(logits)) {
    fraction <- probabilities[j] / remaining
    logits[j] <- stats::qlogis(fraction)
    remaining <- remaining - probabilities[j]
  }
  names(logits) <- paste0("logit", seq_along(logits))
  categoryIndex <- function(x) match(as.character(x), labels)
  probability <- function(par) copulaCategoricalProbabilities(par)
  cdfValue <- function(x, par, left = FALSE) {
    index <- categoryIndex(x); p <- probability(par); cumulative <- cumsum(p)
    out <- rep(NA_real_, length(x)); valid <- !is.na(index)
    out[valid] <- if (left) c(0, cumulative[-length(cumulative)])[index[valid]] else
      cumulative[index[valid]]
    out
  }
  numericLabels <- suppressWarnings(as.numeric(labels))
  scaleFun <- function(par) {
    score <- if (anyNA(numericLabels)) seq_along(labels) else numericLabels
    p <- probability(par); mu <- sum(p * score)
    sqrt(sum(p * (score - mu)^2))
  }
  margin <- copulaMargin(name, logits, rep(-25, length(logits)),
    rep(25, length(logits)), free = rep(isTRUE(free), length(logits)),
    log_density = function(x, par) {
      index <- categoryIndex(x); out <- rep(-Inf, length(x)); valid <- !is.na(index)
      out[valid] <- log(probability(par)[index[valid]]); out
    },
    cdf = function(x, par) cdfValue(x, par, FALSE),
    cdf_left = function(x, par) cdfValue(x, par, TRUE),
    quantile = function(u, par) {
      p <- cumsum(probability(par)); index <- findInterval(u, c(0, p),
        left.open = TRUE, all.inside = TRUE)
      value <- labels[pmin(index, length(labels))]
      if (all(!is.na(numericLabels))) as.numeric(value) else value
    },
    random = function(n, par) {
      value <- sample(labels, n, replace = TRUE, prob = probability(par))
      if (all(!is.na(numericLabels))) as.numeric(value) else value
    },
    type = "discrete", scale = scaleFun, centered = FALSE,
    roles = rep("shape", length(logits)), scale_is_sd = TRUE,
    metadata = list(variable_role = "conditioning", categorical = TRUE,
      categorical_kind = "ordinal", labels = labels,
      parameterization = "stick-breaking logits"), support_fixed = TRUE)
  margin
}

copulaMarginBernoulli <- function(prob = .5, labels = c(0, 1), free = TRUE) {
  if (length(prob) != 1L || !is.finite(prob) || prob <= 0 || prob >= 1)
    stop("Bernoulli probability must be in (0,1)")
  margin <- copulaMarginOrdinal(c(1 - prob, prob), labels, free,
    name = "bernoulli")
  margin$metadata$categorical_kind <- "binary"
  margin
}

## Natural-scale continuous covariate margins.  These are deliberately not
## centred: unlike additive eta coordinates, observed conditioning variables
## own their locations.  `scale` reports the actual marginal standard
## deviation so Gaussian-copula initialization and diagnostics have a finite
## surrogate scale without imposing a Normal marginal.
copulaMarginCovariateNormal <- function(mean, sd) {
  copulaMargin("normal", c(mean = unname(mean), sd = unname(sd)),
    c(-1e8, 1e-8), c(1e8, 1e8),
    log_density = function(x, par)
      stats::dnorm(x, mean = par["mean"], sd = par["sd"], log = TRUE),
    cdf = function(x, par)
      stats::pnorm(x, mean = par["mean"], sd = par["sd"]),
    quantile = function(u, par)
      stats::qnorm(u, mean = par["mean"], sd = par["sd"]),
    random = function(n, par)
      stats::rnorm(n, mean = par["mean"], sd = par["sd"]),
    scale = function(par) par["sd"], centered = FALSE,
    roles = c("location", "scale"), scale_is_sd = TRUE,
    metadata = list(variable_role = "conditioning",
                    scale_definition = "standard deviation"),
    support_fixed = TRUE)
}

copulaMarginCovariateLognormal <- function(meanlog, sdlog) {
  marginal_sd <- function(par) {
    s2 <- par["sdlog"]^2
    sqrt((exp(s2) - 1) * exp(2 * par["meanlog"] + s2))
  }
  copulaMargin("lognormal",
    c(meanlog = unname(meanlog), sdlog = unname(sdlog)),
    c(-50, 1e-8), c(50, 10),
    log_density = function(x, par)
      stats::dlnorm(x, par["meanlog"], par["sdlog"], log = TRUE),
    cdf = function(x, par)
      stats::plnorm(x, par["meanlog"], par["sdlog"]),
    quantile = function(u, par)
      stats::qlnorm(u, par["meanlog"], par["sdlog"]),
    random = function(n, par)
      stats::rlnorm(n, par["meanlog"], par["sdlog"]),
    scale = marginal_sd, centered = FALSE,
    roles = c("location", "scale"), scale_is_sd = TRUE,
    metadata = list(variable_role = "conditioning",
                    scale_definition = "standard deviation"),
    support_fixed = TRUE)
}

copulaMarginCovariateGamma <- function(shape, scale) {
  copulaMargin("gamma", c(shape = unname(shape), scale = unname(scale)),
    c(1e-4, 1e-8), c(1e6, 1e8),
    log_density = function(x, par)
      stats::dgamma(x, shape = par["shape"], scale = par["scale"], log = TRUE),
    cdf = function(x, par)
      stats::pgamma(x, shape = par["shape"], scale = par["scale"]),
    quantile = function(u, par)
      stats::qgamma(u, shape = par["shape"], scale = par["scale"]),
    random = function(n, par)
      stats::rgamma(n, shape = par["shape"], scale = par["scale"]),
    scale = function(par) sqrt(par["shape"]) * par["scale"],
    centered = FALSE, roles = c("shape", "scale"), scale_is_sd = TRUE,
    metadata = list(variable_role = "conditioning",
                    scale_definition = "standard deviation"),
    support_fixed = TRUE)
}

copulaMarginCovariateWeibull <- function(shape, scale) {
  marginal_sd <- function(par) {
    a <- gamma(1 + 2 / par["shape"])
    b <- gamma(1 + 1 / par["shape"])^2
    par["scale"] * sqrt(max(a - b, 0))
  }
  copulaMargin("weibull", c(shape = unname(shape), scale = unname(scale)),
    c(1e-4, 1e-8), c(1e4, 1e8),
    log_density = function(x, par)
      stats::dweibull(x, shape = par["shape"], scale = par["scale"], log = TRUE),
    cdf = function(x, par)
      stats::pweibull(x, shape = par["shape"], scale = par["scale"]),
    quantile = function(u, par)
      stats::qweibull(u, shape = par["shape"], scale = par["scale"]),
    random = function(n, par)
      stats::rweibull(n, shape = par["shape"], scale = par["scale"]),
    scale = marginal_sd, centered = FALSE,
    roles = c("shape", "scale"), scale_is_sd = TRUE,
    metadata = list(variable_role = "conditioning",
                    scale_definition = "standard deviation"),
    support_fixed = TRUE)
}

## Fit a small support-aware candidate set to one observed continuous
## covariate.  Discrete family selection is intentionally kept separate: an
## integer-valued variable may be a count, an ordinal score, or an identifier,
## which cannot be inferred safely from values alone.
copulaFitCovariateMargin <- function(x,
    candidates = c("auto", "normal", "lognormal", "gamma", "weibull")) {
  x <- as.numeric(x)
  if (any(is.infinite(x)))
    stop("copulaFitCovariateMargin does not accept infinite observations")
  x <- x[!is.na(x)]
  if (!length(x) || length(unique(x)) < 3L)
    stop("copulaFitCovariateMargin requires finite, non-degenerate observations")
  candidates <- unique(match.arg(candidates, several.ok = TRUE))
  if ("auto" %in% candidates)
    candidates <- if (all(x > 0))
      c("lognormal", "gamma", "weibull", "normal") else "normal"
  candidates <- unique(candidates)
  n <- length(x)
  normal <- copulaMarginCovariateNormal(mean(x),
    sqrt(mean((x - mean(x))^2)))
  fitted <- list(normal = normal)
  if (all(x > 0)) {
    lx <- log(x)
    fitted$lognormal <- copulaMarginCovariateLognormal(mean(lx),
      sqrt(mean((lx - mean(lx))^2)))
    gamma_start <- c(shape = mean(x)^2 / stats::var(x),
                     scale = stats::var(x) / mean(x))
    weibull_start <- c(shape = 2, scale = mean(x) / gamma(1 + 1 / 2))
    fit_positive <- function(start, density, lower, upper) {
      objective <- function(logpar) {
        par <- exp(logpar)
        value <- suppressWarnings(density(x, par))
        if (any(!is.finite(value))) return(.Machine$double.xmax / 100)
        -sum(value)
      }
      exploratory <- try(stats::optim(log(start), objective,
        method = "Nelder-Mead", control = list(maxit = 500L)), silent = TRUE)
      polishedStart <- if (!inherits(exploratory, "try-error") &&
          is.finite(exploratory$value)) exploratory$par else log(start)
      polishedStart <- pmin(log(upper), pmax(log(lower), polishedStart))
      ans <- stats::optim(polishedStart, objective, method = "L-BFGS-B",
        lower = log(lower), upper = log(upper), control = list(maxit = 300L))
      exp(ans$par)
    }
    gp <- fit_positive(gamma_start,
      function(y, p) stats::dgamma(y, shape = p[1], scale = p[2], log = TRUE),
      c(1e-4, 1e-8), c(1e6, 1e8))
    fitted$gamma <- copulaMarginCovariateGamma(gp[1], gp[2])
    wp <- fit_positive(weibull_start,
      function(y, p) stats::dweibull(y, shape = p[1], scale = p[2], log = TRUE),
      c(1e-4, 1e-8), c(1e4, 1e8))
    fitted$weibull <- copulaMarginCovariateWeibull(wp[1], wp[2])
  }
  fitted <- fitted[intersect(candidates, names(fitted))]
  if (!length(fitted)) stop("no requested covariate-margin family supports the data")
  score <- vapply(fitted, function(m) {
    ll <- sum(m$log_density(x, m$parameters))
    2 * sum(m$free) - 2 * ll
  }, numeric(1))
  winner <- fitted[[which.min(score)]]
  winner$metadata$selection <- list(criterion = "AIC", scores = score,
                                    n = n, candidates = names(fitted))
  copulaMarginValidate(winner)
  winner
}

copulaFitCovariateMargins <- function(x, candidates = "auto") {
  x <- as.matrix(x)
  out <- lapply(seq_len(ncol(x)), function(j)
    copulaFitCovariateMargin(x[, j], candidates = candidates))
  names(out) <- colnames(x) %||% paste0("conditioning", seq_len(ncol(x)))
  out
}

copulaMarginLayout <- function(margins) {
  if (!is.list(margins) || !length(margins)) stop("margins must be a non-empty list")
  lapply(margins, function(margin) if (
      inherits(margin, "saemix_natural_parameter_margin"))
      copulaNaturalMarginValidate(margin, probe = FALSE) else
      copulaMarginValidate(margin, probe = FALSE))
  index <- vector("list", length(margins)); par <- lower <- upper <- numeric()
  for (j in seq_along(margins)) {
    m <- margins[[j]]
    index[[j]] <- if (any(m$free)) length(par) + seq_len(sum(m$free)) else integer()
    par <- c(par, m$parameters[m$free])
    lower <- c(lower, m$lower[m$free]); upper <- c(upper, m$upper[m$free])
  }
  list(index = index, par = par, lower = lower, upper = upper)
}

copulaMarginsWithParameters <- function(margins, layout, parameters) {
  lapply(seq_along(margins), function(j) {
    margin <- margins[[j]]; index <- layout$index[[j]]
    native <- margin$parameters
    if (length(index)) native[margin$free] <- parameters[index]
    if (inherits(margin, "saemix_natural_parameter_margin"))
      copulaNaturalMarginWithParameters(margin, native) else
      copulaMarginWithParameters(margin, native)
  })
}

copulaMarginScales <- function(margins) vapply(margins, function(m) {
  z <- m$scale(m$parameters)
  if (length(z) != 1L || !is.finite(z) || z <= 0) NA_real_ else z
}, numeric(1))

copulaMarginsQuantile <- function(u, margins) {
  u <- as.matrix(u)
  if (ncol(u) != length(margins)) stop("quantile input has wrong dimension")
  out <- u
  for (j in seq_along(margins))
    out[, j] <- margins[[j]]$quantile(u[, j], margins[[j]]$parameters)
  out
}
