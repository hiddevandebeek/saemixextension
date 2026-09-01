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
  supportFixed <- margin$metadata$parameter_independent_support
  if (!is.null(supportFixed) &&
      (length(supportFixed) != 1L || is.na(supportFixed)))
    stop("invalid parameter-independent-support declaration")
  if (identical(supportFixed, FALSE) &&
      (!identical(margin$type, "continuous") ||
       !isTRUE(margin$metadata$fixed_reference_quantile)))
    stop("moving-support margins require a continuous fixed-reference quantile map")
  if (isTRUE(probe)) {
    q <- margin$quantile(c(.2, .5, .8), p)
    if (length(q) != 3L || any(!is.finite(q)))
      stop("margin quantile probe returned invalid values")
    Fq <- margin$cdf(q, p)
    ld <- margin$log_density(q, p)
    if (length(Fq) != 3L || any(!is.finite(Fq)) || any(Fq < 0 | Fq > 1) ||
        length(ld) != 3L || any(!is.finite(ld)))
      stop("margin density/CDF probe returned invalid values")
    if (identical(supportFixed, FALSE) &&
        max(abs(Fq - c(.2, .5, .8))) > 1e-6)
      stop("moving-support margin CDF/quantile round trip is not accurate")
  }
  invisible(TRUE)
}

copulaMarginHasMovingSupport <- function(margin) {
  identical(margin$metadata$parameter_independent_support, FALSE)
}

## Centered moving-support margins for latent eta coordinates.  Their means
## remain zero, so they do not duplicate the saemix population location.
copulaMarginUniformCentered <- function(sd = 1) {
  half <- sqrt(3) * unname(sd)
  copulaMargin("uniform-centered", c(half_width = half), 1e-6, 1e6,
    log_density = function(x, par) {
      h <- par[1L]; ifelse(x > -h & x < h, -log(2 * h), -Inf)
    },
    cdf = function(x, par) {
      h <- par[1L]; pmin(1, pmax(0, (x + h) / (2 * h)))
    },
    quantile = function(u, par) par[1L] * (2 * u - 1),
    random = function(n, par) stats::runif(n, -par[1L], par[1L]),
    scale = function(par) par[1L] / sqrt(3), centered = TRUE,
    set_scale = function(par, value) { par[1L] <- sqrt(3) * value; par },
    roles = "scale", scale_is_sd = TRUE,
    metadata = list(scale_definition = "standard deviation"),
    support_fixed = FALSE)
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

copulaMarginCovariateUniform <- function(lower, upper) {
  if (!is.finite(lower) || !is.finite(upper) || lower >= upper)
    stop("uniform covariate bounds must be finite and increasing")
  copulaMargin("uniform", c(lower = unname(lower), log_width = log(upper - lower)),
    c(-1e12, log(1e-8)), c(1e12, log(1e12)),
    log_density = function(x, par) {
      lo <- par["lower"]; hi <- lo + exp(par["log_width"])
      ifelse(x > lo & x < hi, -par["log_width"], -Inf)
    },
    cdf = function(x, par) {
      lo <- par["lower"]; pmin(1, pmax(0,
        (x - lo) / exp(par["log_width"])))
    },
    quantile = function(u, par)
      par["lower"] + exp(par["log_width"]) * u,
    random = function(n, par)
      par["lower"] + exp(par["log_width"]) * stats::runif(n),
    scale = function(par) exp(par["log_width"]) / sqrt(12),
    centered = FALSE, roles = c("location", "scale"), scale_is_sd = TRUE,
    metadata = list(variable_role = "conditioning",
      scale_definition = "standard deviation"), support_fixed = FALSE)
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

copulaMarginWithFreeParameters <- function(margin, parameters) {
  if (length(parameters) != sum(margin$free))
    stop("wrong number of free marginal parameters")
  p <- margin$parameters
  p[margin$free] <- parameters
  copulaMarginWithParameters(margin, p)
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

## Proper nominal categorical margin with an explicitly declared latent order.
## A scalar copula CDF cannot be permutation invariant for unordered labels;
## requiring latent_order prevents an accidental alphabetical/numeric order
## from silently becoming part of the dependence model.
copulaMarginCategorical <- function(probabilities, labels,
                                    latent_order, free = TRUE) {
  if (missing(latent_order))
    stop("nominal categorical margins require an explicit permutation latent_order")
  labels <- as.character(labels); latent_order <- as.character(latent_order)
  if (length(latent_order) != length(labels) ||
      !setequal(latent_order, labels))
    stop("nominal categorical margins require an explicit permutation latent_order")
  index <- match(latent_order, labels)
  margin <- copulaMarginOrdinal(as.numeric(probabilities)[index],
    labels = latent_order, free = free, name = "categorical")
  margin$metadata$categorical_kind <- "nominal-with-declared-latent-order"
  margin$metadata$latent_order <- latent_order
  margin
}

## Adapter for any distribution exposing compatible d/p/q functions.  The
## caller explicitly chooses native parameters, bounds, and centring; this
## prevents unsafe guessing from the heterogeneous stats::Distributions APIs.
copulaMarginDistribution <- function(distr, parameters, lower, upper,
                                     free = rep(TRUE, length(parameters)),
                                     type = c("continuous", "discrete"),
                                     centered = FALSE, scale = NULL,
                                     set_scale = NULL, cdf_left = NULL,
                                     roles = NULL, scale_is_sd = FALSE,
                                     support_fixed = NULL) {
  type <- match.arg(type)
  d <- get(paste0("d", distr), mode = "function")
  p <- get(paste0("p", distr), mode = "function")
  q <- get(paste0("q", distr), mode = "function")
  r <- get0(paste0("r", distr), mode = "function")
  call_with <- function(fun, x, par, extra = list())
    do.call(fun, c(list(x), as.list(stats::setNames(par, names(parameters))), extra))
  copulaMargin(distr, parameters, lower, upper, free = free,
    log_density = function(x, par) call_with(d, x, par, list(log = TRUE)),
    cdf = function(x, par) call_with(p, x, par),
    quantile = function(u, par) call_with(q, u, par),
    random = if (is.null(r)) NULL else function(n, par)
      do.call(r, c(list(n), as.list(stats::setNames(par, names(parameters))))),
    type = type, cdf_left = cdf_left, scale = scale,
    set_scale = set_scale,
    centered = centered, roles = roles, scale_is_sd = scale_is_sd,
    metadata = list(adapter = "stats-distribution"),
    support_fixed = support_fixed)
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

## Pre-index columns that are unchanged throughout an optimizer call.  The
## values of observed conditioning covariates repeat once for every retained
## latent particle; their integer expansion map is therefore fixed even while
## the covariate-margin parameters change.  This stores no likelihood value,
## only the exact unique values and row map.
copulaMarginsCompression <- function(x, columns = seq_len(ncol(as.matrix(x)))) {
  x <- as.matrix(x); out <- vector("list", ncol(x))
  for (j in intersect(seq_len(ncol(x)), as.integer(columns))) {
    values <- x[, j]
    if (nrow(x) >= 4L && anyDuplicated(values)) {
      uniqueValues <- unique(values)
      if (length(uniqueValues) * 2L <= nrow(x))
        out[[j]] <- list(values = uniqueValues,
                         map = match(values, uniqueValues))
    }
  }
  out
}

copulaMarginsEvaluate <- function(x, margins, eps = NULL,
                                  numericalPolicy = c("exact", "clip"),
                                  compression = NULL) {
  numericalPolicy <- match.arg(numericalPolicy)
  if (identical(numericalPolicy, "exact") && !is.null(eps))
    stop("eps is only valid with numericalPolicy='clip'")
  x <- as.matrix(x)
  if (ncol(x) != length(margins))
    stop("margin evaluation dimension does not match the eta matrix")
  n <- nrow(x); d <- ncol(x)
  u <- matrix(NA_real_, n, d); log_margin <- matrix(NA_real_, n, d)
  discrete <- which(vapply(margins, function(m) m$type == "discrete", logical(1)))
  u_left <- if (length(discrete)) matrix(NA_real_, n, length(discrete)) else NULL
  for (j in seq_len(d)) {
    m <- margins[[j]]; p <- m$parameters
    ## Conditioning values repeat once per retained latent particle. Evaluate a
    ## repeated column on its exact unique values and expand by integer index;
    ## this is an algebraic compression of the same marginal likelihood terms.
    ## The 2:1 threshold avoids paying unique()/match() overhead for latent
    ## columns whose draws are effectively all distinct.
    values <- x[, j]
    cached <- if (length(compression) >= j) compression[[j]] else NULL
    compressed <- !is.null(cached)
    if (compressed) {
      uniqueValues <- cached$values; map <- cached$map
      if (length(map) != n)
        stop("margin compression map has the wrong number of rows")
    } else if (n >= 4L && anyDuplicated(values)) {
      uniqueValues <- unique(values)
      compressed <- length(uniqueValues) * 2L <= n
      if (compressed) map <- match(values, uniqueValues)
    }
    if (compressed) {
      log_margin[, j] <- m$log_density(uniqueValues, p)[map]
      u[, j] <- m$cdf(uniqueValues, p)[map]
    } else {
      log_margin[, j] <- m$log_density(values, p)
      u[, j] <- m$cdf(values, p)
    }
    if (j %in% discrete) {
      jj <- match(j, discrete)
      u_left[, jj] <- if (compressed)
        m$cdf_left(uniqueValues, p)[map] else m$cdf_left(values, p)
    }
  }
  if (any(is.na(log_margin)) || any(is.nan(log_margin)) ||
      any(log_margin == Inf) || any(!is.finite(u)))
    stop("marginal density/mass or CDF returned an undefined value")
  if (any(u < 0 | u > 1)) stop("a marginal CDF returned a value outside [0,1]")
  zero_density <- apply(log_margin == -Inf, 1L, any)
  rawRange <- range(u)
  continuous <- setdiff(seq_len(d), discrete)
  clipped <- matrix(FALSE, n, d)
  if (identical(numericalPolicy, "exact")) {
    interiorRows <- !zero_density
    if (length(continuous) && any((u[, continuous, drop = FALSE] <= 0 |
                                  u[, continuous, drop = FALSE] >= 1) &
                                 interiorRows))
      stop(paste0("an exact continuous marginal PIT reached 0 or 1 in floating-point ",
                  "arithmetic; the likelihood was not clipped"))
  } else {
    if (is.null(eps)) eps <- 1e-10
    if (length(eps) != 1L || !is.finite(eps) || eps <= 0 || eps >= 0.5)
      stop("eps must be in (0, 0.5) for numericalPolicy='clip'")
    if (length(continuous)) {
      clipped[, continuous] <- u[, continuous, drop = FALSE] <= eps |
        u[, continuous, drop = FALSE] >= 1 - eps
      u[, continuous] <- pmin(pmax(u[, continuous, drop = FALSE], eps), 1 - eps)
    }
  }
  if (length(discrete)) {
    if (any(!is.finite(u_left))) stop("marginal left-limit CDF returned non-finite values")
    if (any(u_left < 0 | u_left > 1))
      stop("a discrete marginal left-limit CDF returned a value outside [0,1]")
    if (any(u_left > u[, discrete, drop = FALSE] + 1e-14))
      stop("a discrete marginal left-limit CDF exceeds its CDF")
  }
  list(log_margin = rowSums(log_margin), u = u, u_left = u_left,
       vine_u = if (length(discrete)) cbind(u, u_left) else u,
       discrete = discrete, clipped = clipped, n_clipped = sum(clipped),
       zero_density = zero_density, raw_range = rawRange,
       evaluated_range = range(u))
}

copulaMarginsLogDensity <- function(x,margins) {
  x<-as.matrix(x)
  if(ncol(x)!=length(margins))stop("margin density dimension mismatch")
  z<-vapply(seq_along(margins),function(j)
    margins[[j]]$log_density(x[,j],margins[[j]]$parameters),numeric(nrow(x)))
  z<-matrix(z,nrow=nrow(x),ncol=length(margins))
  if(any(is.na(z))||any(is.nan(z))||any(z==Inf))
    stop("marginal density/mass returned an undefined value")
  rowSums(z)
}

copulaLogDensity <- function(u, vine, cores = 1L,
                             numericalPolicy = c("exact", "floor"),
                             floor = 1e-300) {
  numericalPolicy <- match.arg(numericalPolicy)
  density <- rvinecopulib::dvinecop(as.matrix(u), vine,
                                    cores = as.integer(cores))
  invalid <- !is.finite(density) | density <= 0
  if (identical(numericalPolicy, "exact")) {
    if (any(invalid))
      stop(paste0("the exact vine density was non-positive or non-finite in ",
                  "floating-point arithmetic; the likelihood was not floored"))
    out <- log(density)
  } else {
    if (length(floor) != 1L || !is.finite(floor) || floor <= 0)
      stop("floor must be finite and positive")
    if (any(is.na(density)) || any(is.nan(density)) || any(density < 0) ||
        any(density == Inf))
      stop("the vine density returned an undefined value; a lower floor was not applied")
    out <- log(pmax(density, floor))
  }
  attr(out, "density_floored") <- if (identical(numericalPolicy, "floor"))
    sum(density <= floor | !is.finite(density)) else 0L
  out
}

copulaMarginsQuantile <- function(u, margins) {
  u <- as.matrix(u)
  if (ncol(u) != length(margins)) stop("quantile input has wrong dimension")
  out <- u
  for (j in seq_along(margins))
    out[, j] <- margins[[j]]$quantile(u[, j], margins[[j]]$parameters)
  out
}
