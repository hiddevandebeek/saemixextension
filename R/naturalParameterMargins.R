## Natural-scale individual-parameter margins. The distribution is independent
## of the saemix working-coordinate transform; mapping and Jacobians are applied
## separately by copulaNaturalToWorking()/copulaWorkingToNatural().

## `score` and `inverse_score` are optional closed forms for the Gaussian score
## xi = Phi^-1{F(x)} and its inverse. They default to the probability integral
## transform written out with pnorm/qnorm. Families whose score is affine in a
## monotone transform of x -- Normal, lognormal, logit-normal -- can supply the
## composition directly: the two special functions cancel exactly, so the hook
## is faster and does not lose the far tails to underflow in the uniform.

copulaNaturalMargin <- function(name, support, anchor, parameters, lower, upper,
    log_density, cdf, quantile, random = NULL, score = NULL,
    inverse_score = NULL, evaluate = NULL,
    free = rep(TRUE, length(parameters)), roles = NULL, metadata = list()) {
  support <- match.arg(support, c("real", "positive", "unit"))
  anchor <- match.arg(anchor, c("mean", "median", "geometric-mean", "custom"))
  parameters <- as.numeric(parameters) |>
    stats::setNames(names(parameters) %||% paste0("par", seq_along(parameters)))
  if (is.null(roles)) roles <- rep("shape", length(parameters))
  closedScore <- !is.null(score) && !is.null(inverse_score)
  if (is.null(score)) score <- function(x, typical, par)
    stats::qnorm(cdf(x, typical, par))
  if (is.null(inverse_score)) inverse_score <- function(z, typical, par)
    quantile(stats::pnorm(z), typical, par)
  if (is.null(evaluate)) evaluate <- function(x, typical, par)
    list(log_density = log_density(x, typical, par),
      score = score(x, typical, par))
  margin <- structure(list(name = as.character(name), support = support,
    anchor = anchor, type = "continuous", parameters = parameters,
    lower = as.numeric(lower), upper = as.numeric(upper),
    free = as.logical(free), roles = as.character(roles),
    log_density = log_density, cdf = cdf, quantile = quantile,
    random = random, score = score, inverse_score = inverse_score,
    evaluate = evaluate, closed_score = closedScore, metadata = metadata),
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
    ## A supplied closed form must agree with the transform it replaces.
    combined <- margin$evaluate(q, typical, p)
    if (max(abs(combined$log_density - margin$log_density(q, typical, p))) >
          1e-8 ||
        max(abs(combined$score - margin$score(q, typical, p))) > 1e-8)
      stop("natural parameter margin combined evaluation is inconsistent")
    if (isTRUE(margin$closed_score)) {
      zeta <- stats::qnorm(u)
      if (max(abs(margin$score(q, typical, p) - zeta)) > 1e-8 ||
          max(abs(margin$inverse_score(zeta, typical, p) - q)) >
            1e-8 * max(1, max(abs(q))))
        stop("natural parameter margin closed-form score disagrees with its CDF")
    }
  }
  invisible(TRUE)
}

## Freeze parameters at their current values, so the recursion does not
## estimate them.
##
## Supplying a margin object to a specification fixes the *family* and gives a
## starting value; every parameter stays free and is estimated like any other.
## That is usually what is wanted. When it is not -- a scale known from an
## external study, a shape fixed by regulatory precedent -- this pins the named
## parameters, or all of them if none are named.
copulaFixMargin <- function(margin, which = NULL) {
  if (is.null(margin$free))
    stop("copulaFixMargin expects a margin object")
  named <- names(margin$parameters)
  if (is.null(which)) {
    margin$free <- rep(FALSE, length(margin$parameters))
  } else {
    which <- as.character(which)
    unknown <- setdiff(which, named)
    if (length(unknown))
      stop("this margin has no parameter called ",
        paste(unknown, collapse = ", "), ". It has: ",
        paste(named, collapse = ", "))
    margin$free[named %in% which] <- FALSE
  }
  if (!any(margin$free))
    margin$metadata$fullyFixed <- TRUE
  margin
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
    random = function(n, typical, par) stats::rnorm(n, typical, par["sd"]),
    score = function(x, typical, par) (x - typical) / par["sd"],
    inverse_score = function(z, typical, par) typical + par["sd"] * z,
    evaluate = function(x, typical, par) {
      zeta <- (x - typical) / par["sd"]
      list(log_density = -log(par["sd"]) - .5 * zeta^2 - .5 * log(2 * pi),
        score = zeta)
    })
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
      stats::rlnorm(n, meanlog(typical, par), par["sdlog"]),
    score = function(x, typical, par)
      (log(x) - meanlog(typical, par)) / par["sdlog"],
    inverse_score = function(z, typical, par)
      exp(meanlog(typical, par) + par["sdlog"] * z),
    evaluate = function(x, typical, par) {
      logx <- log(x)
      zeta <- (logx - meanlog(typical, par)) / par["sdlog"]
      list(log_density = -logx - log(par["sdlog"]) - .5 * zeta^2 -
        .5 * log(2 * pi), score = zeta)
    })
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

copulaNaturalMarginLogitNormal <- function(sdlogit = .5,
    anchor = c("median", "custom")) {
  anchor <- match.arg(anchor)
  location <- function(typical) stats::qlogis(typical)
  copulaNaturalMargin("logitnormal", "unit", anchor, c(sdlogit = sdlogit),
    1e-8, 1e4, roles = "scale",
    log_density = function(x, typical, par)
      stats::dnorm(stats::qlogis(x), location(typical), par["sdlogit"],
        log = TRUE) - log(x) - log1p(-x),
    cdf = function(x, typical, par)
      stats::pnorm(stats::qlogis(x), location(typical), par["sdlogit"]),
    quantile = function(u, typical, par)
      stats::plogis(stats::qnorm(u, location(typical), par["sdlogit"])),
    random = function(n, typical, par)
      stats::plogis(stats::rnorm(n, location(typical), par["sdlogit"])),
    score = function(x, typical, par)
      (stats::qlogis(x) - location(typical)) / par["sdlogit"],
    inverse_score = function(z, typical, par)
      stats::plogis(location(typical) + par["sdlogit"] * z),
    evaluate = function(x, typical, par) {
      zeta <- (stats::qlogis(x) - location(typical)) / par["sdlogit"]
      list(log_density = -log(x) - log1p(-x) - log(par["sdlogit"]) -
        .5 * zeta^2 - .5 * log(2 * pi), score = zeta)
    })
}

copulaNaturalMarginBeta <- function(precision = 10,
    anchor = c("mean", "median")) {
  anchor <- match.arg(anchor)
  ## The median anchor has no closed form, so it is solved numerically. Solve
  ## only once per distinct typical value: a parameter without covariate
  ## effects has a single one, and the hot loop then costs one root solve.
  shapes <- function(typical, par) {
    kappa <- par["precision"]
    if (anchor == "mean") return(cbind(typical * kappa, (1 - typical) * kappa))
    key <- unique(typical)
    solved <- vapply(key, function(target) stats::uniroot(function(m)
      stats::qbeta(.5, m * kappa, (1 - m) * kappa) - target,
      c(1e-6, 1 - 1e-6), tol = 1e-12)$root, numeric(1))
    mean <- solved[match(typical, key)]
    cbind(mean * kappa, (1 - mean) * kappa)
  }
  copulaNaturalMargin("beta", "unit", anchor, c(precision = precision),
    .05, 1e6, roles = "shape",
    log_density = function(x, typical, par) {
      ab <- shapes(typical, par)
      stats::dbeta(x, ab[, 1L], ab[, 2L], log = TRUE)
    },
    cdf = function(x, typical, par) {
      ab <- shapes(typical, par)
      stats::pbeta(x, ab[, 1L], ab[, 2L])
    },
    quantile = function(u, typical, par) {
      ab <- shapes(typical, par)
      stats::qbeta(u, ab[, 1L], ab[, 2L])
    },
    random = function(n, typical, par) {
      ab <- shapes(typical, par)
      stats::rbeta(n, ab[, 1L], ab[, 2L])
    })
}

## Generalized gamma on the positive line in Prentice's parameterisation, with
## free shape parameters (sigma, Q) and the location fixed by the anchor:
##   w = (log x - mu) / sigma,   u = Q^-2 exp(Q w) ~ Gamma(Q^-2, 1).
## The three positive families the registry ships are points of this family --
## Q = 1 is Weibull, Q = sigma is Gamma, and Q -> 0 is lognormal -- all at
## finite parameter values, so a discrete choice between them can instead be
## estimated continuously. Stacy's (shape, power) parameterisation puts
## lognormal at shape -> Inf, where an optimiser simply runs away; Prentice's
## puts it at a bounded edge of a compact box, which is why it is used here.
## The lognormal edge is still an edge: near it sigma and Q are only weakly
## separated, so read the pair together rather than individually.

## Prentice's parameterisation, in which the shape Q runs over the whole real
## line and the lognormal sits at Q = 0 as an INTERIOR point. That interiority
## is the reason for this form rather than Stacy's: it makes "is this margin
## lognormal?" an ordinary hypothesis about one coordinate of theta, testable
## with the standard error the score recursion already produces, instead of a
## test on the boundary of the parameter space with non-standard asymptotics.
##
## The sign of Q is the direction of the departure: Q > 0 gives log-scale left
## skew, Q < 0 a heavier upper tail on the log scale. The named special cases
## sit where X = exp(mu) (Q^2 g)^(sigma/Q), with g ~ Gamma(Q^-2), collapses --
## the gamma when the exponent sigma/Q equals one, that is sigma = Q, and then
## the shape is Q^-2; the Weibull when Q^-2 = 1. The gamma is not at Q = 1 for
## a general sigma, which is easy to assume and wrong. Restricting to Q > 0, as
## an earlier version did, made half of the possible departures unreachable and
## -- worse -- silently so, since a truth on the wrong side simply pinned the
## estimate at the lower bound and reported something indistinguishable from a
## lognormal.
##
## Writing w = (log x - location) / sigma and k = Q^-2, the standard density
## follows from g = k exp(Q w) being Gamma(k, 1): the Jacobian is |Q| g, so the
## |Q| below is not a typo for Q but what makes the density integrate to one on
## both sides of zero. The CDF and quantile invert in the usual way, with the
## tail reflected when Q < 0 because w is then decreasing in g.
##
## Near Q = 0 the exact expressions are unusable: k diverges, and the density
## subtracts two quantities of size k to leave a remainder of order one. The
## series below is the expansion of the same density,
##
##   log f(w) = -log(2 pi)/2 - w^2/2 - Q w^3/6 - Q^2 (w^4/24 + 1/12) + O(Q^3),
##
## which is used inside the switch radius. It matters that this is a series and
## not simply the lognormal: the leading correction -w^3/6 is exactly the score
## for Q at Q = 0, so the derivative the whole test rests on stays right at the
## null. Substituting a flat lognormal branch would zero that score and quietly
## destroy the test it exists to support.
copulaNaturalMarginGeneralizedGamma <- function(sigma = .3, Q = 1,
    anchor = c("geometric-mean", "median")) {
  anchor <- match.arg(anchor)
  ## Radius inside which the series replaces the exact form. At |Q| = 1e-3 the
  ## exact expression still carries about ten digits and the series truncation
  ## is of order Q^3, so the two agree far beyond what any optimiser resolves.
  small <- 1e-3
  ## E[w] under the standard density, used to pin the geometric mean.
  ## digamma(k) - log(k) -> 0 as k -> Inf, so the ratio is a 0/0 at Q = 0 and
  ## needs its own expansion: -Q/2 - Q^3/12.
  meanW <- function(par) {
    q <- par[["Q"]]
    if (abs(q) < small) return(-q / 2 - q^3 / 12)
    (2 * log(abs(q)) + digamma(q^-2)) / q
  }
  medianW <- function(par) {
    q <- par[["Q"]]
    if (abs(q) < small) return(-q / 2)
    (2 * log(abs(q)) + log(stats::qgamma(.5, q^-2))) / q
  }
  location <- function(typical, par) switch(anchor,
    `geometric-mean` = log(typical) - par[["sigma"]] * meanW(par),
    median = log(typical) - par[["sigma"]] * medianW(par))
  copulaNaturalMargin("generalizedgamma", "positive", anchor,
    c(sigma = sigma, Q = Q), c(1e-2, -3), c(5, 3),
    roles = c("scale", "shape"),
    log_density = function(x, typical, par) {
      q <- par[["Q"]]; sg <- par[["sigma"]]
      w <- (log(x) - location(typical, par)) / sg
      if (abs(q) < small)
        return(-.5 * log(2 * pi) - .5 * w^2 - q * w^3 / 6 -
          q^2 * (w^4 / 24 + 1 / 12) - log(sg) - log(x))
      k <- q^-2
      log(abs(q)) - log(sg) - log(x) - lgamma(k) + k * log(k) +
        k * (q * w - exp(q * w))
    },
    cdf = function(x, typical, par) {
      q <- par[["Q"]]
      w <- (log(x) - location(typical, par)) / par[["sigma"]]
      if (abs(q) < small) return(stats::pnorm(w))
      k <- q^-2
      p <- stats::pgamma(k * exp(q * w), k)
      if (q > 0) p else 1 - p
    },
    quantile = function(u, typical, par) {
      q <- par[["Q"]]
      if (abs(q) < small)
        return(exp(location(typical, par) + par[["sigma"]] * stats::qnorm(u)))
      k <- q^-2
      w <- log(stats::qgamma(if (q > 0) u else 1 - u, k) / k) / q
      exp(location(typical, par) + par[["sigma"]] * w)
    },
    random = function(n, typical, par) {
      q <- par[["Q"]]
      if (abs(q) < small)
        return(exp(location(typical, par) + par[["sigma"]] * stats::rnorm(n)))
      k <- q^-2
      w <- log(stats::rgamma(n, k) / k) / q
      exp(location(typical, par) + par[["sigma"]] * w)
    },
    ## The Gaussian score and its inverse in closed form, computed through the
    ## log of whichever tail is the smaller.
    ##
    ## Without these the generic fallback composes quantile(pnorm(z)), and that
    ## round trip destroys the far tail: the reporting quadrature evaluates at
    ## Gauss-Hermite nodes out to |z| = 17, where pnorm saturates to exactly 1
    ## in double precision and qgamma(1, k) is infinite. The margin is then
    ## perfectly usable for fitting and yet every derived summary -- the median,
    ## the spread of log psi, the coefficient of variation, which are the
    ## numbers actually reported for a random effect -- comes back NA. The
    ## quantile itself is not extreme (for k = 4 the node maps to about 160);
    ## only the probability is, so carrying it as a log tail is enough.
    score = function(x, typical, par) {
      q <- par[["Q"]]
      w <- (log(x) - location(typical, par)) / par[["sigma"]]
      if (abs(q) < small) return(w)
      k <- q^-2; t <- k * exp(q * w)
      lower <- stats::pgamma(t, k, log.p = TRUE)
      upper <- stats::pgamma(t, k, lower.tail = FALSE, log.p = TRUE)
      logF <- if (q > 0) lower else upper
      logS <- if (q > 0) upper else lower
      ifelse(logF <= log(.5), stats::qnorm(logF, log.p = TRUE),
        -stats::qnorm(logS, log.p = TRUE))
    },
    inverse_score = function(z, typical, par) {
      q <- par[["Q"]]
      if (abs(q) < small)
        return(exp(location(typical, par) + par[["sigma"]] * z))
      k <- q^-2
      ## Q < 0 reflects the CDF, so it reflects the normal score with it.
      zz <- if (q > 0) z else -z
      g <- ifelse(zz <= 0,
        stats::qgamma(stats::pnorm(zz, log.p = TRUE), k, log.p = TRUE),
        stats::qgamma(stats::pnorm(zz, lower.tail = FALSE, log.p = TRUE), k,
          lower.tail = FALSE, log.p = TRUE))
      exp(location(typical, par) + par[["sigma"]] * log(g / k) / q)
    })
}

## Two-component lognormal mixture on the positive line: the smallest margin
## that can be bimodal, which no single-shape family can represent and which is
## the usual reason a pharmacokinetic parameter is not lognormal. Free
## parameters are the second-component weight, the log separation of the two
## component medians, and a shared log-scale spread. Component locations are
## placed so that E[log X] equals the anchor for every parameter value, which
## keeps the population location identified:
##   mu1 = m - w * delta,   mu2 = m + (1 - w) * delta.
## Setting the weight or the separation to zero returns a plain lognormal, so
## the family is a strict extension of the shipped one and the extra parameters
## are exactly what a likelihood criterion has to pay for.
##
## The quantile function has no closed form and is obtained by bisection on the
## log scale, so this margin is slower than the analytic families.

copulaNaturalMarginLognormalMixture <- function(weight = .3, separation = .8,
    sdlog = .25, anchor = c("geometric-mean")) {
  anchor <- match.arg(anchor)
  parts <- function(typical, par) {
    w <- par[["weight"]]; d <- par[["separation"]]; m <- log(typical)
    list(w = w, mu1 = m - w * d, mu2 = m + (1 - w) * d, s = par[["sdlog"]])
  }
  cdfOf <- function(x, typical, par) {
    p <- parts(typical, par)
    (1 - p$w) * stats::plnorm(x, p$mu1, p$s) +
      p$w * stats::plnorm(x, p$mu2, p$s)
  }
  log_density_mixture <- function(x, typical, par) {
    p <- parts(typical, par)
    a <- log1p(-p$w) + stats::dlnorm(x, p$mu1, p$s, log = TRUE)
    b <- log(p$w) + stats::dlnorm(x, p$mu2, p$s, log = TRUE)
    hi <- pmax(a, b)
    hi + log(exp(a - hi) + exp(b - hi))
  }
  copulaNaturalMargin("lognormalmixture", "positive", anchor,
    c(weight = weight, separation = separation, sdlog = sdlog),
    c(1e-3, 0, 1e-2), c(.999, 5, 3), roles = c("shape", "shape", "scale"),
    log_density = log_density_mixture,
    cdf = cdfOf,
    quantile = function(u, typical, par) {
      ## Both the CDF and its derivative are closed form, so the inverse is
      ## found by safeguarded Newton on the log scale rather than by bisection:
      ## quadratic convergence reaches machine precision in a handful of passes
      ## where bisection needs dozens. The bracket is kept and any Newton step
      ## that leaves it falls back to a bisection step, so the iteration cannot
      ## diverge on the flat region between two well-separated components.
      ## The mixture CDF is a convex combination of the two component CDFs, so
      ## its quantile always lies between the two component quantiles: at the
      ## smaller one the mixture CDF is at most u, at the larger one at least u.
      ## That bracket has width equal to the component separation, independent
      ## of u, and is far tighter than any fixed multiple of the spread. Newton
      ## then does the work, falling back to bisection in the low-density
      ## region between well-separated modes where the slope underflows.
      p <- parts(typical, par)
      z <- stats::qnorm(u)
      q1 <- p$mu1 + p$s * z; q2 <- p$mu2 + p$s * z
      lo <- pmin(q1, q2); hi <- pmax(q1, q2)
      ## Starting point. When the components are well separated the first one
      ## carries the mass below u = 1 - w and the second carries the rest, so
      ## inverting the responsible component alone is almost exact -- and that
      ## is precisely the configuration where a blend start is worst and the
      ## density between the modes is too flat for Newton to recover quickly.
      eps <- 1e-12
      first <- p$mu1 + p$s * stats::qnorm(pmin(u / max(1 - p$w, eps), 1 - eps))
      second <- p$mu2 + p$s * stats::qnorm(pmax(
        (u - (1 - p$w)) / max(p$w, eps), eps))
      t <- pmin(pmax(ifelse(u < 1 - p$w, first, second), lo), hi)
      ## Most points converge in three or four Newton passes; only those in the
      ## low-density region between well-separated modes need the bisection
      ## fallback, so the loop stops as soon as every point is converged.
      for (i in seq_len(24L)) {
        value <- cdfOf(exp(t), typical, par) - u
        if (max(abs(value)) < 1e-13) break
        below <- value < 0
        lo <- ifelse(below, t, lo); hi <- ifelse(below, hi, t)
        slope <- exp(log_density_mixture(exp(t), typical, par) + t)
        step <- ifelse(slope > 0, value / slope, 0)
        proposal <- t - step
        ## strict, so a converged step landing on a bracket endpoint is kept
        outside <- !is.finite(proposal) | proposal < lo | proposal > hi
        t <- ifelse(outside, (lo + hi) / 2, proposal)
      }
      exp(t)
    },
    random = function(n, typical, par) {
      p <- parts(typical, par)
      second <- stats::runif(n) < p$w
      stats::rlnorm(n, ifelse(second, p$mu2, p$mu1), p$s)
    })
}

## Semi-nonparametric (SNP) marginal density, after Gallant and Nychka (1987)
## and its use for nonlinear mixed effects by Davidian and Gallant (1993).
##
## On a standardised latent scale the density is a squared polynomial times a
## standard normal,
##
##   f(z) = P(z)^2 phi(z) / (a' M a),   P(z) = sum_{k=0}^{K} a_k z^k,
##
## with a_0 fixed at 1 for identifiability and M_{kl} = E[Z^{k+l}] the standard
## normal moments. Squaring makes the density non-negative without constraints,
## and the normaliser is a quadratic form in the coefficients, so everything is
## closed form -- including the CDF, because the partial moments of the normal
## satisfy I_r(z) = -z^{r-1} phi(z) + (r-1) I_{r-2}(z).
##
## Two properties matter here. Degree zero is exactly the base family --
## lognormal on the positive line, Normal on the real line -- so the standard
## model is nested at K = 0 and the degree is the only discrete choice, one
## number per coordinate rather than a family per coordinate. And the score in
## the coefficients is analytic, so this needs no finite differences: for an
## exponential-family-like density the observed-data score reduces by Fisher's
## identity to a pooled-posterior moment minus a model moment, which is exactly
## the stationarity condition the individual conditional distributions define.

copulaSnpMoments <- function(order) {
  m <- numeric(order + 1L); m[1L] <- 1
  if (order >= 2L) for (r in 2:order)
    m[r + 1L] <- if (r %% 2L == 1L) 0 else (r - 1L) * m[r - 1L]
  m
}

## int_{-inf}^{z} t^r phi(t) dt, for r = 0 .. order
copulaSnpPartialMoments <- function(z, order) {
  out <- matrix(0, length(z), order + 1L)
  out[, 1L] <- stats::pnorm(z)
  if (order < 1L) return(out)
  ## the normal density and the running power of z are shared by every order,
  ## so both are carried through the recursion rather than recomputed in it
  density <- stats::dnorm(z)
  out[, 2L] <- -density
  if (order >= 2L) {
    power <- rep(1, length(z))
    for (r in 2:order) {
      power <- power * z
      out[, r + 1L] <- (r - 1L) * out[, r - 1L] - power * density
    }
  }
  out
}

copulaNaturalMarginSNP <- function(degree = 2L, scale = .3,
    coefficients = NULL, support = c("positive", "real")) {
  support <- match.arg(support)
  degree <- as.integer(degree)
  if (degree < 0L || degree > 6L) stop("SNP degree must be between 0 and 6")
  if (is.null(coefficients)) coefficients <- rep(0, degree)
  if (length(coefficients) != degree) stop("coefficients must have length degree")
  moments <- copulaSnpMoments(2L * degree + 1L)
  full <- function(par) c(1, if (degree) par[paste0("a", seq_len(degree))] else
    numeric())
  ## coefficients of P(z)^2, so every quantity below is a single pass over
  ## 2K+1 terms rather than a (K+1)^2 double sum
  square <- function(a) if (length(a) == 1L) a^2 else
    stats::convolve(a, rev(a), type = "open")
  horner <- function(a, z) {
    out <- rep(a[length(a)], length(z))
    if (length(a) > 1L) for (k in rev(seq_len(length(a) - 1L)))
      out <- out * z + a[k]
    out
  }
  quad <- function(a, shift = 0L) {
    c2 <- square(a)
    sum(c2 * moments[seq_along(c2) + shift])
  }
  centre <- function(par) quad(full(par), 1L) / quad(full(par))
  location <- function(typical, par) {
    base <- if (identical(support, "positive")) log(typical) else typical
    base - par[["scale"]] * centre(par)
  }
  standardise <- function(x, typical, par) {
    inner <- if (identical(support, "positive")) log(x) else x
    (inner - location(typical, par)) / par[["scale"]]
  }
  parameters <- c(scale = scale)
  if (degree) parameters <- c(parameters,
    stats::setNames(coefficients, paste0("a", seq_len(degree))))
  copulaNaturalMargin("snp", support,
    if (identical(support, "positive")) "geometric-mean" else "mean",
    parameters,
    c(1e-3, rep(-25, degree)), c(10, rep(25, degree)),
    roles = c("scale", rep("shape", degree)),
    log_density = function(x, typical, par) {
      z <- standardise(x, typical, par); a <- full(par)
      poly <- horner(a, z)
      value <- 2 * log(abs(poly)) + stats::dnorm(z, log = TRUE) -
        log(quad(a)) - log(par[["scale"]])
      if (identical(support, "positive")) value - log(x) else value
    },
    cdf = function(x, typical, par) {
      z <- standardise(x, typical, par); a <- full(par)
      c2 <- square(a)
      partial <- copulaSnpPartialMoments(z, 2L * degree)
      pmin(pmax(as.numeric(partial[, seq_along(c2), drop = FALSE] %*% c2) /
        quad(a), 0), 1)
    },
    quantile = function(u, typical, par) {
      ## closed-form CDF and density, so safeguarded Newton on the latent scale
      a <- full(par); c2 <- square(a); norm <- quad(a)
      lo <- rep(-12, length(u)); hi <- rep(12, length(u))
      z <- stats::qnorm(u)
      inner <- function(t) as.numeric(
        copulaSnpPartialMoments(t, 2L * degree)[, seq_along(c2), drop = FALSE]
          %*% c2) / norm
      ## Most points converge in a handful of Newton passes while a few tail
      ## points need two or three times as many. Iterating the whole vector
      ## until the slowest one is done wastes the majority of the work, so the
      ## loop carries an active set and shrinks it as points converge.
      active <- seq_along(z)
      for (i in seq_len(30L)) {
        value <- inner(z[active]) - u[active]
        done <- abs(value) < 1e-12
        if (all(done)) break
        keep <- active[!done]
        value <- value[!done]
        below <- value < 0
        lo[keep] <- ifelse(below, z[keep], lo[keep])
        hi[keep] <- ifelse(below, hi[keep], z[keep])
        slope <- horner(a, z[keep])^2 * stats::dnorm(z[keep]) / norm
        proposal <- z[keep] - ifelse(slope > 0, value / slope, 0)
        outside <- !is.finite(proposal) | proposal < lo[keep] |
          proposal > hi[keep]
        z[keep] <- ifelse(outside, (lo[keep] + hi[keep]) / 2, proposal)
        active <- keep
      }
      value <- location(typical, par) + par[["scale"]] * z
      if (identical(support, "positive")) exp(value) else value
    },
    random = function(n, typical, par) {
      ## rejection against the base normal: P^2 phi <= max(P^2) phi
      out <- numeric(0); a <- full(par)
      bound <- max(horner(a, seq(-8, 8, length.out = 4000L))^2)
      while (length(out) < n) {
        z <- stats::rnorm(2L * n)
        keep <- stats::runif(length(z)) < horner(a, z)^2 / bound
        out <- c(out, z[keep])
      }
      z <- out[seq_len(n)]
      value <- location(typical, par) + par[["scale"]] * z
      if (identical(support, "positive")) exp(value) else value
    })
}

## Exponentially tilted marginal density: a base family on the latent scale
## multiplied by exp of an orthonormal expansion in its own probability
## integral transform,
##
##   f(z) = exp{ sum_k beta_k L_k(2 Phi(z) - 1) } phi(z) / Z(beta),
##
## with L_k the shifted Legendre polynomials scaled to be orthonormal under the
## base. This is the same idea as the squared-polynomial SNP density but in the
## form that makes the estimation problem convex. Writing the tilt in the
## exponent puts the margin in an exponential family with sufficient statistics
## L_k, so
##
##   d log f / d beta_k  =  L_k(z) - E_beta[L_k]
##   d2 log f / d beta^2 = -Var_beta(L)
##
## AT A FIXED LOCATION the log-likelihood is therefore concave with a
## negative-definite Hessian, the maximum is unique, and combined with Fisher's
## identity the stationarity condition is exactly
##
##   (1/N) sum_i E_{pi_0i}[ L_k(z_ij) ]  =  E_beta[ L_k ]
##
## -- the fitted marginal reproduces the pooled-posterior moments of the basis.
## Anchoring ties the location to beta: the shift is scale * E_beta[Z], so that
## the geometric mean stays at the typical value, and every evaluation point
## slides when beta moves. Both the gradient and the Hessian then pick up a
## rank-one term through dE_beta[Z]/dbeta = Cov_beta(Z, L). The moment
## condition above is the stationarity condition of the fixed-location problem,
## not of the anchored one, and treating it as if it were is what stopped the
## solver converging above degree two. What survives anchoring is that the
## gradient is analytic and the objective is concave near its maximum -- the
## numerical Hessian stayed negative definite out to ||beta - beta*|| = 0.6 --
## so the fit is a small quasi-Newton ascent from beta = 0 with no restarts and
## no starting value to choose, not a general search over shapes.
##
## The basis is bounded, |L_k| <= sqrt(2k+1), which is what makes the family
## regular: Z(beta) is finite and smooth for every beta in R^degree, so the
## natural parameter space is the whole of R^degree and there is nothing to
## constrain. A Hermite expansion in the exponent would not be: it is
## integrable only when the leading term has even degree and a negative
## coefficient, and truncating it to fix that makes the tilt pile mass on the
## truncation edge. Bounded statistics also leave the tails alone -- the tilted
## density is the base density times a factor in
## [exp(-sum|beta_k| sqrt(2k+1)), exp(sum|beta_k| sqrt(2k+1))] -- so a
## lognormal base keeps lognormal tails and no tail behaviour is extrapolated
## from a region the data never visit.
##
## The expansion starts at L_1, so the tilt carries location and spread as well
## as shape and everything except the base scale is estimated inside the one
## concave problem. An expansion starting at L_3 was tried, leaving location to
## the anchor and spread to `scale`: it is better conditioned in principle but
## measurably weaker in practice, because a bounded tilt cannot reshape a base
## whose spread is wrong. At a matched number of free parameters the L_1 start
## was better on every smooth target tried. beta = 0 is exactly the base
## family, so degree 0 is lognormal on the positive line and Normal on the real
## line, and the sequence in degree is nested.
##
## Everything is evaluated on one fixed Gauss-weighted grid: the normaliser,
## the moments, the Gaussian score as the running integral put through
## qnorm, and the quantile by interpolating that score the other way. No root
## finding is needed anywhere. The score rather than the probability is the
## interpolation variable because a probability grid is useless in the tails,
## where successive values differ by less than a double can hold.

## How far the coefficients may travel. The bound is not a convenience: the
## profile over the base scale visits scales at which the unconstrained
## maximiser is far out, so the Newton step is projected into this box, and the
## same box is what the score-SA recursion is allowed to search when the margin
## is left free during a joint refit. Widening it to 50 so the fit could never
## touch it made the joint refit fail outright -- the recursion wandered and
## attained a LOWER observed log-likelihood (median -1102.7 against -1071.6
## with the margin frozen). At 12 the free refit matches the frozen likelihood
## to within importance-sampling error and improves the fitted margin
## (median KL 0.0117 against 0.0164). Tightening further to 4 costs accuracy
## again (0.0298).
copulaTiltLimit <- 12

copulaTiltBasis <- function(z, degree) {
  ## orthonormal shifted Legendre polynomials L_1 .. L_degree, evaluated in the
  ## base probability integral transform of z
  out <- matrix(0, length(z), degree)
  if (!degree) return(out)
  t <- 2 * stats::pnorm(z) - 1
  previous <- rep(1, length(z)); current <- t              # P_0, P_1
  out[, 1L] <- current * sqrt(3)
  if (degree >= 2L) for (k in 2:degree) {
    nextTerm <- ((2 * k - 1) * t * current - (k - 1) * previous) / k
    previous <- current; current <- nextTerm
    out[, k] <- current * sqrt(2 * k + 1)
  }
  out
}

## The tilt itself, beta' L(z), accumulated in place. The density path wants
## only this linear combination, and building the whole basis matrix to
## multiply it away costs an allocation of length(z) times degree and a matrix
## product on every evaluation -- and the density is the most-called thing in
## the fit.
copulaTiltShape <- function(z, beta) {
  degree <- length(beta)
  if (!degree) return(rep(0, length(z)))
  t <- 2 * stats::pnorm(z) - 1
  previous <- rep(1, length(z)); current <- t              # P_0, P_1
  shape <- current * (beta[1L] * sqrt(3))
  if (degree >= 2L) for (k in 2:degree) {
    nextTerm <- ((2 * k - 1) * t * current - (k - 1) * previous) / k
    previous <- current; current <- nextTerm
    shape <- shape + current * (beta[k] * sqrt(2 * k + 1))
  }
  shape
}

## Basis and its derivative in z together. The derivative is needed because the
## anchor shift depends on the coefficients: the standardised argument is
## z = residual / scale + E_beta[Z], so moving beta moves every evaluation point
## as well as the model moments. Both read the same normal CDF and the same
## recursion, and that CDF is the single most expensive thing in the fit, so
## they are produced in one pass.
copulaTiltBasisBoth <- function(z, degree) {
  empty <- matrix(0, length(z), degree)
  if (!degree) return(list(basis = empty, derivative = empty))
  basis <- derivative <- empty
  t <- 2 * stats::pnorm(z) - 1
  slope <- 2 * stats::dnorm(z)
  denominator <- t * t - 1
  flat <- abs(denominator) < 1e-12          # |t| = 1, where the slope vanishes
  denominator[flat] <- 1
  previous <- rep(1, length(z)); current <- t              # P_0, P_1
  for (k in seq_len(degree)) {
    if (k > 1L) {
      nextTerm <- ((2 * k - 1) * t * current - (k - 1) * previous) / k
      previous <- current; current <- nextTerm
    }
    root <- sqrt(2 * k + 1)
    basis[, k] <- current * root
    slopeTerm <- k * (t * current - previous) / denominator
    slopeTerm[flat] <- 0
    derivative[, k] <- slopeTerm * slope * root
  }
  list(basis = basis, derivative = derivative)
}

copulaTiltBasisDerivative <- function(z, degree)
  copulaTiltBasisBoth(z, degree)$derivative

copulaTiltGrid <- local({
  cache <- new.env(parent = emptyenv())
  function(degree) {
    key <- as.character(degree)
    if (!is.null(cache[[key]])) return(cache[[key]])
    limit <- 12
    z <- seq(-limit, limit, length.out = 2048L)
    step <- z[2L] - z[1L]
    weight <- rep(step, length(z))
    weight[c(1L, length(z))] <- step / 2                  # trapezoid
    value <- list(z = z, limit = limit, weight = weight * stats::dnorm(z),
      basis = copulaTiltBasis(z, degree))
    assign(key, value, envir = cache)
    value
  }
})

## The bottom rung, written out. Identical in law to the degree-0 tilt but with
## no grid behind it: the score is affine in log x, so density, score, CDF,
## quantile and sampler are all closed form, and the family keeps the name and
## the parameter layout the ladder expects.
copulaNaturalMarginTiltedBase <- function(scale, support, anchor) {
  positive <- identical(support, "positive")
  centre <- if (positive) function(typical) log(typical) else
    function(typical) typical
  copulaNaturalMargin("tilted", support, anchor, c(scale = scale), 1e-3, 10,
    roles = "scale",
    log_density = function(x, typical, par) {
      value <- stats::dnorm((if (positive) log(x) else x) - centre(typical),
        0, par[["scale"]], log = TRUE)
      if (positive) value - log(x) else value
    },
    cdf = function(x, typical, par)
      stats::pnorm((if (positive) log(x) else x) - centre(typical), 0,
        par[["scale"]]),
    quantile = function(u, typical, par) {
      value <- centre(typical) + par[["scale"]] * stats::qnorm(u)
      if (positive) exp(value) else value
    },
    random = function(n, typical, par) {
      value <- centre(typical) + par[["scale"]] * stats::rnorm(n)
      if (positive) exp(value) else value
    },
    score = function(x, typical, par)
      ((if (positive) log(x) else x) - centre(typical)) / par[["scale"]],
    inverse_score = function(z, typical, par) {
      value <- centre(typical) + par[["scale"]] * z
      if (positive) exp(value) else value
    },
    evaluate = function(x, typical, par) {
      inner <- if (positive) log(x) else x
      zeta <- (inner - centre(typical)) / par[["scale"]]
      density <- -log(par[["scale"]]) - .5 * zeta^2 - .5 * log(2 * pi)
      list(log_density = if (positive) density - inner else density,
        score = zeta)
    })
}

copulaNaturalMarginTilted <- function(degree = 2L, scale = .3, beta = NULL,
    support = c("positive", "real")) {
  support <- match.arg(support)
  degree <- as.integer(degree)
  if (degree < 0L || degree > 6L) stop("tilt degree must be between 0 and 6")
  if (is.null(beta)) beta <- rep(0, degree)
  if (length(beta) != degree) stop("beta must have length degree")
  ## Degree 0 is the base family itself, and the base family has closed forms
  ## for all of it. Going through the grid would give the same answer to
  ## quadrature error, but the bottom rung of the ladder is also the one the
  ## selection keeps whenever a coordinate shows no departure from lognormal,
  ## so it is the rung most likely to be carried into a fit -- and there it
  ## sits in the MCMC inner loop.
  if (!degree) {
    if (identical(support, "positive"))
      return(copulaNaturalMarginTiltedBase(scale, "positive",
        "geometric-mean"))
    return(copulaNaturalMarginTiltedBase(scale, "real", "mean"))
  }
  grid <- copulaTiltGrid(degree)
  names <- paste0("b", seq_len(degree))
  ## Every grid quantity -- the normaliser, the tilt centre, the CDF and its
  ## inverse -- is a function of beta alone, and the optimiser and the sampler
  ## ask for them repeatedly at the same beta. Compute them once per distinct
  ## coefficient vector. The key is the coefficient vector itself, so a miss
  ## only costs a recomputation.
  state <- local({
    lastPar <- NULL; last <- NULL; value <- NULL
    function(par) {
      ## The common case is the very same parameter vector again -- the density,
      ## the score and the anchor each ask for this within one evaluation -- so
      ## that is tested first, before anything is allocated. Positional rather
      ## than by name for the fallback: the layout is scale first and then the
      ## coefficients, and character subsetting is not free either.
      if (!is.null(lastPar) && identical(lastPar, par)) return(value)
      b <- unname(par[-1L])
      if (!is.null(last) && length(last) == length(b) && all(last == b)) {
        lastPar <<- par
        return(value)
      }
      lw <- as.numeric(grid$basis %*% b)
      peak <- max(lw); raw <- grid$weight * exp(lw - peak)
      total <- sum(raw); w <- raw / total
      ## The interpolation variable is the Gaussian score, not the probability.
      ## Both are equivalent in the bulk, but a probability grid is useless in
      ## the tails -- successive values differ by less than a double can hold,
      ## so most tail nodes collapse and the inverse loses them -- whereas the
      ## score grid stays close to the identity and interpolates evenly across
      ## the whole range. Each tail is accumulated from its own end so that
      ## neither is formed as one minus the other.
      below <- cumsum(w) - w / 2
      above <- rev(cumsum(rev(w))) - w / 2
      ## `ifelse` evaluates both branches in full, so the tail that is not
      ## selected is still passed through qnorm -- twice the work on 2048
      ## nodes. Filling each half separately is the same answer at 1.8 times
      ## the speed.
      lowerHalf <- below <= .5
      scoreGrid <- numeric(length(below))
      scoreGrid[lowerHalf] <- stats::qnorm(below[lowerHalf])
      scoreGrid[!lowerHalf] <- -stats::qnorm(above[!lowerHalf])
      keep <- c(TRUE, diff(scoreGrid) > 0) & is.finite(scoreGrid)
      last <<- b; lastPar <<- par
      value <<- list(logNormaliser = peak + log(total), tilt = lw,
        centre = sum(w * grid$z), moments = drop(w %*% grid$basis),
        scoreGrid = scoreGrid, keptScore = scoreGrid[keep],
        keptZ = grid$z[keep])
      value
    }
  })
  positive <- identical(support, "positive")
  ## The anchor is a logarithm of the typical values, and in a fit those are
  ## the same vector for every proposal of every iteration. Comparing the
  ## vector is a memcmp; taking its logarithm is not.
  anchor <- local({
    last <- NULL; value <- NULL
    function(typical) {
      if (!positive) return(typical)
      if (!is.null(last) && identical(last, typical)) return(value)
      last <<- typical; value <<- log(typical)
      value
    }
  })
  location <- function(typical, par)
    anchor(typical) - par[[1L]] * state(par)$centre
  standardise <- function(x, typical, par)
    ((if (positive) log(x) else x) - location(typical, par)) / par[[1L]]
  ## The two directions of one monotone piecewise-linear map, so composing them
  ## returns the argument exactly.
  ## The two interpolations are the hottest thing in a tilted fit, and `approx`
  ## re-sorts and de-duplicates its grid on every call. Both grids are fixed
  ## once beta is, so the regularisation is hoisted into a closure and reused.
  ##
  ## Lazily, though, and that matters. Hoisting pays only above about two calls
  ## per beta, which is true while sampling at a fixed parameter and false
  ## while fitting the margin, where the optimiser visits many trial beta and
  ## asks for one density each -- and needs neither of these maps. Building
  ## them eagerly in the memo made the construction itself ten per cent of a
  ## fit.
  interpolator <- function(makeGrid) {
    lastBeta <- NULL; f <- NULL
    function(par, current) {
      b <- unname(par[-1L])
      if (is.null(lastBeta) || !identical(lastBeta, b)) {
        f <<- makeGrid(current); lastBeta <<- b
      }
      f
    }
  }
  scoreMap <- interpolator(function(current)
    stats::approxfun(grid$z, current$scoreGrid, rule = 2))
  zMap <- interpolator(function(current)
    stats::approxfun(current$keptScore, current$keptZ, rule = 2))
  gaussianScore <- function(z, par) scoreMap(par, state(par))(z)
  inverseScore <- function(zeta, typical, par) {
    kept <- state(par)
    z <- zMap(par, kept)(zeta)
    value <- location(typical, par) + par[[1L]] * z
    if (positive) exp(value) else value
  }
  ## The density evaluates the basis exactly rather than interpolating the
  ## tabulated tilt. Interpolation would be much cheaper, but it costs 4.3e-05
  ## per row, and the selection criterion sums over every subject and every
  ## posterior draw -- several log-likelihood units on a pooled sample, against
  ## a BIC penalty of a few. The basis needs one normal CDF for the whole
  ## vector and a recursion in K after that, so the exact route is not the
  ## expensive part of the hot path anyway.
  logDensity <- function(x, typical, par, z = NULL) {
    if (is.null(z)) z <- standardise(x, typical, par)
    current <- state(par)
    value <- copulaTiltShape(z, unname(par[-1L])) +
      stats::dnorm(z, log = TRUE) - current$logNormaliser - log(par[[1L]])
    value[!is.finite(z)] <- -Inf
    if (positive) value - log(x) else value
  }
  parameters <- c(scale = scale)
  if (degree) parameters <- c(parameters,
    stats::setNames(beta, paste0("b", seq_len(degree))))
  copulaNaturalMargin("tilted", support,
    if (identical(support, "positive")) "geometric-mean" else "mean",
    parameters, c(1e-3, rep(-copulaTiltLimit, degree)),
    c(10, rep(copulaTiltLimit, degree)),
    roles = c("scale", rep("shape", degree)),
    ## Density, score, CDF and quantile all read the one cached grid. The score
    ## is the primitive and the CDF is Phi of it, so the closed-form score and
    ## the transform it replaces agree by construction rather than to within a
    ## tolerance, and the copula path never forms the uniform at all.
    log_density = logDensity,
    score = function(x, typical, par)
      gaussianScore(standardise(x, typical, par), par),
    inverse_score = function(zeta, typical, par)
      inverseScore(zeta, typical, par),
    evaluate = function(x, typical, par) {
      z <- standardise(x, typical, par)
      list(log_density = logDensity(x, typical, par, z = z),
        score = gaussianScore(z, par))
    },
    cdf = function(x, typical, par)
      stats::pnorm(gaussianScore(standardise(x, typical, par), par)),
    quantile = function(u, typical, par)
      inverseScore(stats::qnorm(u), typical, par),
    ## Sampling goes through the score directly: a standard normal draw mapped
    ## by the inverse score is the same thing as inverting the CDF at a
    ## uniform, without the round trip through the probability scale.
    random = function(n, typical, par)
      inverseScore(stats::rnorm(n), typical, par))
}

## Both of these sit in the MCMC inner loop, and both were coercing their
## argument twice -- once to count its columns and once to pass it on -- when
## every caller already holds a matrix.
copulaWorkingToNatural <- function(phi, transform) {
  if (!is.matrix(phi)) phi <- as.matrix(phi)
  transphi(phi, rep_len(as.integer(transform), ncol(phi)))
}

copulaNaturalToWorking <- function(psi, transform) {
  if (!is.matrix(psi)) psi <- as.matrix(psi)
  transpsi(psi, rep_len(as.integer(transform), ncol(psi)))
}

copulaWorkingLogJacobian <- function(phi, transform) {
  if (!is.matrix(phi)) phi <- as.matrix(phi)
  transform <- rep_len(as.integer(transform), ncol(phi))
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
    column <- copulaMarginColumn(margins[[j]], psi[, j], typical[, j])
    valid <- valid & column$valid
    z[column$valid, j] <- column$z[column$valid]
    logMargin[column$valid, j] <- column$logMargin[column$valid]
  }
  list(z = z, logMargin = logMargin, valid = valid)
}

## One margin against one column of values. Kept separate so callers that hold
## a single coordinate -- the bridge gradient, which perturbs one margin at a
## time -- can use it without wrapping vectors into matrices and subsetting
## them back out.
copulaMarginColumn <- function(margin, x, typical) {
  if (isTRUE(margin$closed_score)) {
    ## Density and Gaussian score share one affine transform, so they are
    ## produced together: the anchor is formed once and the uniform is never
    ## built, so it cannot underflow in the tails.
    both <- margin$evaluate(x, typical, margin$parameters)
    density <- both$log_density; value <- both$score
    return(list(z = value, logMargin = density,
      valid = is.finite(density) & is.finite(value)))
  }
  density <- margin$log_density(x, typical, margin$parameters)
  probability <- margin$cdf(x, typical, margin$parameters)
  valid <- is.finite(density) & is.finite(probability) &
    probability > 0 & probability < 1
  value <- rep(NA_real_, length(valid))
  value[valid] <- stats::qnorm(probability[valid])
  list(z = value, logMargin = density, valid = valid)
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

## Map Gaussian scores straight to the natural scale. Families with a closed
## form skip the uniform entirely; the rest fall back to quantile(pnorm(z)).
copulaNaturalMarginsFromScore <- function(z, typical, margins) {
  z <- as.matrix(z); typical <- as.matrix(typical)
  if (any(dim(z) != dim(typical)) || ncol(z) != length(margins))
    stop("scores, typical values, and margins do not align")
  answer <- z
  for (j in seq_along(margins)) answer[, j] <- margins[[j]]$inverse_score(
    z[, j], typical[, j], margins[[j]]$parameters)
  answer
}

## Per-column margin evaluation with the previous call remembered.
##
## The second MCMC kernel proposes one eta coordinate at a time, so between
## consecutive prior evaluations exactly one column of the batch has moved and
## the rest give the same transform, the same margin density and the same
## Jacobian as before. Comparing a column is a memcmp; recomputing one is an
## exp, a log, a density and a Gaussian score over every subject times every
## chain. Only the copula quadratic form still needs all of them.
copulaNaturalColumnEvaluator <- function(margins, typical, transform) {
  d <- length(margins)
  memo <- new.env(parent = emptyenv())
  function(phi) {
    previous <- memo$phi
    reuse <- !is.null(previous) && identical(dim(previous), dim(phi))
    if (reuse) {
      moved <- which(vapply(seq_len(d), function(j)
        !identical(previous[, j], phi[, j]), logical(1)))
      z <- memo$z; logMargin <- memo$logMargin
      jacobian <- memo$jacobian; valid <- memo$valid
    } else {
      moved <- seq_len(d)
      z <- logMargin <- jacobian <- matrix(NA_real_, nrow(phi), d)
      valid <- matrix(TRUE, nrow(phi), d)
    }
    if (length(moved)) {
      psi <- copulaWorkingToNatural(phi[, moved, drop = FALSE],
        transform[moved])
      jacobian[, moved] <- copulaWorkingLogJacobian(
        phi[, moved, drop = FALSE], transform[moved])
      for (k in seq_along(moved)) {
        j <- moved[k]
        column <- copulaMarginColumn(margins[[j]], psi[, k], typical[, j])
        z[, j] <- column$z; logMargin[, j] <- column$logMargin
        valid[, j] <- column$valid
      }
    }
    memo$phi <- phi; memo$z <- z; memo$logMargin <- logMargin
    memo$jacobian <- jacobian; memo$valid <- valid
    list(z = z, logMargin = logMargin, jacobian = jacobian,
      valid = if (d == 1L) valid[, 1L] else
        matrixStats_rowAll(valid))
  }
}

## rowSums on a logical matrix allocates a double vector; this stays integer
## and short-circuits the common single-coordinate case.
matrixStats_rowAll <- function(x) {
  answer <- x[, 1L]
  if (ncol(x) > 1L) for (j in 2:ncol(x)) answer <- answer & x[, j]
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
  decomposition <- copulaGaussianCopulaDecomposition(correlation)
  columns <- copulaNaturalColumnEvaluator(margins, typical, transform)
  negative <- function(eta) {
    eta <- if (is.matrix(eta)) eta else as.matrix(eta)
    if (any(dim(eta) != dim(predictor)))
      stop("natural-parameter eta batch does not align with its predictor")
    phi <- predictor + eta
    evaluated <- columns(phi)
    ## Every row valid is the ordinary case, and taking a row subset of three
    ## matrices to say so copies the whole batch for nothing.
    if (all(evaluated$valid))
      return(-(copulaGaussianCopulaLogDensity(evaluated$z,
        decomposition = decomposition) + rowSums(evaluated$logMargin) +
        rowSums(evaluated$jacobian)))
    answer <- rep(Inf, nrow(eta)); rows <- which(evaluated$valid)
    if (length(rows)) {
      logDensity <- copulaGaussianCopulaLogDensity(
        evaluated$z[rows, , drop = FALSE], decomposition = decomposition) +
        rowSums(evaluated$logMargin[rows, , drop = FALSE]) +
        rowSums(evaluated$jacobian[rows, , drop = FALSE])
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

## Joint natural-parameter/covariate density in saemix working coordinates.
## The first dEta margins act on natural individual parameters; remaining
## margins act directly on observed continuous covariates.
## The observed-covariate probability integral transform depends only on the
## covariate values and their margins, never on the proposed eta, yet the
## natural FREM prior is evaluated many times per iteration with both held
## fixed. Memoise the last result. The key is compared with identical(), so a
## miss is always safe: it costs a recomputation, never a stale answer.

.copulaCovariateTransformCache <- new.env(parent = emptyenv())

copulaCovariateTransform <- function(conditioning, margins) {
  ## The key is compared as numbers and identical strings rather than through
  ## format(): formatting doubles to 17 digits on every call cost more than the
  ## transform it was meant to protect.
  key <- list(names = vapply(margins, `[[`, character(1), "name"),
    types = vapply(margins, `[[`, character(1), "type"),
    parameters = unlist(lapply(margins, `[[`, "parameters")))
  entry <- .copulaCovariateTransformCache$entry
  if (!is.null(entry) && identical(entry$key, key) &&
      identical(entry$conditioning, conditioning)) return(entry$value)
  value <- copulaGaussianFremEvaluateMargins(conditioning, margins)
  .copulaCovariateTransformCache$entry <- list(conditioning = conditioning,
    key = key, value = value)
  value
}

copulaNaturalFremLogPrior <- function(eta, conditioning, vine, margins, dEta,
    predictor, transform, likelihoodTarget = c("joint", "conditional")) {
  likelihoodTarget <- match.arg(likelihoodTarget)
  eta <- as.matrix(eta); conditioning <- as.matrix(conditioning)
  predictor <- as.matrix(predictor); dEta <- as.integer(dEta)
  dConditioning <- ncol(conditioning); d <- dEta + dConditioning
  if (ncol(eta) != dEta || nrow(eta) != nrow(conditioning) ||
      any(dim(predictor) != dim(eta)) || length(margins) != d ||
      !copulaIsFullGaussianVine(vine, d))
    stop("invalid natural-parameter FREM prior dimensions")
  if (anyNA(eta) || any(!is.finite(eta)) || anyNA(conditioning) ||
      any(!is.finite(conditioning)))
    stop("natural-parameter FREM prior requires complete finite rows")
  transform <- rep_len(as.integer(transform), dEta)
  phi <- predictor + eta
  typical <- copulaWorkingToNatural(predictor, transform)
  psi <- copulaWorkingToNatural(phi, transform)
  parameter <- copulaNaturalMarginsEvaluate(psi, typical,
    margins[seq_len(dEta)])
  covariate <- if (dConditioning)
    copulaCovariateTransform(conditioning,
      margins[dEta + seq_len(dConditioning)]) else
    list(z = matrix(numeric(), nrow(eta), 0L),
      logMargin = matrix(numeric(), nrow(eta), 0L), valid = rep(TRUE, nrow(eta)))
  valid <- parameter$valid & covariate$valid
  result <- rep(-Inf, nrow(eta)); rows <- which(valid)
  if (!length(rows)) return(result)
  z <- cbind(parameter$z, covariate$z)
  R <- copulaGaussianRvineCor(vine, d)
  logJoint <- copulaGaussianCopulaLogDensity(z[rows, , drop = FALSE], R) +
    rowSums(parameter$logMargin[rows, , drop = FALSE]) +
    rowSums(covariate$logMargin[rows, , drop = FALSE]) +
    rowSums(copulaWorkingLogJacobian(phi[rows, , drop = FALSE], transform))
  if (identical(likelihoodTarget, "conditional") && dConditioning) {
    Rc <- R[dEta + seq_len(dConditioning), dEta + seq_len(dConditioning),
      drop = FALSE]
    zc <- covariate$z[rows, , drop = FALSE]
    logCov <- copulaGaussianCopulaLogDensity(zc, Rc) +
      rowSums(covariate$logMargin[rows, , drop = FALSE])
    logJoint <- logJoint - logCov
  }
  result[rows] <- logJoint
  result
}

## Exact conditional-prior kernel for natural parameters given fully observed
## continuous covariates. Conditioning is Gaussian in latent-score space;
## quantile maps and the exact working-coordinate Jacobian recover the fitted
## natural parameter law.
copulaNaturalFremConditionalKernel <- function(conditioning, vine, margins,
    dEta, predictor, transform) {
  conditioning <- as.matrix(conditioning); predictor <- as.matrix(predictor)
  dEta <- as.integer(dEta); d <- length(margins)
  if (nrow(conditioning) != nrow(predictor) || ncol(predictor) != dEta ||
      ncol(conditioning) != d - dEta || anyNA(conditioning))
    stop("invalid natural-parameter conditional-kernel dimensions")
  conditional <- copulaGaussianFremConditional(conditioning, vine, margins, dEta)
  patternState <- conditional$patternState
  for (key in names(patternState))
    patternState[[key]]$chol <- chol(patternState[[key]]$covariance)
  transform <- rep_len(as.integer(transform), dEta)
  typical <- copulaWorkingToNatural(predictor, transform)
  parameterMargins <- margins[seq_len(dEta)]
  ## Everything the conditional prior needs that does not depend on eta is
  ## formed once here rather than on every random-walk proposal: the typical
  ## values, the covariate transform, the correlation and, because the
  ## conditioning is fixed for the whole fit, the covariate copula density that
  ## the conditional target subtracts. The eta margins themselves are evaluated
  ## one column at a time and remembered, since the second MCMC kernel moves
  ## one coordinate per proposal.
  dConditioning <- d - dEta
  correlation <- copulaGaussianRvineCor(vine, d)
  covariate <- if (dConditioning)
    copulaGaussianFremEvaluateMargins(conditioning,
      margins[dEta + seq_len(dConditioning)]) else
    list(z = matrix(numeric(), nrow(predictor), 0L),
      logMargin = matrix(numeric(), nrow(predictor), 0L),
      valid = rep(TRUE, nrow(predictor)))
  if (!all(covariate$valid))
    stop("natural-parameter conditional kernel has an invalid covariate row")
  covariateLogMargin <- if (dConditioning) rowSums(covariate$logMargin) else
    rep(0, nrow(predictor))
  conditionalOffset <- if (dConditioning) {
    index <- dEta + seq_len(dConditioning)
    copulaGaussianCopulaLogDensity(covariate$z,
      correlation[index, index, drop = FALSE]) + covariateLogMargin
  } else rep(0, nrow(predictor))
  decomposition <- copulaGaussianCopulaDecomposition(correlation)
  ## The quadratic form splits into a block that moves with eta, a block that
  ## is linear in it, and a block that does not involve it at all:
  ##
  ##   z' E z  =  z_eta' E_ee z_eta + 2 z_eta' E_ec z_c + z_c' E_cc z_c
  ##
  ## and the conditioning is fixed for the whole fit, so the second block is a
  ## fixed vector against z_eta and the third is a constant. Splitting it this
  ## way keeps the per-proposal matrix product down to the eta coordinates and
  ## avoids rebuilding the full z matrix. Everything constant -- including the
  ## covariate margins and the conditional target's own offset -- is folded
  ## into one vector here.
  etaIndex <- seq_len(dEta)
  covariateIndex <- dEta + seq_len(dConditioning)
  etaBlock <- decomposition$excess[etaIndex, etaIndex, drop = FALSE]
  crossTerm <- if (dConditioning)
    2 * covariate$z %*% t(decomposition$excess[etaIndex, covariateIndex,
      drop = FALSE]) else matrix(0, nrow(predictor), dEta)
  constantTerm <- if (dConditioning)
    rowSums((covariate$z %*% decomposition$excess[covariateIndex,
      covariateIndex, drop = FALSE]) * covariate$z) else
    rep(0, nrow(predictor))
  fixedPart <- -.5 * (decomposition$logDet + constantTerm) +
    covariateLogMargin - conditionalOffset
  columns <- copulaNaturalColumnEvaluator(parameterMargins, typical, transform)
  negative <- function(eta) {
    eta <- if (is.matrix(eta)) eta else as.matrix(eta)
    if (any(dim(eta) != dim(predictor)) || anyNA(eta) || any(!is.finite(eta)))
      return(-copulaNaturalFremLogPrior(eta, conditioning, vine, margins,
        dEta, predictor, transform, "conditional"))
    phi <- predictor + eta
    evaluated <- columns(phi)
    if (all(evaluated$valid)) {
      z <- evaluated$z
      return(-(fixedPart - .5 * (rowSums((z %*% etaBlock) * z) +
        rowSums(z * crossTerm)) + rowSums(evaluated$logMargin) +
        rowSums(evaluated$jacobian)))
    }
    answer <- rep(Inf, nrow(eta)); rows <- which(evaluated$valid)
    if (length(rows)) {
      logJoint <- copulaGaussianCopulaLogDensity(
        cbind(evaluated$z, covariate$z)[rows, , drop = FALSE],
        decomposition = decomposition) +
        rowSums(evaluated$logMargin[rows, , drop = FALSE]) +
        covariateLogMargin[rows] +
        rowSums(evaluated$jacobian[rows, , drop = FALSE])
      answer[rows] <- -(logJoint - conditionalOffset[rows])
    }
    answer
  }
  random <- function() {
    z <- conditional$mean
    for (state in patternState) {
      rows <- state$rows
      z[rows, ] <- z[rows, , drop = FALSE] +
        matrix(stats::rnorm(length(rows) * dEta), ncol = dEta) %*% state$chol
    }
    psi <- copulaNaturalMarginsFromScore(z, typical,
      parameterMargins)
    copulaNaturalToWorking(psi, transform) - predictor
  }
  list(negative = negative, random = random, conditional = conditional,
    patternState = patternState,
    method = "natural-parameter-conditional-gaussian-score-kernel")
}

## Maximum-likelihood fit of the tilted family to a pooled sample, with the
## anchor applied row by row so the typical value may vary between subjects.
##
## Writing r_i for the residual of the sample against its own typical value on
## the latent scale, the standardised argument is z_i = r_i / scale + E_beta[Z]
## and the log-likelihood in beta is, up to terms free of beta,
##
##   sum_i { beta' L(z_i) + log phi(z_i) } - N log Z(beta).
##
## Were the location fixed this would be an exponential family with gradient
## mean_i L(z_i) - E_beta[L] and Hessian -N Var_beta(L). The anchor makes z_i
## itself depend on beta, so the gradient carries an extra rank-one term; it is
## still analytic, and it is the one supplied to the solver. Only `scale` is
## left to a search, and only in one dimension.

## Only two summaries of the sample enter the Newton step, and both are
## averages of a bounded smooth function of the residual. Replacing the sample
## by a fine histogram of it therefore costs nothing measurable and makes the
## iteration cost independent of how many posterior draws were pooled: without
## it every iteration touches all N rows, and N is the number of subjects times
## the number of draws each. The profile over the base scale is an average of
## the same kind, so it reads the same histogram.
##
## Two thousand bins, not eight: against eight thousand the fitted margin came
## out at the same KL to five decimal places on every target tried, for a third
## of the time. Warm-starting the sweep from the previous scale was tried at
## the same time and rejected -- it cut the evaluations by a third but landed
## the Weibull fit at KL 0.0117 against 0.00025, so the sweep starts every
## solve from beta = 0.
copulaTiltHistogram <- function(residual, bins = 2048L) {
  if (length(residual) <= bins)
    return(list(point = residual, share = rep(1 / length(residual),
      length(residual))))
  edge <- seq(min(residual), max(residual), length.out = bins + 1L)
  index <- pmin(pmax(findInterval(residual, edge, all.inside = TRUE), 1L), bins)
  count <- tabulate(index, bins)
  keep <- count > 0L
  list(point = ((edge[-1L] + edge[-(bins + 1L)]) / 2)[keep],
    share = count[keep] / length(residual))
}

copulaNaturalMarginTiltedNewton <- function(residual, scale, degree,
    tolerance = 1e-10, maxit = 60L, limit = copulaTiltLimit,
    bins = 2048L, start = NULL, histogram = NULL) {
  if (!degree) return(list(beta = numeric(0), iterations = 0L))
  grid <- copulaTiltGrid(degree)
  if (is.null(histogram)) histogram <- copulaTiltHistogram(residual, bins)
  point <- histogram$point; share <- histogram$share
  ## The objective is the anchored log-likelihood in beta, up to terms free of
  ## it. Writing c = E_beta[Z] for the anchor shift and z_b = p_b/scale + c,
  ##
  ##   l(beta) = sum_b share_b [ beta' L(z_b) + log phi(z_b) ] - log Z(beta).
  ##
  ## At a FIXED location the gradient is exactly the moment discrepancy
  ## sum_b share_b L(z_b) - E_beta[L] and the Hessian is -Var_beta(L). Under
  ## the anchor the location is tied to beta, every evaluation point slides
  ## with it, and both pick up a rank-one term through dc/dbeta =
  ## Cov_beta(Z, L). Treating the moment discrepancy as the gradient anyway is
  ## what stopped the iteration converging above degree two: it cycled at the
  ## iteration limit for every tolerance tried, at degree 4 and degree 6.
  evaluate <- function(beta) {
    tilt <- as.numeric(grid$basis %*% beta)
    peak <- max(tilt); raw <- grid$weight * exp(tilt - peak)
    total <- sum(raw); w <- raw / total
    centre <- sum(w * grid$z)
    z <- point / scale + centre
    moment <- drop(w %*% grid$basis)
    both <- copulaTiltBasisBoth(z, degree)
    ## dc/dbeta, and the chain-rule factor that the sliding points contribute:
    ## d log phi / dz = -z, and the tilt itself moves with the points.
    shift <- drop(w %*% (grid$basis * grid$z)) - centre * moment
    chain <- sum(share * (drop(both$derivative %*% beta) - z))
    list(objective = sum(share * (drop(both$basis %*% beta) +
        stats::dnorm(z, log = TRUE))) - (peak + log(total)),
      gradient = drop(share %*% both$basis) - moment + chain * shift)
  }
  ## Preconditioning the gradient by the exponential-family information gives a
  ## direction that is always uphill, but only linear convergence: it needed
  ## 130 to 200 iterations at degree 4 and above. L-BFGS-B on the same analytic
  ## gradient reaches the same maximum in a fraction of that, and the box the
  ## family declares is imposed directly rather than by projection.
  count <- 0L
  memo <- new.env(parent = emptyenv())
  at <- function(beta) {
    if (!is.null(memo$beta) && all(memo$beta == beta)) return(memo$value)
    count <<- count + 1L
    memo$beta <- beta; memo$value <- evaluate(beta)
    memo$value
  }
  fit <- try(stats::optim(if (is.null(start)) rep(0, degree) else start,
    function(beta) -at(beta)$objective,
    function(beta) -at(beta)$gradient,
    method = "L-BFGS-B", lower = rep(-limit, degree),
    upper = rep(limit, degree),
    control = list(maxit = maxit, pgtol = tolerance)), silent = TRUE)
  if (inherits(fit, "try-error"))
    return(list(beta = if (is.null(start)) rep(0, degree) else start,
      iterations = count))
  list(beta = fit$par, iterations = count)
}

copulaNaturalMarginTiltedFit <- function(x, typical, degree = 2L,
    support = c("positive", "real"), scalePoints = 9L) {
  support <- match.arg(support)
  degree <- as.integer(degree)
  inner <- if (identical(support, "positive")) log(x) else x
  centre <- if (identical(support, "positive")) log(typical) else typical
  residual <- inner - centre
  finite <- is.finite(residual)
  if (!any(finite)) stop("tilted margin fit needs finite observations")
  residual <- residual[finite]
  iterations <- 0L
  histogram <- copulaTiltHistogram(residual)
  grid <- copulaTiltGrid(degree)
  ## At degree 0 the base scale has a closed-form maximum likelihood estimate
  ## and there is nothing else to fit, so the sweep is skipped entirely.
  scale <- sqrt(mean(residual^2))
  if (degree) {
    ## The profile in the base scale is smooth but not concave -- a coarse
    ## sweep followed by a local refinement is what keeps a bimodal target from
    ## settling on the wrong side of it, and each evaluation is only a Newton
    ## solve on the concave inner problem. The profile is evaluated on the same
    ## histogram: against the full sample it differs by the term -sum(log x),
    ## which does not involve any parameter.
    floor <- -.Machine$double.xmax / 100
    warm <- NULL
    evaluate <- function(v) {
      trial <- exp(v)
      fit <- copulaNaturalMarginTiltedNewton(residual, trial, degree,
        start = warm, histogram = histogram)
      ## The sweep visits neighbouring scales, so the previous solution is a
      ## good starting point; the inner problem is concave, so a warm start
      ## changes how fast it converges and not where it converges to.
      warm <<- fit$beta
      iterations <<- iterations + fit$iterations
      tilt <- as.numeric(grid$basis %*% fit$beta)
      peak <- max(tilt); w <- grid$weight * exp(tilt - peak)
      logNormaliser <- peak + log(sum(w)); w <- w / sum(w)
      z <- histogram$point / trial + sum(w * grid$z)
      value <- sum(histogram$share * (stats::approx(grid$z, tilt, z,
        rule = 2)$y + stats::dnorm(z, log = TRUE))) - logNormaliser - v
      if (is.finite(value)) value else floor
    }
    ## Nine points, not fifteen: the sweep only has to bracket the maximum for
    ## the refinement below, and a coarser bracket leaves that refinement a
    ## wider window to search. Against fifteen points it was never worse on the
    ## targets tried and sometimes better, at two thirds of the cost.
    sweep <- seq(log(scale) - 1.5, log(scale) + 1.5,
      length.out = max(3L, as.integer(scalePoints)))
    values <- vapply(sweep, evaluate, numeric(1))
    best <- which.max(values)
    window <- sweep[c(max(1L, best - 1L), min(length(sweep), best + 1L))]
    warm <- NULL
    scale <- exp(if (window[1L] < window[2L])
      stats::optimize(evaluate, window, maximum = TRUE, tol = 1e-4)$maximum else
      sweep[best])
  }
  fit <- copulaNaturalMarginTiltedNewton(residual, scale, degree,
    histogram = histogram)
  margin <- copulaNaturalMarginTilted(degree, scale,
    if (degree) fit$beta else NULL, support)
  list(margin = margin, scale = scale, beta = fit$beta,
    newtonIterations = iterations + fit$iterations,
    logLik = sum(margin$log_density(x, typical, margin$parameters)))
}
