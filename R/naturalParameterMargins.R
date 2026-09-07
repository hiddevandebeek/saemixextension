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
