## Standard errors and a covariance matrix on the scale a reader reports.
##
## The recursion accumulates its Fisher information in the coordinates it steps
## in, which are not the ones anybody wants to see: population locations sit on
## the working scale, so a location standard error is for log V rather than V,
## and the dependence is carried as Cholesky angles, so a "correlation"
## standard error is for an angle in radians. Handing that matrix out invites
## exactly the mistake of reading a working-scale error next to a natural-scale
## estimate, which is a factor of the typical value.
##
## So the map to reportable coordinates is applied here, once, with its
## Jacobian: typical values on their own scale, marginal parameters under their
## own names, pairwise correlations rather than angles, and the residual error.
## The delta method is exact for a linear map and the standard first-order
## approximation otherwise, and the Jacobian is taken numerically because the
## parameter vector is short enough that the cost is nothing.

## Gauss-Hermite nodes and weights for the standard normal, by Golub and
## Welsch: the probabilists' Hermite recursion is the symmetric tridiagonal
## matrix with zero diagonal and sqrt(k) off it, whose eigenvalues are the
## nodes and whose first eigenvector components square to the weights.
copulaNormalQuadrature <- local({
  cache <- list()
  function(nodes = 160L) {
    key <- as.character(nodes)
    if (!is.null(cache[[key]])) return(cache[[key]])
    off <- sqrt(seq_len(nodes - 1L))
    jacobi <- matrix(0, nodes, nodes)
    jacobi[cbind(seq_len(nodes - 1L), 2:nodes)] <- off
    jacobi[cbind(2:nodes, seq_len(nodes - 1L))] <- off
    decomposed <- eigen(jacobi, symmetric = TRUE)
    order <- order(decomposed$values)
    answer <- list(node = decomposed$values[order],
      weight = (decomposed$vectors[1L, order])^2)
    cache[[key]] <<- answer
    answer
  }
})

## The quantities a reader actually reports for a random effect, for any margin
## in this family.
##
## The classical parameterisation writes psi = theta * exp(eta) with eta normal,
## and there theta is at once the median and the geometric mean, and omega is at
## once the spread of log psi and, near enough, the coefficient of variation.
## A tilted margin separates all of those. What survives exactly is the
## geometric mean: the location is anchored so that E[log psi] = log theta for
## every tilt coefficient, which is what makes the degrees of the ladder
## comparable. What does not survive is the reading of the scale parameter as
## an inter-individual spread -- against a true spread it can be a third out in
## either direction once a coordinate takes a tilt.
##
## So the derived quantities are computed rather than assumed. Every margin
## here is built so that psi = inverse_score(zeta) with zeta standard normal,
## which makes each of these an ordinary Gaussian integral and lets one
## quadrature serve every family.
copulaMarginDerived <- function(margin, typical, parameters, nodes = 160L) {
  quadrature <- copulaNormalQuadrature(nodes)
  value <- try(margin$inverse_score(quadrature$node, typical, parameters),
    silent = TRUE)
  if (inherits(value, "try-error") || any(!is.finite(value)))
    return(c(median = NA_real_, spread = NA_real_, cv = NA_real_))
  w <- quadrature$weight
  positive <- identical(margin$support, "positive")
  ## zeta = 0 is the median by construction, since the margin's Gaussian score
  ## is defined so that its CDF is Phi of it.
  median <- try(margin$inverse_score(0, typical, parameters), silent = TRUE)
  if (inherits(median, "try-error")) median <- NA_real_
  inner <- if (positive) log(value) else value
  centre <- sum(w * inner)
  spread <- sqrt(max(0, sum(w * (inner - centre)^2)))
  mean <- sum(w * value)
  variance <- max(0, sum(w * (value - mean)^2))
  c(median = as.numeric(median), spread = spread,
    cv = if (is.finite(mean) && abs(mean) > 0) sqrt(variance) / abs(mean) else
      NA_real_)
}

copulaReportedParameters <- function(object) {
  state <- copulaGet(object)
  covariance <- state$fisherCovariance
  if (is.null(covariance) || is.null(rownames(covariance)))
    stop("this fit carries no Fisher information; it is accumulated for ",
      "every score-sa fit unless options(saemix.fisherInformation = FALSE) ",
      "was set", call. = FALSE)
  labels <- rownames(covariance)
  index <- as.integer(state$etaIndex)
  transform <- as.integer(object["model"]["transform.par"][index])
  parameterNames <- as.character(object["model"]["name.modpar"])[index]
  if (!length(parameterNames) || anyNA(parameterNames))
    parameterNames <- paste0("parameter", seq_along(index))
  margins <- state$margins
  free <- lapply(margins, function(margin) which(margin$free))
  angles <- copulaGaussianCorrelationAngles(
    copulaGaussianRvineCor(state$vine, state$d))
  residual <- as.numeric(object["results"]["respar"])

  ## The value of every coordinate the covariance is indexed by, in its order.
  native <- vapply(labels, function(label) {
    part <- strsplit(label, ".", fixed = TRUE)[[1L]]
    switch(part[1L],
      location = object["results"]["mean.phi"][1L, index[as.integer(part[2L])]],
      margin = unname(margins[[as.integer(part[2L])]]$parameters[
        free[[as.integer(part[2L])]][as.integer(part[3L])]]),
      correlation = angles[as.integer(part[2L])],
      residual = residual[as.integer(part[2L])],
      NA_real_)
  }, numeric(1))

  ## The same coordinates, mapped to what gets reported.
  report <- function(value) {
    typical <- numeric(0); typicalNames <- character(0)
    marginValue <- numeric(0); marginNames <- character(0)
    correlation <- numeric(0); correlationNames <- character(0)
    residualValue <- numeric(0); residualNames <- character(0)
    currentAngles <- angles
    ## A coordinate's typical value and its margin coefficients arrive in
    ## different blocks of the vector, and a derived quantity needs both, so
    ## they are gathered first. A coordinate whose location is not estimated
    ## keeps its fitted value, which is constant and so contributes nothing to
    ## the Jacobian, which is the right answer rather than a missing one.
    currentTypical <- vapply(seq_along(margins), function(j)
      if (j <= length(index))
        object["results"]["mean.phi"][1L, index[j]] else NA_real_, numeric(1))
    currentTypical <- copulaWorkingToNatural(
      matrix(currentTypical, 1L), transform)[1L, ]
    currentMargin <- lapply(free, function(f) numeric(length(f)))
    for (k in seq_along(labels)) {
      part <- strsplit(labels[k], ".", fixed = TRUE)[[1L]]
      if (identical(part[1L], "location")) {
        j <- as.integer(part[2L])
        currentTypical[j] <- copulaWorkingToNatural(
          matrix(value[k], 1L, 1L), transform[j])[1L, 1L]
      } else if (identical(part[1L], "margin"))
        currentMargin[[as.integer(part[2L])]][as.integer(part[3L])] <- value[k]
    }
    for (k in seq_along(labels)) {
      part <- strsplit(labels[k], ".", fixed = TRUE)[[1L]]
      if (identical(part[1L], "correlation"))
        currentAngles[as.integer(part[2L])] <- value[k]
    }
    for (k in seq_along(labels)) {
      part <- strsplit(labels[k], ".", fixed = TRUE)[[1L]]
      if (identical(part[1L], "location")) {
        j <- as.integer(part[2L])
        typical <- c(typical, copulaWorkingToNatural(
          matrix(value[k], 1L, 1L), transform[j])[1L, 1L])
        typicalNames <- c(typicalNames, parameterNames[j])
      } else if (identical(part[1L], "margin")) {
        j <- as.integer(part[2L]); l <- as.integer(part[3L])
        marginValue <- c(marginValue, value[k])
        ## Named by the coordinate the margin belongs to, not by the family
        ## alone: two coordinates can carry the same family, and two rows both
        ## called "tilted.scale" say nothing about which parameter is which.
        owner <- if (j <= length(parameterNames)) parameterNames[j] else
          paste0("covariate", j - length(parameterNames))
        marginNames <- c(marginNames, paste(owner, margins[[j]]$name,
          names(margins[[j]]$parameters)[free[[j]][l]], sep = "."))
      } else if (identical(part[1L], "residual")) {
        residualValue <- c(residualValue, value[k])
        residualNames <- c(residualNames, paste0("residual.", part[2L]))
      }
    }
    if (length(currentAngles)) {
      R <- copulaGaussianCorrelationFromAngles(currentAngles, state$d)$R
      pairs <- which(lower.tri(R), arr.ind = TRUE)
      correlation <- R[lower.tri(R)]
      correlationNames <- paste0("correlation.", pairs[, 2L], ".",
        pairs[, 1L])
    }
    ## The derived quantities. They are functions of the typical value and the
    ## margin coefficients together, so they belong inside `report` rather than
    ## beside it: the numerical Jacobian below then carries their dependence on
    ## every coordinate at once, and the delta method needs no special case.
    derivedValue <- numeric(0); derivedNames <- character(0)
    for (j in seq_along(margins)) {
      owner <- if (j <= length(parameterNames)) parameterNames[j] else
        paste0("covariate", j - length(parameterNames))
      here <- currentTypical[j]
      par <- margins[[j]]$parameters
      if (length(free[[j]])) par[free[[j]]] <- currentMargin[[j]]
      derived <- copulaMarginDerived(margins[[j]], here, par)
      positive <- identical(margins[[j]]$support, "positive")
      derivedValue <- c(derivedValue, derived[["median"]],
        derived[["spread"]], derived[["cv"]])
      derivedNames <- c(derivedNames, paste0(owner, ".median"),
        paste0(owner, if (positive) ".sd.log" else ".sd"),
        paste0(owner, ".cv"))
    }
    stats::setNames(c(typical, marginValue, derivedValue, correlation,
      residualValue),
      c(typicalNames, marginNames, derivedNames, correlationNames,
        residualNames))
  }

  value <- report(native)
  jacobian <- matrix(0, length(value), length(native))
  for (k in seq_along(native)) {
    step <- 1e-6 * max(1, abs(native[k]))
    up <- down <- native; up[k] <- up[k] + step; down[k] <- down[k] - step
    jacobian[, k] <- (report(up) - report(down)) / (2 * step)
  }
  reported <- jacobian %*% covariance %*% t(jacobian)
  dimnames(reported) <- list(names(value), names(value))
  list(estimate = value, covariance = reported)
}

copulaStandardErrors <- function(object, level = .95) {
  got <- copulaReportedParameters(object)
  error <- sqrt(pmax(diag(got$covariance), 0))
  quantile <- stats::qnorm(1 - (1 - level) / 2)
  data.frame(parameter = names(got$estimate),
    estimate = as.numeric(got$estimate), se = as.numeric(error),
    lower = as.numeric(got$estimate) - quantile * error,
    upper = as.numeric(got$estimate) + quantile * error,
    row.names = NULL, stringsAsFactors = FALSE)
}
