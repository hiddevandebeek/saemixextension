## Exact parameterization of one Gaussian correlation matrix on any full
## R-vine structure.  The structure changes coordinates, not the model class.

copulaRVineEdgeSets <- function(structure) {
  if (!inherits(structure, "rvine_structure"))
    stop("structure must inherit from rvine_structure")
  order <- as.integer(structure$order)
  array <- structure$struct_array
  lapply(seq_len(structure$trunc_lvl), function(tree) {
    lapply(seq_len(structure$d - tree), function(edge) {
      conditioned <- c(order[edge], order[array[[tree]][edge]])
      conditioning <- if (tree == 1L) integer() else
        order[as.integer(vapply(seq_len(tree - 1L), function(t)
          array[[t]][edge], numeric(1)))]
      list(conditioned = conditioned, conditioning = conditioning)
    })
  })
}

copulaValidateCorrelation <- function(correlation) {
  correlation <- as.matrix(correlation)
  d <- nrow(correlation)
  if (ncol(correlation) != d || d < 2L || any(!is.finite(correlation)) ||
      max(abs(correlation - t(correlation))) > 1e-10 ||
      any(abs(diag(correlation) - 1) > 1e-10) ||
      min(eigen(correlation, symmetric = TRUE,
        only.values = TRUE)$values) <= 0)
    stop("correlation must be a positive-definite correlation matrix")
  correlation
}

copulaGaussianRvineFromCor <- function(correlation, structure = NULL) {
  correlation <- copulaValidateCorrelation(correlation)
  d <- nrow(correlation)
  if (is.null(structure))
    structure <- rvinecopulib::dvine_structure(seq_len(d))
  if (!inherits(structure, "rvine_structure") || structure$d != d)
    stop("correlation and structure dimensions differ")
  if (as.integer(structure$trunc_lvl) != d - 1L)
    stop("exact Gaussian correlation mapping requires a full R-vine")
  edges <- copulaRVineEdgeSets(structure)
  pairs <- lapply(edges, function(tree) lapply(tree, function(edge) {
    a <- edge$conditioned[1L]
    b <- edge$conditioned[2L]
    given <- edge$conditioning
    rho <- if (!length(given)) correlation[a, b] else {
      precision <- solve(correlation[c(a, b, given),
        c(a, b, given), drop = FALSE])
      -precision[1L, 2L] /
        sqrt(precision[1L, 1L] * precision[2L, 2L])
    }
    rvinecopulib::bicop_dist("gaussian",
      parameters = max(-1 + 1e-10, min(1 - 1e-10, as.numeric(rho))))
  }))
  rvinecopulib::vinecop_dist(pairs, structure)
}

copulaGaussianRvineUpdateCor <- function(vine, correlation) {
  correlation <- copulaValidateCorrelation(correlation)
  d <- nrow(correlation)
  if (!copulaIsFullGaussianVine(vine, d))
    stop("correlation update requires a full Gaussian R-vine")
  updated <- vine
  edges <- copulaRVineEdgeSets(vine$structure)
  for (tree in seq_along(edges)) for (edge in seq_along(edges[[tree]])) {
    info <- edges[[tree]][[edge]]
    a <- info$conditioned[1L]
    b <- info$conditioned[2L]
    given <- info$conditioning
    rho <- if (!length(given)) correlation[a, b] else {
      precision <- solve(correlation[c(a, b, given),
        c(a, b, given), drop = FALSE])
      -precision[1L, 2L] /
        sqrt(precision[1L, 1L] * precision[2L, 2L])
    }
    pair <- updated$pair_copulas[[tree]][[edge]]
    pair$family <- "gaussian"
    pair$rotation <- 0
    pair$parameters <- matrix(
      max(-1 + 1e-10, min(1 - 1e-10, as.numeric(rho))), 1L, 1L)
    pair$npars <- 1L
    updated$pair_copulas[[tree]][[edge]] <- pair
  }
  updated$npars <- as.integer(d * (d - 1L) / 2L)
  updated$loglik <- NA_real_
  if (isTRUE(getOption("saemix.copula.auditFastGaussianUpdate", FALSE)) &&
      max(abs(copulaGaussianRvineCor(updated, d) - correlation)) > 2e-10)
    stop("fast Gaussian R-vine correlation update failed reconstruction")
  updated
}

copulaGaussianRvineCor <- function(vine,
                                   d = as.integer(vine$structure$d)) {
  if (!inherits(vine, "vinecop_dist") || vine$structure$d != d)
    stop("vine and requested dimensions differ")
  if (as.integer(vine$structure$trunc_lvl) != d - 1L)
    stop("an unrestricted Gaussian correlation requires a full R-vine")
  edges <- copulaRVineEdgeSets(vine$structure)
  correlation <- diag(d)
  for (tree in seq_along(edges)) for (edge in seq_along(edges[[tree]])) {
    pair <- vine$pair_copulas[[tree]][[edge]]
    if (!pair$family %in% c("gaussian", "indep"))
      stop("Gaussian mapping requires Gaussian or independence edges")
    partial <- if (identical(pair$family, "indep")) 0 else
      as.numeric(pair$parameters)[1L]
    info <- edges[[tree]][[edge]]
    a <- info$conditioned[1L]
    b <- info$conditioned[2L]
    given <- info$conditioning
    value <- if (!length(given)) partial else {
      inverse <- solve(correlation[given, given, drop = FALSE])
      left <- correlation[a, given, drop = FALSE]
      right <- correlation[b, given, drop = FALSE]
      projected <- as.numeric(left %*% inverse %*% t(right))
      varianceA <- as.numeric(1 - left %*% inverse %*% t(left))
      varianceB <- as.numeric(1 - right %*% inverse %*% t(right))
      if (!is.finite(varianceA) || !is.finite(varianceB) ||
          varianceA <= 0 || varianceB <= 0)
        stop("invalid partial-correlation reconstruction")
      projected + partial * sqrt(varianceA * varianceB)
    }
    correlation[a, b] <- correlation[b, a] <- value
  }
  copulaValidateCorrelation(correlation)
}

copulaGaussianDvineFromCor <- function(correlation) {
  correlation <- copulaValidateCorrelation(correlation)
  copulaGaussianRvineFromCor(correlation,
    rvinecopulib::dvine_structure(seq_len(nrow(correlation))))
}

copulaGaussianFremEnsurePd <- function(S, relativeFloor = 1e-10) {
  S <- (S + t(S)) / 2
  ee <- eigen(S, symmetric = TRUE)
  floor <- max(max(ee$values), 1) * relativeFloor
  if (min(ee$values) < floor)
    S <- ee$vectors %*% (pmax(ee$values, floor) * t(ee$vectors))
  (S + t(S)) / 2
}

## A Gaussian observed-data log density is quadratic within each missingness
## pattern.  With a linear location design it is quadratic jointly in the
## observed transformed coordinates and the design row.  Total weight, first
## moments and cross-moments of that augmented vector are therefore a complete
## finite statistic; design rows need not be identical.  The symmetric-point
## representation below preserves those moments exactly and prevents the SA
## empirical-particle history from making the MVN M-step grow with k.

copulaGaussianCorrelationAngles <- function(R) {
  R <- as.matrix(R); d <- nrow(R)
  if (ncol(R) != d || d < 1L) stop("R must be a square correlation matrix")
  L <- t(chol(R)); answer <- numeric(d * (d - 1L) / 2L); cursor <- 0L
  if (d > 1L) for (i in 2:d) {
    product <- 1
    for (j in seq_len(i - 1L)) {
      ratio <- L[i, j] / product
      ratio <- max(-1, min(1, ratio))
      angle <- acos(ratio)
      cursor <- cursor + 1L; answer[cursor] <- angle
      product <- product * sin(angle)
    }
  }
  answer
}

copulaGaussianCorrelationFromAngles <- function(angles, d) {
  d <- as.integer(d); angles <- as.numeric(angles)
  if (length(angles) != d * (d - 1L) / 2L)
    stop("wrong number of hyperspherical correlation angles")
  L <- diag(d); cursor <- 0L
  if (d > 1L) for (i in 2:d) {
    product <- 1
    for (j in seq_len(i - 1L)) {
      cursor <- cursor + 1L; angle <- angles[cursor]
      L[i, j] <- product * cos(angle)
      product <- product * sin(angle)
    }
    L[i, i] <- product
  }
  list(R = tcrossprod(L), L = L)
}

copulaGaussianCorrelationAngleObjective <- function(angles, S, d = nrow(S)) {
  S <- as.matrix(S); built <- copulaGaussianCorrelationFromAngles(angles, d)
  U <- chol(built$R); Omega <- chol2inv(U)
  value <- .5 * (2 * sum(log(diag(U))) + sum(Omega * S))
  G <- .5 * (Omega - Omega %*% S %*% Omega)
  dLd <- 2 * G %*% built$L
  gradient <- numeric(length(angles)); cursor <- 0L
  if (d > 1L) for (i in 2:d) {
    rowAngles <- angles[cursor + seq_len(i - 1L)]
    for (k in seq_len(i - 1L)) {
      derivative <- numeric(i)
      prefix <- if (k == 1L) 1 else prod(sin(rowAngles[seq_len(k - 1L)]))
      derivative[k] <- -prefix * sin(rowAngles[k])
      if (k < i - 1L) derivative[(k + 1L):(i - 1L)] <-
        built$L[i, (k + 1L):(i - 1L)] / tan(rowAngles[k])
      derivative[i] <- built$L[i, i] / tan(rowAngles[k])
      gradient[cursor + k] <- sum(dLd[i, seq_len(i)] * derivative)
    }
    cursor <- cursor + i - 1L
  }
  list(value = value, gradient = gradient,
    R = built$R, L = built$L, Omega = Omega)
}
