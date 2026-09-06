## Misspecification-robust intervals.
##
## The score-covariance estimate is the outer-product form. It equals the
## Fisher information when the model is right, and then its inverse is the
## covariance of the estimate; when the model is wrong the two come apart and
## the honest interval is the sandwich
##
##   H^-1 J H^-1 / N,
##
## with J the covariance of the per-subject scores -- already accumulated -- and
## H the observed information per subject. That distinction is not academic
## here: this workflow exists to relax a distributional assumption, so the case
## where the assumption is wrong is the case it is for. Measured on a fit whose
## clearance was Gamma while the model said Normal, the outer-product interval
## for one coordinate was 0.17 of the spread it should have covered; with a
## correctly specified truth the same coordinate gave 0.39, and with the model
## right and converged, the ratios sat between 0.81 and 1.25.
##
## The bread comes from Louis (1982), which writes the observed information of
## a latent-variable model using the complete-data density alone:
##
##   -grad^2 log g_i = E[-grad^2 log f_i | y_i] - Var[grad log f_i | y_i].
##
## Both terms are posterior expectations over the same draws the score already
## needs. The second is the conditional variance of the complete-data score,
## which is a by-product of computing the score itself; the first needs second
## derivatives of the complete-data log density, and those are cheap because
## that density is a closed-form expression over subjects rather than anything
## requiring the MCMC.

copulaScoreSandwich <- function(object, draws = 200L, max.iter = NULL,
    seed = 935001L, step = 1e-4, hessianStep = 1e-3, thin = 1L) {
  context <- copulaFisherContext(object)
  perSubject <- context$perSubject
  theta <- context$theta
  n <- context$subjects
  p <- length(theta)

  pool <- copulaPosteriorEtaDraws(object, draws, max.iter, as.integer(seed),
    thin = thin)
  samples <- dim(pool)[3L]
  if (samples < 4L)
    stop("the sandwich needs at least four posterior draws", call. = FALSE)

  gradientAt <- function(value, eta) {
    answer <- matrix(0, n, p)
    for (k in seq_len(p)) {
      h <- step * max(1, abs(value[k]))
      up <- down <- value; up[k] <- up[k] + h; down[k] <- down[k] - h
      answer[, k] <- (perSubject(up, eta) - perSubject(down, eta)) / (2 * h)
    }
    answer
  }

  half <- rep(c(1L, 2L), length.out = samples)
  score <- matrix(0, n, p)
  halves <- list(matrix(0, n, p), matrix(0, n, p))
  ## The conditional variance of the complete-data score, and the conditional
  ## expectation of its second derivative: the two halves of Louis' identity.
  outer <- array(0, c(n, p, p))
  curvature <- array(0, c(n, p, p))
  for (r in seq_len(samples)) {
    eta <- pool[, , r, drop = FALSE]
    dim(eta) <- dim(pool)[1:2]
    gradient <- gradientAt(theta, eta)
    score <- score + gradient / samples
    halves[[half[r]]] <- halves[[half[r]]] +
      gradient / sum(half == half[r])
    for (a in seq_len(p)) for (b in seq_len(a))
      outer[, a, b] <- outer[, a, b] + gradient[, a] * gradient[, b] / samples

    base <- perSubject(theta, eta)
    for (a in seq_len(p)) {
      ha <- hessianStep * max(1, abs(theta[a]))
      for (b in seq_len(a)) {
        hb <- hessianStep * max(1, abs(theta[b]))
        shift <- function(sa, sb) {
          value <- theta
          value[a] <- value[a] + sa * ha
          value[b] <- value[b] + sb * hb
          perSubject(value, eta)
        }
        second <- if (a == b)
          (shift(1, 0) - 2 * base + shift(-1, 0)) / (ha * ha) else
          (shift(1, 1) - shift(1, -1) - shift(-1, 1) + shift(-1, -1)) /
            (4 * ha * hb)
        curvature[, a, b] <- curvature[, a, b] + second / samples
      }
    }
  }
  symmetrise <- function(x) {
    for (a in seq_len(p)) for (b in seq_len(a)) x[, b, a] <- x[, a, b]
    x
  }
  outer <- symmetrise(outer); curvature <- symmetrise(curvature)

  ## Louis, per subject, then averaged.
  bread <- matrix(0, p, p)
  for (a in seq_len(p)) for (b in seq_len(p))
    bread[a, b] <- mean(-curvature[, a, b] -
      (outer[, a, b] - score[, a] * score[, b]))
  bread <- (bread + t(bread)) / 2
  cross <- crossprod(halves[[1L]], halves[[2L]]) / n
  meat <- (cross + t(cross)) / 2

  inverse <- try(solve(bread), silent = TRUE)
  covariance <- if (inherits(inverse, "try-error")) NULL else {
    answer <- inverse %*% meat %*% inverse / n
    dimnames(answer) <- list(names(theta), names(theta))
    (answer + t(answer)) / 2
  }
  dimnames(bread) <- dimnames(meat) <- list(names(theta), names(theta))
  list(bread = bread, meat = meat, covariance = covariance,
    standardErrors = if (is.null(covariance)) NULL else
      stats::setNames(sqrt(pmax(diag(covariance), 0)), names(theta)),
    modelBased = local({
      solved <- try(solve(meat), silent = TRUE)
      if (inherits(solved, "try-error")) NULL else
        stats::setNames(sqrt(pmax(diag(solved) / n, 0)), names(theta))
    }),
    parameters = theta, draws = samples,
    method = "louis-bread, score-covariance meat")
}
