## Did the recursion converge?
##
## The fit records what is needed to answer this on every iteration -- the
## score, the parameter vector, the gain, and how often the truncation and the
## step-halving had to intervene -- but nothing reads it back, so the question
## has been answered by looking at the estimates and hoping. It is worth asking
## directly, because a metric-dependent difference in the estimates turned out
## to mean exactly this and nothing else: two runs that must share a fixed
## point had not reached it.
##
## Four things are checked, each of which can fail on its own.
##
##   The score has to be zero. It is a Monte Carlo estimate, so it never is,
##   and comparing it against its own past cannot separate "still moving" from
##   "arrived at the noise floor" -- that comparison called converged fits
##   unsettled here, on fits whose estimates were identical at 800 and 2500
##   iterations. The per-subject scores answer it properly, and are tested
##   below; the window ratio is kept only for fits that lack them.
##
##   The steps have to shrink. Stochastic approximation converges because the
##   gain decays, so the distance moved per iteration should decay with it; a
##   displacement that is flat means the recursion is still travelling.
##
##   The truncation should be idle. Andrieu et al.'s framework tolerates
##   excursions from the compact set, but the convergence result describes a
##   tail that stays inside one. Reprojections after the adaptation freeze mean
##   the realised trajectory is not that tail.
##
##   The metric has to have settled, since the theory is written for a
##   preconditioner that is eventually fixed.

copulaConvergence <- function(object, tail = .2) {
  state <- copulaGet(object)
  trace <- state$trace
  if (!length(trace))
    stop("this fit recorded no score-sa trace", call. = FALSE)
  iterations <- length(trace)
  window <- max(2L, as.integer(iterations * tail))
  if (iterations < 4L * window)
    warning("too few iterations to judge convergence from the tail alone",
      call. = FALSE)

  parameters <- t(vapply(trace, function(step) c(step$beta, step$margin,
    step$copula, step$residual), numeric(length(c(trace[[1L]]$beta,
      trace[[1L]]$margin, trace[[1L]]$copula, trace[[1L]]$residual)))))
  score <- vapply(trace, function(step) step$scoreMax %||% NA_real_,
    numeric(1))
  gain <- vapply(trace, function(step) step$gamma %||% NA_real_, numeric(1))
  displacement <- c(NA_real_, sqrt(rowSums(diff(parameters)^2)))

  last <- seq.int(iterations - window + 1L, iterations)
  previous <- seq.int(iterations - 2L * window + 1L, iterations - window)
  scoreRatio <- mean(score[last], na.rm = TRUE) /
    mean(score[previous], na.rm = TRUE)

  ## The stationarity test, when the fit carried the per-subject scores. They
  ## are independent across subjects, so their average has covariance V/N and
  ## under the hypothesis that the observed score is zero
  ##   T = N gbar' I^-1 gbar ~ chi-square on p degrees of freedom.
  ## That asks the right question, where comparing a Monte Carlo score against
  ## its own past cannot separate "still moving" from "arrived at the floor".
  delta <- state$scoreState$deltaSubject
  information <- state$scoreState$fisher
  statistic <- pValue <- NA_real_
  if (!is.null(delta) && !is.null(information)) {
    average <- colMeans(delta)
    solved <- try(solve(information, average), silent = TRUE)
    if (!inherits(solved, "try-error")) {
      statistic <- nrow(delta) * sum(average * solved)
      if (is.finite(statistic) && statistic >= 0)
        pValue <- stats::pchisq(statistic, df = length(average),
          lower.tail = FALSE)
    }
  }
  ## The step should shrink at least as fast as the gain that drives it.
  displacementRatio <- (mean(displacement[last], na.rm = TRUE) /
    mean(displacement[previous], na.rm = TRUE)) /
    (mean(gain[last], na.rm = TRUE) / mean(gain[previous], na.rm = TRUE))

  joint <- state$lastJoint
  projections <- joint$postFreezeProjectionCount %||% 0L
  backtracks <- joint$postFreezeBacktrackCount %||% 0L
  noMove <- joint$postFreezeNoMoveCount %||% 0L
  conditioning <- if (is.null(joint$preconditionerMin) ||
      is.null(joint$preconditionerMax)) NA_real_ else
    joint$preconditionerMax / joint$preconditionerMin

  ## The test decides when it is available; the window ratio is a fallback for
  ## fits that did not carry the per-subject scores.
  stationary <- if (is.finite(pValue)) pValue > .01 else
    isTRUE(scoreRatio <= 1.1)
  settled <- stationary && isTRUE(displacementRatio <= 1.5) &&
    projections == 0L && noMove == 0L
  answer <- list(iterations = iterations, window = window,
    statistic = statistic, pValue = pValue,
    scoreRatio = scoreRatio, displacementRatio = displacementRatio,
    projections = projections, backtracks = backtracks, noMove = noMove,
    preconditionerConditioning = conditioning,
    finalScore = score[iterations], settled = settled)
  class(answer) <- c("saemixScoreConvergence", "list")
  answer
}

print.saemixScoreConvergence <- function(x, ...) {
  cat("Score-SA convergence\n")
  cat(sprintf("  iterations                         %d (tail window %d)\n",
    x$iterations, x$window))
  if (is.finite(x$pValue))
    cat(sprintf("  score is zero: chi-square %.1f, p = %.3f %s\n",
      x$statistic, x$pValue,
      if (x$pValue > .01) "" else "  <- score still detectable"))
  else
    cat(sprintf("  score, last window over previous   %.3f %s\n", x$scoreRatio,
      if (isTRUE(x$scoreRatio <= 1.1)) "" else "  <- not falling"))
  cat(sprintf("  step decay relative to the gain    %.3f %s\n",
    x$displacementRatio,
    if (isTRUE(x$displacementRatio <= 1.5)) "" else "  <- still travelling"))
  cat(sprintf("  reprojections after the freeze     %d %s\n", x$projections,
    if (x$projections == 0L) "" else "  <- outside the compact set"))
  cat(sprintf("  step halvings after the freeze     %d\n", x$backtracks))
  cat(sprintf("  refused steps after the freeze     %d %s\n", x$noMove,
    if (x$noMove == 0L) "" else "  <- proposals rejected as invalid"))
  if (is.finite(x$preconditionerConditioning))
    cat(sprintf("  preconditioner conditioning        %.3g\n",
      x$preconditionerConditioning))
  cat(sprintf("  verdict                            %s\n",
    if (isTRUE(x$settled)) "settled" else "NOT settled"))
  invisible(x)
}
