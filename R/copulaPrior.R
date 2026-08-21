## Copula prior on the random effects (research extension).
##
## MODEL
##   eta_j ~ N(0, sd_j^2)                      marginals stay Gaussian
##   (u_1..u_d) ~ R-vine copula                dependence is the new part
##   u_j = pnorm(eta_j / sd_j)
##
## With Gaussian pair copulas everywhere this is EXACTLY MVN(0, Omega)
## (Bedford-Cooke), so the stock saemix model is the nested null.  estep() and
## mstep() take the copula path only when copulaSet() has been called; otherwise
## they are bit-identical to stock saemix.
##
## WHAT IS AND IS NOT JUSTIFIED
## The complete-data likelihood is NOT in the curved exponential family once any
## pair copula is non-Gaussian, so assumption (M1) of Delyon et al. (1999) fails
## and there is no finite sufficient statistic.  Q-hat is instead carried as a
## WEIGHTED PARTICLE POOL -- which is Delyon's general recursion (their eq. 6)
## with Q-hat represented by a weighted empirical measure.  Correct algorithm,
## no theorem.  The M-step is IFM (margins by weighted moments, copula by
## weighted MLE), not the exact joint maximiser; convergence is checked
## empirically, never assumed.
##
## ORDERING CONTRACT: vine variable j == eta column j == varList$ind.eta[j].

.cop <- new.env(parent = emptyenv())

copulaClear <- function() rm(list = ls(.cop), envir = .cop)

#' @param vine       rvinecopulib vinecop_dist on d = length(ind.eta) variables
#' @param familySet  families the M-step may select from; NULL = keep structure
#'                   and families fixed, refit parameters only
#' @param poolMax    max SA particle-pool size (iterations retained)
#' @param refitEvery refit the vine every k iterations (the M-step is the
#'                   expensive part; k>1 is a deliberate ECM-style delay)
copulaSet <- function(vine, sd, familySet = NULL, poolMax = 40L,
                      refitEvery = 1L, fitFrom = 1L, verbose = FALSE,
                      freezeSd = FALSE, freezeVine = FALSE,
                      mode = c("pool", "sa", "joint"), sdFromSS = TRUE, truncLvl = Inf) {
  mode <- match.arg(mode)
  copulaClear()
  .cop$vine <- vine
  .cop$sd <- sd
  .cop$sdPrev <- sd          # anneal relative to the INITIAL sd, as stock does to omega.init
  .cop$d <- length(sd)
  .cop$familySet <- familySet
  .cop$poolMax <- poolMax
  .cop$refitEvery <- refitEvery
  .cop$fitFrom <- fitFrom
  .cop$verbose <- verbose
  .cop$freezeSd <- freezeSd
  .cop$freezeVine <- freezeVine
  .cop$mode <- mode
  .cop$jointMaxit <- 40L
  .cop$sdFromSS <- sdFromSS
  .cop$truncLvl <- truncLvl
  .cop$famFixed <- NULL
  .cop$poolEta <- NULL
  .cop$poolW <- NULL
  .cop$trace <- list()
  .cop$diag <- list()
  invisible(TRUE)
}

copulaActive <- function() !is.null(.cop$vine)
copulaGet <- function() as.list(.cop)

## -log p(eta) up to an additive constant, matching saemix's U.eta convention.
copulaUeta <- function(etaM) {
  eta <- as.matrix(etaM)
  if (ncol(eta) != .cop$d)
    stop("copula prior: eta has ", ncol(eta), " columns, vine has ", .cop$d)
  sdv <- .cop$sd
  u <- pnorm(sweep(eta, 2, sdv, "/"))
  u <- pmin(pmax(u, 1e-10), 1 - 1e-10)
  lm <- rowSums(dnorm(sweep(eta, 2, sdv, "/"), log = TRUE)) - sum(log(sdv))
  lc <- log(pmax(rvinecopulib::dvinecop(u, .cop$vine), 1e-300))
  -(lm + lc)
}

## Draw n etas from the prior (kernel 1's independence proposal).
copulaRandEta <- function(n) {
  u <- rvinecopulib::rvinecop(n, .cop$vine)
  sweep(qnorm(u), 2, .cop$sd, "*")
}

## Copula-implied covariance of eta.  Used only where a Gaussian surrogate is
## still needed: the betas GLS block and the random-walk step scaling.
##
## eta_j = sd_j * qnorm(u_j), so the implied CORRELATION depends on the vine
## alone, not on sd.  Cache it per vine (refit is rare) and rescale for free;
## fixed CRN keeps the map deterministic so no simulation noise enters the loop.
copulaImpliedCor <- function(nsim = 20000L) {
  if (!is.null(.cop$corVine) && identical(.cop$corVine, .cop$vine)) return(.cop$corCache)
  ## Fixed CRN keeps sd -> Omega deterministic, but set.seed() here would reset
  ## R's GLOBAL stream inside the SAEM loop and make the E-step replay identical
  ## draws every iteration.  Save and restore around it.
  R <- withSeed(9781L, stats::cor(qnorm(rvinecopulib::rvinecop(nsim, .cop$vine))))
  .cop$corVine <- .cop$vine; .cop$corCache <- R
  R
}

## Evaluate expr under a fixed seed WITHOUT disturbing the caller's RNG stream.
withSeed <- function(seed, expr) {
  hasOld <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (hasOld) old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit({
    if (hasOld) assign(".Random.seed", old, envir = globalenv())
    else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
      rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  set.seed(seed)
  expr
}

copulaOmega <- function() {
  R <- copulaImpliedCor()
  diag(.cop$sd) %*% R %*% diag(.cop$sd)
}

## ---- SA particle pool: Delyon eq. (6) with a weighted empirical measure -----
## poolW carries the stochastic-approximation forgetting factor; a draw enters
## with weight gamma_k/nchains and every older draw is scaled by (1 - gamma_k).
copulaPoolUpdate <- function(etaM, gamma, nchains) {
  .cop$curEta <- etaM
  if (identical(.cop$mode, "sa")) return(invisible(NULL))   # SA mode keeps no pool
  n <- nrow(etaM)
  wNew <- rep(gamma / nchains, n)
  if (is.null(.cop$poolEta)) {
    .cop$poolEta <- etaM; .cop$poolW <- wNew
  } else {
    .cop$poolEta <- rbind(.cop$poolEta * 1, etaM)
    .cop$poolW <- c(.cop$poolW * (1 - gamma), wNew)
  }
  ## prune negligible particles and cap the pool
  keep <- .cop$poolW > max(.cop$poolW) * 1e-6
  if (sum(keep) > .cop$poolMax * n) {
    ord <- order(.cop$poolW, decreasing = TRUE)[seq_len(.cop$poolMax * n)]
    keep <- seq_along(.cop$poolW) %in% ord
  }
  .cop$poolEta <- .cop$poolEta[keep, , drop = FALSE]
  .cop$poolW <- .cop$poolW[keep]
  invisible(NULL)
}

## ---- M-step (IFM): margins by weighted moments, copula by weighted MLE ------
etaVineFromFlat <- function(flat, d, structure) {
  k <- 1L; pcs <- vector("list", d - 1)
  for (t in seq_len(d - 1)) {
    ed <- vector("list", d - t)
    for (e in seq_len(d - t)) { ed[[e]] <- flat[[k]]; k <- k + 1L }
    pcs[[t]] <- ed
  }
  rvinecopulib::vinecop_dist(pcs, structure)
}


## ---- joint (MLE) M-step helpers --------------------------------------------
##
## The IFM M-step solves M_j = 0 for sd_j (the Gaussian-marginal score) and
## drops D_j, the copula's dependence on sd through u = pnorm(eta/sd).  E[D_j]
## = 0 at the truth, so IFM is consistent -- but it solves an unbiased
## ESTIMATING EQUATION, not the score, so the fixed point is not the MLE.
## Consequences: no EM ascent guarantee, and standard errors / likelihood-ratio
## tests from the observed information are invalid (they need the sandwich).
##
## mode = "joint" maximises Q over (sd, pair-copula taus) TOGETHER, restoring
## MLE status for the copula model.  Families and structure stay frozen after
## the first selection, so Q is a fixed function of a fixed-length parameter
## vector -- which is also what Delyon's (M1)/(M5) require.

copulaPadFlat <- function(vine, d) {
  f <- unlist(vine$pair_copulas, recursive = FALSE)
  nFull <- d * (d - 1) / 2
  if (length(f) < nFull)
    f <- c(f, rep(list(rvinecopulib::bicop_dist("indep")), nFull - length(f)))
  f
}

## Admissible Kendall-tau range for a fixed family+rotation, unconstrained.
## Archimedean families are one-sided: rotation 0/180 -> tau > 0, 90/270 -> < 0.
copulaTauTwoSided <- function(fam) fam %in% c("gaussian", "frank", "t", "indep")
copulaTauToX <- function(tau, fam, rot) {
  if (copulaTauTwoSided(fam)) return(atanh(max(min(tau / 0.98, 0.999), -0.999)))
  stats::qlogis(max(min(abs(tau) / 0.95, 0.999), 0.001))
}
copulaXToTau <- function(x, fam, rot) {
  if (copulaTauTwoSided(fam)) return(0.98 * tanh(x))
  (if (rot %in% c(90, 270)) -1 else 1) * 0.95 * stats::plogis(x)
}

copulaBicopFromTau <- function(tmpl, tau) {
  if (tmpl$family == "indep") return(tmpl)
  if (tmpl$family == "t")
    return(rvinecopulib::bicop_dist("t", rotation = tmpl$rotation,
             parameters = c(sin(pi * tau / 2), tmpl$parameters[2])))
  rvinecopulib::bicop_dist(tmpl$family, rotation = tmpl$rotation,
                           parameters = rvinecopulib::ktau_to_par(tmpl, tau))
}

copulaLogPrior <- function(E, vine, sdv) {
  z <- sweep(E, 2, sdv, "/")
  u <- pmin(pmax(pnorm(z), 1e-10), 1 - 1e-10)
  rowSums(dnorm(z, log = TRUE)) - sum(log(sdv)) +
    log(pmax(rvinecopulib::dvinecop(u, vine), 1e-300))
}

## Warm-started, limited-iteration joint maximisation.  An increase in Q is all
## a generalised-EM step needs; iterating across SAEM iterations converges to
## the maximiser (Delyon et al. 1999, Sec. 8.1).
copulaMaximiseJoint <- function(E, w, sd0, vine0, d, maxit = 40L) {
  flat <- copulaPadFlat(vine0, d)
  fams <- vapply(flat, function(b) b$family, character(1))
  rots <- vapply(flat, function(b) b$rotation, numeric(1))
  free <- which(fams != "indep")
  tau0 <- vapply(seq_along(flat), function(i)
    if (fams[i] == "indep") 0 else rvinecopulib::par_to_ktau(flat[[i]]), numeric(1))
  build <- function(tt) {
    fl <- lapply(seq_along(flat), function(i) copulaBicopFromTau(flat[[i]], tt[i]))
    k <- 1L; pcs <- vector("list", d - 1)
    for (t in seq_len(d - 1)) {
      ed <- vector("list", d - t)
      for (e in seq_len(d - t)) { ed[[e]] <- fl[[k]]; k <- k + 1L }
      pcs[[t]] <- ed
    }
    rvinecopulib::vinecop_dist(pcs, vine0$structure)
  }
  if (!length(free)) return(list(sd = sd0, vine = vine0, conv = NA_integer_))
  par0 <- c(log(sd0), vapply(free, function(i)
    copulaTauToX(tau0[i], fams[i], rots[i]), numeric(1)))
  negQ <- function(pv) {
    sdv <- exp(pv[seq_len(d)])
    if (any(!is.finite(sdv)) || any(sdv < 1e-6) || any(sdv > 1e3)) return(1e10)
    tt <- tau0
    tt[free] <- vapply(seq_along(free), function(k)
      copulaXToTau(pv[d + k], fams[free[k]], rots[free[k]]), numeric(1))
    v <- try(build(tt), silent = TRUE)
    if (inherits(v, "try-error")) return(1e10)
    val <- -sum(w * copulaLogPrior(E, v, sdv))
    if (!is.finite(val)) 1e10 else val
  }
  base <- negQ(par0)
  op <- try(stats::optim(par0, negQ, method = "BFGS", control = list(maxit = maxit)),
            silent = TRUE)
  if (inherits(op, "try-error") || !is.finite(op$value) || op$value >= base)
    return(list(sd = sd0, vine = vine0, conv = -1L))
  tt <- tau0
  tt[free] <- vapply(seq_along(free), function(k)
    copulaXToTau(op$par[d + k], fams[free[k]], rots[free[k]]), numeric(1))
  list(sd = exp(op$par[seq_len(d)]), vine = build(tt), conv = op$convergence)
}

copulaMstep <- function(kiter, nbiterSa = 0L, alpha1Sa = 1, sdSS = NULL, gamma = 1) {
  if (is.null(.cop$curEta)) return(invisible(NULL))
  if (kiter < .cop$fitFrom) return(invisible(NULL))
  if (kiter %% .cop$refitEvery != 0L) return(invisible(NULL))
  if (identical(.cop$mode, "joint")) {
    if (is.null(.cop$poolEta)) return(invisible(NULL))
    E <- .cop$poolEta; w <- .cop$poolW / sum(.cop$poolW)
    if (is.null(.cop$famFixed) && !isTRUE(.cop$freezeVine)) {
      u0 <- pmin(pmax(pnorm(sweep(E, 2, .cop$sd, "/")), 1e-6), 1 - 1e-6)
      fs0 <- if (is.null(.cop$familySet)) "parametric" else .cop$familySet
      f0 <- try(rvinecopulib::vinecop(u0, structure = .cop$vine$structure,
                                      family_set = fs0, weights = w * length(w),
                                      selcrit = "aic", trunc_lvl = .cop$truncLvl),
                silent = TRUE)
      if (!inherits(f0, "try-error")) {
        .cop$vine <- f0
        .cop$famFixed <- unique(vapply(unlist(f0$pair_copulas, recursive = FALSE),
                                       function(b) b$family, character(1)))
      }
    }
    jm <- copulaMaximiseJoint(E, w, .cop$sd, .cop$vine, .cop$d,
            maxit = if (is.null(.cop$jointMaxit)) 40L else .cop$jointMaxit)
    sdNew <- if (isTRUE(.cop$freezeSd)) .cop$sd else jm$sd
    if (kiter <= nbiterSa && !is.null(.cop$sdPrev))
      sdNew <- pmax(sdNew, .cop$sdPrev * sqrt(alpha1Sa))
    .cop$sdPrev <- sdNew; .cop$sd <- sdNew
    if (!isTRUE(.cop$freezeVine)) .cop$vine <- jm$vine
    if (isTRUE(.cop$verbose))
      .cop$trace[[length(.cop$trace) + 1L]] <- list(
        kiter = kiter, sd = sdNew, npool = nrow(E), conv = jm$conv,
        tau = vapply(copulaPadFlat(.cop$vine, .cop$d), function(b)
          if (b$family == "indep") 0 else rvinecopulib::par_to_ktau(b), numeric(1)),
        fam = vapply(copulaPadFlat(.cop$vine, .cop$d), function(b) b$family, character(1)))
    return(invisible(NULL))
  }
  if (identical(.cop$mode, "sa")) {
    E <- .cop$curEta
    w <- rep(1 / nrow(E), nrow(E))
  } else {
    if (is.null(.cop$poolEta)) return(invisible(NULL))
    E <- .cop$poolEta
    w <- .cop$poolW / sum(.cop$poolW)
  }

  ## block 1: marginal SDs (exact maximiser of the Gaussian marginal term;
  ## ignores the copula term's sd dependence -- this is the IFM approximation)
  ## The Gaussian margins STILL have an exact sufficient statistic even when the
  ## copula does not, so take sd from saemix's own SA-accumulated second moment
  ## rather than from the truncated pool.  Only the copula term needs the pool.
  sdNew <- if (isTRUE(.cop$freezeSd)) .cop$sd
           else if (isTRUE(.cop$sdFromSS) && !is.null(sdSS)) sdSS
           else sqrt(colSums(w * E^2))
  ## Simulated annealing, same rule saemix applies to diag(omega): during
  ## burn-in a variance may fall by at most a factor alpha1.sa per iteration.
  ## This is LOAD-BEARING, not cosmetic.  Without it a weakly identified eta
  ## (here V2) collapses: a tighter prior shrinks its posterior draws, which
  ## shrinks sd, which tightens the prior -- and the copula absorbs the lost
  ## marginal scale as near-perfect dependence (tau -> 0.84, implied rho 0.91).
  ## The vine parametrisation makes this worse than the plain covariance one,
  ## because the copula is fitted on normal SCORES and is therefore invariant to
  ## the marginal scale it is competing with.
  if (kiter <= nbiterSa && !is.null(.cop$sdPrev))
    sdNew <- pmax(sdNew, .cop$sdPrev * sqrt(alpha1Sa))
  .cop$sdPrev <- sdNew
  .cop$sd <- sdNew

  ## block 2: pair-copula parameters on the pinned structure
  u <- pnorm(sweep(E, 2, sdNew, "/"))
  u <- pmin(pmax(u, 1e-6), 1 - 1e-6)
  fs <- if (is.null(.cop$familySet)) unique(vapply(
          unlist(.cop$vine$pair_copulas, recursive = FALSE),
          function(b) b$family, character(1))) else .cop$familySet
  if (!isTRUE(.cop$freezeVine)) {
    if (identical(.cop$mode, "sa")) {
      ## SA-on-parameter: fit the vine to the CURRENT draws only, then take a
      ## Robbins-Monro step on the pair-copula parameters in Kendall-tau space.
      ## No pool, no truncation, O(1) memory -- the same optimise-then-damp
      ## pattern saemix already uses for the residual-error parameters.
      ## Families are selected once (at fitFrom) and then held, so that tau is
      ## being averaged within a fixed model rather than across changing ones.
      cur <- .cop$curEta
      uc <- pnorm(sweep(cur, 2, sdNew, "/"))
      uc <- pmin(pmax(uc, 1e-6), 1 - 1e-6)
      fsUse <- if (is.null(.cop$famFixed)) fs else .cop$famFixed
      fit <- try(rvinecopulib::vinecop(uc, structure = .cop$vine$structure,
                                       family_set = fsUse, selcrit = "aic",
                                       trunc_lvl = .cop$truncLvl), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        if (is.null(.cop$famFixed))
          .cop$famFixed <- unique(vapply(unlist(fit$pair_copulas, recursive = FALSE),
                                         function(b) b$family, character(1)))
        ## A truncated fit returns fewer trees than the full vine; pad the
        ## missing higher trees with independence so the flat edge list keeps
        ## its canonical length and order.
        padFlat <- function(v) {
          f <- unlist(v$pair_copulas, recursive = FALSE)
          nFull <- .cop$d * (.cop$d - 1) / 2
          if (length(f) < nFull)
            f <- c(f, rep(list(rvinecopulib::bicop_dist("indep")), nFull - length(f)))
          f
        }
        old <- padFlat(.cop$vine)
        new <- padFlat(fit)
        blended <- lapply(seq_along(new), function(i) {
          tOld <- if (old[[i]]$family == "indep") 0 else rvinecopulib::par_to_ktau(old[[i]])
          tNew <- if (new[[i]]$family == "indep") 0 else rvinecopulib::par_to_ktau(new[[i]])
          tt <- tOld + gamma * (tNew - tOld)
          b <- new[[i]]
          if (b$family == "indep") return(b)
          tt <- sign(tt) * min(abs(tt), 0.95)
          rvinecopulib::bicop_dist(b$family, rotation = b$rotation,
                                   parameters = rvinecopulib::ktau_to_par(b, tt))
        })
        .cop$vine <- etaVineFromFlat(blended, .cop$d, .cop$vine$structure)
      }
    } else {
      ## Freeze families after the first selection, exactly as mode="sa" does.
      ## Re-selecting families at every refit changes the objective between
      ## iterations, so there is no fixed target for the algorithm to converge
      ## to -- the SA recursion is then chasing a moving function.
      fsUse <- if (is.null(.cop$famFixed)) fs else .cop$famFixed
      fit <- try(rvinecopulib::vinecop(u, structure = .cop$vine$structure,
                                       family_set = fsUse, weights = w * length(w),
                                       selcrit = "aic", trunc_lvl = .cop$truncLvl),
                 silent = TRUE)
      if (!inherits(fit, "try-error")) {
        if (is.null(.cop$famFixed))
          .cop$famFixed <- unique(vapply(unlist(fit$pair_copulas, recursive = FALSE),
                                         function(b) b$family, character(1)))
        .cop$vine <- fit
      }
    }
  }
  if (isTRUE(.cop$verbose))
    .cop$trace[[length(.cop$trace) + 1L]] <- list(
      kiter = kiter, sd = sdNew, npool = nrow(E),
      tau = vapply(unlist(.cop$vine$pair_copulas, recursive = FALSE),
                   function(b) if (b$family == "indep") 0 else
                     rvinecopulib::par_to_ktau(b), numeric(1)),
      fam = vapply(unlist(.cop$vine$pair_copulas, recursive = FALSE),
                   function(b) b$family, character(1)))
  invisible(NULL)
}

## Extra per-iteration diagnostics, only when verbose.  etaMean is the smoking
## gun for a mean/eta ridge: it must stay at 0.
copulaDiag <- function(kiter, betas, etaMean, pres, gamma) {
  if (!isTRUE(.cop$verbose)) return(invisible(NULL))
  n <- length(.cop$diag) + 1L
  .cop$diag[[n]] <- list(kiter = kiter, betas = betas, etaMean = etaMean,
                         pres = pres, gamma = gamma, sd = .cop$sd)
  invisible(NULL)
}

## Path-agnostic per-iteration trace, so the copula and stock paths can be
## compared iteration by iteration under the same seed.
.trc <- new.env(parent = emptyenv())
saemixTraceReset <- function() { .trc$L <- list(); options(saemixTrace = TRUE) }
saemixTraceGet <- function() .trc$L
.saemixTracePush <- function(kiter, betas, omdiag, pres, gamma, statrese) {
  .trc$L[[length(.trc$L) + 1L]] <- list(kiter = kiter, betas = betas,
    omdiag = omdiag, pres = pres, gamma = gamma, statrese = statrese)
  invisible(NULL)
}
