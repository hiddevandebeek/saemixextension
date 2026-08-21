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
                      freezeSd = FALSE, freezeVine = FALSE) {
  copulaClear()
  .cop$vine <- vine
  .cop$sd <- sd
  .cop$d <- length(sd)
  .cop$familySet <- familySet
  .cop$poolMax <- poolMax
  .cop$refitEvery <- refitEvery
  .cop$fitFrom <- fitFrom
  .cop$verbose <- verbose
  .cop$freezeSd <- freezeSd
  .cop$freezeVine <- freezeVine
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
copulaImpliedCor <- function(nsim = 60000L) {
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
copulaMstep <- function(kiter) {
  if (is.null(.cop$poolEta)) return(invisible(NULL))
  if (kiter < .cop$fitFrom) return(invisible(NULL))
  if (kiter %% .cop$refitEvery != 0L) return(invisible(NULL))
  E <- .cop$poolEta; w <- .cop$poolW / sum(.cop$poolW)

  ## block 1: marginal SDs (exact maximiser of the Gaussian marginal term;
  ## ignores the copula term's sd dependence -- this is the IFM approximation)
  sdNew <- if (isTRUE(.cop$freezeSd)) .cop$sd else sqrt(colSums(w * E^2))
  .cop$sd <- sdNew

  ## block 2: pair-copula parameters on the pinned structure
  u <- pnorm(sweep(E, 2, sdNew, "/"))
  u <- pmin(pmax(u, 1e-6), 1 - 1e-6)
  fs <- if (is.null(.cop$familySet)) unique(vapply(
          unlist(.cop$vine$pair_copulas, recursive = FALSE),
          function(b) b$family, character(1))) else .cop$familySet
  if (!isTRUE(.cop$freezeVine)) {
    fit <- try(rvinecopulib::vinecop(u, structure = .cop$vine$structure,
                                     family_set = fs, weights = w * length(w),
                                     selcrit = "aic"), silent = TRUE)
    if (!inherits(fit, "try-error")) .cop$vine <- fit
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
