## Route-C diagnostic: can a non-Gaussian eta dependence be SEEN from a
## standard Gaussian-eta SAEM fit?
##
## Statistic: 2*(loglik of a family-selected vine - loglik of the all-Gaussian
## vine) on the same D-vine structure, so only the pair-copula FAMILIES differ.
## The all-Gaussian vine is exactly MVN, i.e. the current saemix model, so this
## is a likelihood-ratio-shaped statistic against the correct nested null.
##
## Its null distribution is NOT chi-square here (family selection + boundary +
## within-subject draw correlation), so it is calibrated empirically from
## replicates simulated under the Gaussian truth.  That empirical calibration is
## what makes the pooled-draw dependence harmless.

source("R/etaCopula.R")

VINE_FAMSET <- c("gaussian", "t", "clayton", "gumbel", "frank", "joe")

## u from etas.  "param" = parametric PIT with Gaussian margins (matches the
## model, and is what an M-step would do).  "rank" = pseudo-observations
## (immune to marginal misspecification; pure dependence).
etaToUnif <- function(eta, sd = NULL, how = c("param", "rank")) {
  how <- match.arg(how)
  eta <- as.matrix(eta)
  if (how == "rank") return(rvinecopulib::pseudo_obs(eta))
  if (is.null(sd)) sd <- apply(eta, 2, stats::sd)
  etaToU(eta, sd)
}

## Fit both vines on the pinned D-vine structure and return the statistic.
vineLRT <- function(u) {
  d  <- ncol(u)
  st <- etaVineStructure(d)
  f0 <- rvinecopulib::vinecop(u, structure = st, family_set = "gaussian")
  f1 <- rvinecopulib::vinecop(u, structure = st, family_set = VINE_FAMSET,
                              selcrit = "aic")
  fams <- vapply(unlist(f1$pair_copulas, recursive = FALSE),
                 function(b) b$family, character(1))
  rots <- vapply(unlist(f1$pair_copulas, recursive = FALSE),
                 function(b) b$rotation, numeric(1))
  list(stat    = 2 * (f1$loglik - f0$loglik),
       nNonGauss = sum(!fams %in% c("gaussian", "indep")),
       families = paste(ifelse(rots == 0, fams, paste0(fams, rots)), collapse = ","),
       ll0 = f0$loglik, ll1 = f1$loglik, npars1 = f1$npars, npars0 = f0$npars)
}

## Empirical upper/lower tail-dependence proxies on a pair, at level q.
tailProxy <- function(u, i, j, q = 0.95) {
  up <- mean(u[, i] > q & u[, j] > q) / (1 - q)
  lo <- mean(u[, i] < 1 - q & u[, j] < 1 - q) / (1 - q)
  c(upper = up, lower = lo)
}

## Conditional draws from a saemix fit, as an N x d x nsamp array of ETAs
## (phi minus its individual mean), restricted to the columns with IIV.
condEtaDraws <- function(fit, nsamp = 5) {
  fit <- saemix::conddist.saemix(fit, nsamp = nsamp, plot = FALSE)
  phi <- fit@results@phi.samp                       # N x nbpar x nsamp
  mp  <- fit@results@mean.phi                       # N x nbpar
  idx <- which(diag(fit@results@omega) > 1e-10)
  arr <- array(0, dim = c(dim(phi)[1], length(idx), dim(phi)[3]))
  for (s in seq_len(dim(phi)[3])) arr[, , s] <- phi[, idx, s] - mp[, idx]
  list(draws = arr, ebe = fit@results@cond.mean.phi[, idx, drop = FALSE] - mp[, idx, drop = FALSE],
       sd = sqrt(diag(fit@results@omega)[idx]), idx = idx, fit = fit)
}

## Average the statistic over the nsamp slices: each slice is one draw per
## subject (N rows), so n stays at N and the within-subject mixture is avoided.
vineLRTdraws <- function(arr, sd = NULL, how = "param") {
  res <- lapply(seq_len(dim(arr)[3]), function(s)
    vineLRT(etaToUnif(arr[, , s], sd, how)))
  list(stat = mean(vapply(res, function(r) r$stat, numeric(1))),
       nNonGauss = mean(vapply(res, function(r) r$nNonGauss, numeric(1))),
       families = res[[1]]$families,
       perSlice = vapply(res, function(r) r$stat, numeric(1)))
}
