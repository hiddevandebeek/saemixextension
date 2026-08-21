## Simulate etas from a Gaussian null and from a non-Gaussian R-vine TWIN that
## is matched to it edge-by-edge on Kendall's tau, with identical margins.
##
## The matching is the whole point: a Gaussian-eta SAEM sees the same marginal
## variances and (near) the same rank correlations in both arms, so anything
## that separates them is dependence SHAPE, not dependence STRENGTH.

source("R/etaCopula.R")

## family: rvinecopulib name ("gumbel", "clayton", "t", "frank", "joe")
## Negative-tau edges get the 90-degree rotation.
makeTwinVines <- function(R, family = "gumbel", dfT = 4) {
  d  <- ncol(R)
  pc <- dvinePartialCor(R)
  gaussFlat <- lapply(pc, function(p) rvinecopulib::bicop_dist("gaussian", parameters = p))
  tau <- vapply(gaussFlat, rvinecopulib::par_to_ktau, numeric(1))
  altFlat <- lapply(tau, function(tt) {
    if (abs(tt) < 1e-10) return(rvinecopulib::bicop_dist("indep"))
    if (family == "t") {
      rho <- sin(pi * tt / 2)
      return(rvinecopulib::bicop_dist("t", parameters = c(rho, dfT)))
    }
    rot  <- if (tt > 0) 0 else 90
    tmpl <- rvinecopulib::bicop_dist(family, rotation = rot,
                                     parameters = rvinecopulib::ktau_to_par(family, 0.5))
    rvinecopulib::bicop_dist(family, rotation = rot,
                             parameters = rvinecopulib::ktau_to_par(tmpl, tt))
  })
  list(gauss  = etaVine(gaussFlat, d),
       alt    = etaVine(altFlat,   d),
       tau    = tau,
       pc     = pc,
       tauAlt = vapply(altFlat, function(b)
                  if (b$family == "indep") 0 else rvinecopulib::par_to_ktau(b), numeric(1)))
}

## Draw N etas with given marginal SDs.
simEtaVine <- function(N, vine, sd) rEtaVine(N, vine, sd)
