## Canonical copula representation of the random-effect distribution.
##
## MODEL
##   eta_j ~ N(0, sd_j^2)                        marginals stay Gaussian
##   (u_1..u_d) ~ R-vine copula C(.; theta_c)    dependence is the only new part
##   u_j = pnorm(eta_j / sd_j)
##
## With Gaussian pair copulas at every edge this is EXACTLY MVN(0, Omega)
## (Bedford-Cooke), so the current saemix model is the nested null.
##
## ORDERING CONTRACT (load-bearing -- see tests/test-ordering.R)
##   vine variable j  ==  eta column j  ==  varList$ind.eta[j] in saemix.
##   Structure is always dvine_structure(1:d); pair_copulas[[t]][[e]] is the
##   copula of (e, e+t) given the variables strictly between them.
##   rvinecopulib is canonical.  VineCopula uses a TRANSPOSED array convention
##   and is used only as an independent oracle; never mix the two objects.

suppressMessages(library(rvinecopulib))

## ---- structure -------------------------------------------------------------

etaVineStructure <- function(d) rvinecopulib::dvine_structure(seq_len(d))

## Assemble a vinecop_dist from a flat list of bicop_dist in canonical order.
## `flat` is length d(d-1)/2, ordered tree 1 edges 1..d-1, tree 2 edges 1..d-2, ...
etaVine <- function(flat, d, order = seq_len(d)) {
  stopifnot(length(flat) == d * (d - 1) / 2)
  k <- 1L; pcs <- vector("list", d - 1)
  for (t in seq_len(d - 1)) {
    ed <- vector("list", d - t)
    for (e in seq_len(d - t)) { ed[[e]] <- flat[[k]]; k <- k + 1L }
    pcs[[t]] <- ed
  }
  rvinecopulib::vinecop_dist(pcs, rvinecopulib::dvine_structure(order))
}

## Relabel a canonical D-vine to the reversed variable order.  Same joint law,
## different inverse-Rosenblatt map -- the ordering trap, made explicit.
etaVineReversed <- function(flat, d) {
  k <- 1L; out <- list()
  for (t in seq_len(d - 1)) {
    n <- d - t
    out <- c(out, rev(flat[k:(k + n - 1)]))
    k <- k + n
  }
  etaVine(out, d, order = rev(seq_len(d)))
}

## Flat canonical edge list of a vinecop_dist built by etaVine().
etaVineFlat <- function(vine) unlist(vine$pair_copulas, recursive = FALSE)

## ---- partial correlations of a correlation matrix, in canonical edge order --

partialCor <- function(R, a, b, cond) {
  if (length(cond) == 0) return(R[a, b])
  idx <- c(a, b, cond)
  P <- solve(R[idx, idx, drop = FALSE])
  -P[1, 2] / sqrt(P[1, 1] * P[2, 2])
}

dvinePartialCor <- function(R) {
  d <- ncol(R)
  out <- numeric(0)
  for (t in seq_len(d - 1)) for (e in seq_len(d - t)) {
    cond <- if (t > 1) (e + 1):(e + t - 1) else integer(0)
    out <- c(out, partialCor(R, e, e + t, cond))
  }
  out
}

## Gaussian vine == MVN(0, R).  Validated exactly against dmvnorm in T1.
etaVineGaussian <- function(R) {
  pc <- dvinePartialCor(R)
  etaVine(lapply(pc, function(p) rvinecopulib::bicop_dist("gaussian", parameters = p)), ncol(R))
}

## ---- density and sampling on the eta scale ---------------------------------

etaToU <- function(eta, sd) {
  u <- pnorm(sweep(eta, 2, sd, "/"))
  ## rvinecopulib is unhappy at exactly 0/1
  pmin(pmax(u, 1e-10), 1 - 1e-10)
}

## log p(eta) = sum_j log dnorm(eta_j; 0, sd_j) + log c(u)
dEtaVine <- function(eta, vine, sd, log = TRUE) {
  eta <- as.matrix(eta)
  lm <- rowSums(dnorm(sweep(eta, 2, sd, "/"), log = TRUE)) - sum(log(sd))
  lc <- log(rvinecopulib::dvinecop(etaToU(eta, sd), vine))
  if (log) lm + lc else exp(lm + lc)
}

rEtaVine <- function(n, vine, sd) {
  u <- rvinecopulib::rvinecop(n, vine)
  sweep(qnorm(u), 2, sd, "*")
}

## Non-centered map eta = T(z), z ~ N(0, I_d), via the inverse Rosenblatt.
## ORDER-DEPENDENT by construction: a different variable order is a different T
## for the same joint law (test-ordering.R T5).  Pin it and never permute.
etaFromZ <- function(z, vine, sd) {
  u <- rvinecopulib::inverse_rosenblatt(pnorm(as.matrix(z)), vine)
  sweep(qnorm(u), 2, sd, "*")
}
zFromEta <- function(eta, vine, sd)
  qnorm(rvinecopulib::rosenblatt(etaToU(as.matrix(eta), sd), vine))
