## Direct test of the integration-by-parts identity behind approximation (A).
##
## Claim: for eta ~ the copula prior itself (Gaussian margins sd_j, vine c),
##   E[ z_j phi(z_j) * dlog c/du_j ]  =  -E[ 1 - z_j^2 ]  =  0,
## because d/du [ z phi(z) ] = 1 - z^2 with z = Phi^{-1}(u), and z_j ~ N(0,1)
## marginally under the copula.  Hence the term the sd-update DROPS,
##   D_j = -(1/sd_j) E[ z_j phi(z_j) dlog c/du_j ],
## has mean ZERO at the true parameter -- the sd moment equation is an unbiased
## estimating function, so the estimator is consistent though not efficient.
##
## This is pure mathematics, no SAEM, so it can be checked to Monte Carlo
## precision at large n.
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R"); source("R/simEta.R")

dlogc_du <- function(u, vine, j, h = 1e-5) {
  up <- u; um <- u
  up[, j] <- pmin(u[, j] + h, 1 - 1e-12); um[, j] <- pmax(u[, j] - h, 1e-12)
  (log(rvinecopulib::dvinecop(up, vine)) - log(rvinecopulib::dvinecop(um, vine))) /
    (up[, j] - um[, j])
}

check <- function(lbl, vine, sdv, n = 400000, seed = 5) {
  set.seed(seed)
  eta <- rEtaVine(n, vine, sdv)
  d <- ncol(eta)
  Z <- sweep(eta, 2, sdv, "/")
  U <- pmin(pmax(pnorm(Z), 1e-9), 1 - 1e-9)
  for (j in seq_len(d)) {
    g <- Z[, j] * dnorm(Z[, j])
    dl <- dlogc_du(U, vine, j)
    lhs <- g * dl                 # should average to 0
    rhs <- 1 - Z[, j]^2           # the identity says E[lhs] = -E[rhs]
    cat(sprintf("%-22s j=%d  E[z phi(z) dlogc]=%+.5f (se %.5f)   -E[1-z^2]=%+.5f (se %.5f)   D_j=%+.5f (se %.5f)\n",
        lbl, j, mean(lhs), sd(lhs) / sqrt(n), -mean(rhs), sd(rhs) / sqrt(n),
        -mean(lhs) / sdv[j], sd(lhs) / sqrt(n) / sdv[j]))
  }
}

d <- 2; sdv <- c(0.30, 0.30)
tw2 <- makeTwinVines(0.55^abs(outer(1:2, 1:2, "-")), "gumbel")
check("d=2 Gaussian vine", tw2$gauss, sdv)
check("d=2 Gumbel vine",   tw2$alt,   sdv)
tw2c <- makeTwinVines(0.55^abs(outer(1:2, 1:2, "-")), "clayton")
check("d=2 Clayton vine",  tw2c$alt,  sdv)

d4 <- 0.55^abs(outer(1:4, 1:4, "-")); sd4 <- rep(0.30, 4)
tw4 <- makeTwinVines(d4, "gumbel")
check("d=4 Gumbel vine",   tw4$alt,   sd4, n = 200000)

## And the counterfactual: if the margins are NOT the ones the PIT assumes
## (sd misspecified), the identity should FAIL -- confirming that it is the
## correctness of the marginal spec at the truth that buys unbiasedness.
cat("\n-- counterfactual: PIT uses a WRONG sd (0.35 instead of 0.30) --\n")
set.seed(5); eta <- rEtaVine(400000, tw2$alt, c(0.30, 0.30))
Zw <- sweep(eta, 2, c(0.35, 0.35), "/")
Uw <- pmin(pmax(pnorm(Zw), 1e-9), 1 - 1e-9)
for (j in 1:2) {
  g <- Zw[, j] * dnorm(Zw[, j]); dl <- dlogc_du(Uw, tw2$alt, j)
  cat(sprintf("wrong-sd PIT j=%d  E[z phi(z) dlogc]=%+.5f (se %.5f)\n",
      j, mean(g * dl), sd(g * dl) / sqrt(400000)))
}
