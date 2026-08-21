## Does the copula-SAEM estimator still target the EXACT marginal likelihood?
##
## Testing the SAEM fixed point directly would confound the IFM approximation
## with SAEM's own Monte Carlo noise.  So the ESTIMATOR DEFINITION is isolated:
## on one dataset, compute
##   theta_IFM  the fixed point of the IFM estimating equations the M-step solves
##   theta_ML   the direct maximiser of the exact marginal likelihood
## both with common random numbers, and compare.  propErr is FIXED at the truth
## in both, so the comparison is over (mu1, mu2, sd1, sd2, tau) only.
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R"); source("R/simEta.R"); source("R/simpleModels.R"); source("llcheck_core.R")

a <- commandArgs(trailingOnly = TRUE)
NREP <- if (length(a) > 0) as.integer(a[1]) else 3
N    <- if (length(a) > 1) as.integer(a[2]) else 100
M    <- if (length(a) > 2) as.integer(a[3]) else 1500

m <- MODELS$iv1; d <- 2; PE <- 0.10
tw <- makeTwinVines(etaR(d), "gumbel")
FAM <- etaVineFlat(tw$alt)[[1]]$family; ROT <- etaVineFlat(tw$alt)[[1]]$rotation
TAU_TRUE <- tw$tauAlt[1]
cat(sprintf("truth: mu=%s sd=%s tau=%.4f family=%s rot=%d\n",
            paste(m$true, collapse = ","), ETA_SD, TAU_TRUE, FAM, ROT))

llAt <- function(p, dat, seed) {
  q <- unpack(p)
  isPass(dat, q$mu, q$sdv, vine2FromTau(q$tau, FAM, ROT), PE, m$f,
         M = M, seed = seed, keepDraws = FALSE)$ll
}

## ---- the IFM M-step, exactly as copulaMstep() does it -----------------------
## mu absorbs the posterior mean of eta; sd is the second moment; the copula is
## fitted to the PIT with that sd -- the sd block ignoring the copula's own
## dependence on sd is precisely approximation (A).
ifmStep <- function(p, dat, seed) {
  q <- unpack(p)
  vn <- vine2FromTau(q$tau, FAM, ROT)
  r <- isPass(dat, q$mu, q$sdv, vn, PE, m$f, M = M, seed = seed, keepDraws = TRUE)
  mbar <- colMeans(r$m1)
  sdNew <- sqrt(pmax(colMeans(r$m2) - mbar^2, 1e-8))
  U <- do.call(rbind, lapply(seq_along(r$E), function(k) pnorm(sweep(r$E[[k]], 2, mbar, "-") %*% diag(1/sdNew))))
  W <- unlist(r$W) / length(r$E)
  U <- pmin(pmax(U, 1e-6), 1 - 1e-6)
  negpl <- function(z) {
    v <- vine2FromTau(tanh(z), FAM, ROT)
    -sum(W * log(pmax(rvinecopulib::dvinecop(U, v), 1e-300)))
  }
  zopt <- optimize(negpl, c(-2.5, 2.5))$minimum
  list(p = pack(q$mu + mbar, sdNew, tanh(zopt)), ll = r$ll, essMin = r$essMin)
}

## ---- the dropped score term (A), evaluated directly -------------------------
## dQ/d sd_j = E[-1/sd_j + eta_j^2/sd_j^3]  +  D_j,   D_j = -(1/sd_j) E[z_j phi(z_j) dlogc/du_j]
## The algebra says E_prior[D_j] = 0 exactly, so by the tower property D_j has
## mean zero at the TRUTH -- consistent but not the efficient score.
droppedTerm <- function(p, dat, seed) {
  q <- unpack(p); vn <- vine2FromTau(q$tau, FAM, ROT)
  r <- isPass(dat, q$mu, q$sdv, vn, PE, m$f, M = M, seed = seed, keepDraws = TRUE)
  Dj <- numeric(d); Mj <- numeric(d)
  for (j in seq_len(d)) {
    acc <- accM <- 0
    for (k in seq_along(r$E)) {
      E <- r$E[[k]]; w <- r$W[[k]] / length(r$E)
      z <- E[, j] / q$sdv[j]
      u <- pnorm(sweep(E, 2, q$sdv, "/")); u <- pmin(pmax(u, 1e-8), 1 - 1e-8)
      h <- 1e-5; up <- u; um <- u
      up[, j] <- pmin(u[, j] + h, 1 - 1e-9); um[, j] <- pmax(u[, j] - h, 1e-9)
      dlc <- (log(rvinecopulib::dvinecop(up, vn)) - log(rvinecopulib::dvinecop(um, vn))) / (up[, j] - um[, j])
      acc  <- acc  + sum(w * (-(1 / q$sdv[j]) * z * dnorm(z) * dlc))
      accM <- accM + sum(w * (-1 / q$sdv[j] + E[, j]^2 / q$sdv[j]^3))
    }
    Dj[j] <- acc; Mj[j] <- accM
  }
  list(dropped = Dj, moment = Mj)
}

out <- list()
for (rep in seq_len(NREP)) {
  set.seed(100 + rep); s <- simModel(m, N, tw$alt)
  dat <- s$data
  pTrue <- pack(log(m$true), rep(ETA_SD, d), TAU_TRUE)

  ## (a) dropped score term at the TRUTH -- tests the algebra
  dt <- droppedTerm(pTrue, dat, seed = 11)
  cat(sprintf("[rep %d] at TRUTH  dropped D=%s   moment-part=%s\n", rep,
              paste(sprintf("%+.3f", dt$dropped), collapse = ","),
              paste(sprintf("%+.3f", dt$moment), collapse = ",")))
  flush.console()

  ## (b) IFM fixed point
  p <- pTrue; tr <- list()
  for (it in 1:25) {
    st <- ifmStep(p, dat, seed = 7)
    delta <- max(abs(st$p - p)); p <- st$p; tr[[it]] <- c(p, delta)
    if (delta < 1e-4) break
  }
  pIFM <- p
  cat(sprintf("[rep %d] IFM converged in %d iters, delta=%.2e\n", rep, it, delta)); flush.console()

  ## (c) direct maximiser of the exact marginal likelihood, CRN
  op <- optim(pIFM, function(pp) -llAt(pp, dat, seed = 7), method = "Nelder-Mead",
              control = list(maxit = 400, reltol = 1e-8))
  pML <- op$par
  cat(sprintf("[rep %d] direct ML done, ll=%.4f (IFM ll=%.4f)\n", rep, -op$value,
              llAt(pIFM, dat, seed = 7))); flush.console()

  ## (d) gradient of the exact LL at both points, error bars over IS seeds
  gradAt <- function(pp, seeds, h = 0.02) {
    G <- sapply(seeds, function(sd_) sapply(seq_along(pp), function(j) {
      pu <- pp; pl <- pp; pu[j] <- pu[j] + h; pl[j] <- pl[j] - h
      (llAt(pu, dat, sd_) - llAt(pl, dat, sd_)) / (2 * h)
    }))
    list(mean = rowMeans(G), sd = apply(G, 1, sd))
  }
  seeds <- c(21, 22, 23)
  gI <- gradAt(pIFM, seeds); gM <- gradAt(pML, seeds)

  nmv <- c("mu1", "mu2", "logSd1", "logSd2", "atanhTau")
  cat(sprintf("[rep %d] grad at IFM : %s\n", rep,
      paste(sprintf("%s=%+.2f(%.2f)", nmv, gI$mean, gI$sd), collapse = " ")))
  cat(sprintf("[rep %d] grad at ML  : %s\n", rep,
      paste(sprintf("%s=%+.2f(%.2f)", nmv, gM$mean, gM$sd), collapse = " ")))
  cat(sprintf("[rep %d] pIFM=%s\n[rep %d] pML =%s\n", rep,
      paste(sprintf("%.4f", pIFM), collapse = ","), rep,
      paste(sprintf("%.4f", pML), collapse = ","))); flush.console()

  out[[rep]] <- list(rep = rep, pTrue = pTrue, pIFM = pIFM, pML = pML,
                     gradIFM = gI, gradML = gM, dropped = dt,
                     llIFM = llAt(pIFM, dat, 7), llML = -op$value, N = N, M = M)
}
saveRDS(out, "out/llcheck_main.rds")
cat("saved out/llcheck_main.rds\n")
