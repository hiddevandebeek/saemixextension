## Does the JOINT (MLE) M-step make the fixed point stationary for the EXACT
## marginal likelihood?  Same harness as llcheck_gradN.R, but the M-step
## maximises Q over (mu-shift, log sd, tau) TOGETHER instead of the IFM split.
##
## Why a mu-shift is part of it: the draws E are etas relative to the current
## mu, so moving mu to mu+delta relabels the same phi as eta = E - delta.  The
## OBSERVATION term is therefore invariant and maximising over delta using the
## same weighted draws is exact.
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R"); source("R/simEta.R"); source("R/simpleModels.R"); source("llcheck_core.R")

a <- commandArgs(trailingOnly = TRUE)
NREP <- if (length(a) > 0) as.integer(a[1]) else 3
N    <- if (length(a) > 1) as.integer(a[2]) else 100
M    <- if (length(a) > 2) as.integer(a[3]) else 1500
STEP <- if (length(a) > 3) a[4] else "joint"

m <- MODELS$iv1; d <- 2; PE <- 0.10
tw <- makeTwinVines(etaR(d), "gumbel")
FAM <- etaVineFlat(tw$alt)[[1]]$family; ROT <- etaVineFlat(tw$alt)[[1]]$rotation
TAU_TRUE <- tw$tauAlt[1]
nmv <- c("mu1", "mu2", "logSd1", "logSd2", "atanhTau")
oneSided <- !(FAM %in% c("gaussian", "frank", "t", "indep"))

llAt <- function(p, dat, seed) {
  q <- unpack(p)
  isPass(dat, q$mu, q$sdv, vine2FromTau(q$tau, FAM, ROT), PE, m$f,
         M = M, seed = seed, keepDraws = FALSE)$ll
}

## IFM step (the incumbent), reproduced for a matched comparison.
ifmStep <- function(p, dat, seed) {
  q <- unpack(p); vn <- vine2FromTau(q$tau, FAM, ROT)
  r <- isPass(dat, q$mu, q$sdv, vn, PE, m$f, M = M, seed = seed, keepDraws = TRUE)
  mbar <- colMeans(r$m1); sdNew <- sqrt(pmax(colMeans(r$m2) - mbar^2, 1e-8))
  U <- do.call(rbind, lapply(seq_along(r$E), function(k)
        pnorm(sweep(sweep(r$E[[k]], 2, mbar, "-"), 2, sdNew, "/"))))
  W <- unlist(r$W) / length(r$E)
  U <- pmin(pmax(U, 1e-6), 1 - 1e-6)
  z <- optimize(function(z) -sum(W * log(pmax(rvinecopulib::dvinecop(U, vine2FromTau(tanh(z), FAM, ROT)), 1e-300))),
                c(-2.5, 2.5))$minimum
  pack(q$mu + mbar, sdNew, tanh(z))
}

## JOINT step: one maximisation of Q over everything.
jointStep <- function(p, dat, seed) {
  q <- unpack(p); vn <- vine2FromTau(q$tau, FAM, ROT)
  r <- isPass(dat, q$mu, q$sdv, vn, PE, m$f, M = M, seed = seed, keepDraws = TRUE)
  EE <- do.call(rbind, r$E); WW <- unlist(r$W) / length(r$E)
  negQ <- function(v) {
    sdv <- exp(v[3:4]); tau <- tanh(v[5])
    if (any(!is.finite(sdv)) || any(sdv < 1e-6)) return(1e10)
    if (oneSided && tau <= 1e-3) return(1e10)
    if (abs(tau) > 0.95) return(1e10)
    val <- -sum(WW * dEtaVine(sweep(EE, 2, v[1:2], "-"),
                              vine2FromTau(tau, FAM, ROT), sdv, log = TRUE))
    if (!is.finite(val)) 1e10 else val
  }
  v0 <- c(0, 0, log(q$sdv), atanh(q$tau))
  op <- optim(v0, negQ, method = "BFGS", control = list(maxit = 100))
  if (op$value >= negQ(v0)) return(p)
  pack(q$mu + op$par[1:2], exp(op$par[3:4]), tanh(op$par[5]))
}
STEPFUN <- if (STEP == "joint") jointStep else ifmStep

res <- list()
for (rep in seq_len(NREP)) {
  set.seed(100 + rep); dat <- simModel(m, N, tw$alt)$data
  p <- pack(log(m$true), rep(ETA_SD, d), TAU_TRUE)
  for (it in 1:8) { pn <- STEPFUN(p, dat, 7); dl <- max(abs(pn - p)); p <- pn; if (dl < 1e-4) break }
  seeds <- c(31, 32, 33); h <- 0.02
  G <- sapply(seeds, function(sd_) sapply(seq_along(p), function(j) {
    pu <- p; pl <- p; pu[j] <- pu[j] + h; pl[j] <- pl[j] - h
    (llAt(pu, dat, sd_) - llAt(pl, dat, sd_)) / (2 * h) }))
  H <- sapply(seq_along(p), function(j) {
    pu <- p; pl <- p; pu[j] <- pu[j] + h; pl[j] <- pl[j] - h
    (llAt(pu, dat, 31) - 2 * llAt(p, dat, 31) + llAt(pl, dat, 31)) / h^2 })
  g <- rowMeans(G); se <- apply(G, 1, sd) / sqrt(length(seeds))
  disp <- ifelse(H < 0, -g / H, NA)
  cat(sprintf("[%s rep %d] fixed point = %s (iters %d)\n", STEP, rep,
      paste(sprintf("%.4f", p), collapse = ","), it))
  for (j in seq_along(p))
    cat(sprintf("   %-9s grad=%+9.3f se=%7.3f |t|=%6.2f  shift=%+.4f\n",
        nmv[j], g[j], se[j], abs(g[j]) / max(se[j], 1e-12), disp[j]))
  flush.console()
  res[[rep]] <- list(rep = rep, p = p, g = g, se = se, H = H, disp = disp)
}
saveRDS(res, sprintf("out/llcheck_%s_N%d.rds", STEP, N))
cat(sprintf("\n== %s, N=%d: mean |t| and mean |shift| over %d reps ==\n", STEP, N, length(res)))
for (j in seq_along(nmv))
  cat(sprintf("   %-9s mean|t|=%7.2f  mean|shift|=%.4f  signs=%s\n", nmv[j],
      mean(vapply(res, function(x) abs(x$g[j]) / max(x$se[j], 1e-12), numeric(1))),
      mean(vapply(res, function(x) abs(x$disp[j]), numeric(1))),
      paste(vapply(res, function(x) if (x$disp[j] > 0) "+" else "-", character(1)), collapse = "")))
