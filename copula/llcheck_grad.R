## Stationarity test: is the IFM fixed point a stationary point of the EXACT
## marginal likelihood?  Gradient with honest Monte Carlo error bars, plus the
## implied parameter displacement g/|H_jj| so the answer is interpretable rather
## than just "nonzero".
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R"); source("R/simEta.R"); source("R/simpleModels.R"); source("llcheck_core.R")

a <- commandArgs(trailingOnly = TRUE)
NREP <- if (length(a) > 0) as.integer(a[1]) else 4
N    <- if (length(a) > 1) as.integer(a[2]) else 100
M    <- if (length(a) > 2) as.integer(a[3]) else 1500
m <- MODELS$iv1; d <- 2; PE <- 0.10
tw <- makeTwinVines(etaR(d), "gumbel")
FAM <- etaVineFlat(tw$alt)[[1]]$family; ROT <- etaVineFlat(tw$alt)[[1]]$rotation
TAU_TRUE <- tw$tauAlt[1]
nmv <- c("mu1", "mu2", "logSd1", "logSd2", "atanhTau")

llAt <- function(p, dat, seed) {
  q <- unpack(p)
  isPass(dat, q$mu, q$sdv, vine2FromTau(q$tau, FAM, ROT), PE, m$f,
         M = M, seed = seed, keepDraws = FALSE)$ll
}
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

res <- list()
for (rep in seq_len(NREP)) {
  set.seed(100 + rep); dat <- simModel(m, N, tw$alt)$data
  p <- pack(log(m$true), rep(ETA_SD, d), TAU_TRUE)
  for (it in 1:20) { pn <- ifmStep(p, dat, 7); dl <- max(abs(pn - p)); p <- pn; if (dl < 1e-4) break }
  pIFM <- p
  seeds <- c(31, 32, 33, 34, 35); h <- 0.02
  G <- sapply(seeds, function(sd_) sapply(seq_along(pIFM), function(j) {
    pu <- pIFM; pl <- pIFM; pu[j] <- pu[j] + h; pl[j] <- pl[j] - h
    (llAt(pu, dat, sd_) - llAt(pl, dat, sd_)) / (2 * h) }))
  H <- sapply(seq_along(pIFM), function(j) {
    pu <- pIFM; pl <- pIFM; pu[j] <- pu[j] + h; pl[j] <- pl[j] - h
    (llAt(pu, dat, 31) - 2 * llAt(pIFM, dat, 31) + llAt(pl, dat, 31)) / h^2 })
  g <- rowMeans(G); se <- apply(G, 1, sd) / sqrt(length(seeds))
  disp <- ifelse(H < 0, -g / H, NA)
  cat(sprintf("[rep %d] IFM p = %s  (iters %d)\n", rep, paste(sprintf("%.4f", pIFM), collapse = ","), it))
  for (j in seq_along(pIFM))
    cat(sprintf("   %-9s grad=%+9.3f  se=%7.3f  |t|=%5.2f  Hjj=%+10.1f  implied shift=%+.4f\n",
        nmv[j], g[j], se[j], abs(g[j]) / max(se[j], 1e-12), H[j], disp[j]))
  flush.console()
  res[[rep]] <- list(rep = rep, pIFM = pIFM, g = g, se = se, H = H, disp = disp, G = G)
}
saveRDS(res, "out/llcheck_grad.rds")
cat("saved out/llcheck_grad.rds\n")
