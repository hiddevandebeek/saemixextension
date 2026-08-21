## Does the tree-1 STRUCTURE matter for a truncated vine?
##
## We hard-wired dvine_structure(1:d).  For a FULL vine that is nearly a free
## choice (all vines on d vars span the same model class); for a TRUNCATED vine
## the tree-1 structure IS the model -- it says which conditional independences
## you impose.  At d=4 there are only d^(d-2) = 16 one-truncated vines (Cayley),
## so this is enumerable, not searchable.
##
## Arms (all 1-truncated, Gumbel truth, same data/seed):
##   fixed  dvine_structure(1:4)  -- what we have been using
##   best   argmax loglik over all 16, chosen ONCE on post-burn-in draws then frozen
##   worst  argmin loglik over all 16 -- the cost of getting it wrong
## Nikoloulopoulos reports the decomposition does not matter at small N; check.
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R"); source("R/simEta.R"); source("R/simpleModels.R")
NREP <- 5; N <- 150
m <- MODELS$iv2; d <- 4; sdv <- rep(ETA_SD, d); truth <- c(m$true, sdv)
tw <- makeTwinVines(etaR(d), "gumbel")
FAMSET <- c("gaussian", "clayton", "gumbel", "frank")
ctl <- function(sd) list(seed = sd, save = FALSE, save.graphs = FALSE, print = FALSE,
                         displayProgress = FALSE, nbiter.saemix = c(150, 100),
                         nbiter.mcmc = c(2, 2, 2, 0), warnings = FALSE)
## all 16 labelled trees on 4 nodes = 12 D-vine paths + 4 C-vine stars
allStructures <- function() {
  out <- list()
  for (p in unique(lapply(combinat::permn(1:4), identity))) {
    s <- try(rvinecopulib::dvine_structure(unlist(p), trunc_lvl = 1), silent = TRUE)
    if (!inherits(s, "try-error")) out[[length(out) + 1]] <- s
  }
  for (r in 1:4) {
    s <- try(rvinecopulib::cvine_structure(c(r, setdiff(1:4, r)), trunc_lvl = 1), silent = TRUE)
    if (!inherits(s, "try-error")) out[[length(out) + 1]] <- s
  }
  out
}
STR <- allStructures()
cat("candidate 1-truncated structures:", length(STR), "\n\n")
corners <- function(vine, q = 0.10, nsim = 1e5) {
  u <- withSeed(31, rvinecopulib::rvinecop(nsim, vine))
  c(pairHigh = mean(u[, 1] > 1 - q & u[, 2] > 1 - q) / q)
}
cTrue <- corners(tw$alt)
res <- list()
for (rr in seq_len(NREP)) {
  set.seed(8000 + rr)
  s <- simModel(m, N, tw$alt); dat <- saemixDataFor(s$data); mod <- saemixModelFor(m)
  ## stock fit -> conditional draws -> score all 16 ONCE, then freeze
  copulaClear()
  fS <- try(saemix::saemix(mod, dat, ctl(rr)), silent = TRUE); if (inherits(fS, "try-error")) next
  sdS <- sqrt(diag(fS@results@omega))
  cd <- saemix::conddist.saemix(fS, nsamp = 3, plot = FALSE)
  idx <- which(diag(fS@results@omega) > 1e-10)
  E <- cd@results@phi.samp[, idx, 1] - cd@results@mean.phi[, idx]
  uu <- pmin(pmax(pnorm(sweep(E, 2, sdS, "/")), 1e-6), 1 - 1e-6)
  ll <- vapply(STR, function(st) {
    f <- try(rvinecopulib::vinecop(uu, structure = st, family_set = FAMSET), silent = TRUE)
    if (inherits(f, "try-error")) -Inf else f$loglik }, numeric(1))
  pick <- list(fixed = 1L, best = which.max(ll), worst = which.min(ll))
  for (a in names(pick)) {
    st <- if (a == "fixed") rvinecopulib::dvine_structure(1:4, trunc_lvl = 1) else STR[[pick[[a]]]]
    v0 <- rvinecopulib::vinecop_dist(list(lapply(1:3, function(e)
            rvinecopulib::bicop_dist("gaussian", parameters = 0.3))), st)
    copulaSet(v0, sdS, familySet = FAMSET, mode = "sa", refitEvery = 3L,
              fitFrom = 30L, truncLvl = 1L)
    fC <- try(saemix::saemix(mod, dat, ctl(rr)), silent = TRUE)
    if (inherits(fC, "try-error")) { copulaClear(); next }
    stt <- copulaGet(); copulaClear()
    res[[length(res) + 1]] <- list(rep = rr, arm = a,
      parErr = mean(abs(c(fC@results@fixed.effects, stt$sd) - truth) / truth),
      pairHigh = corners(stt$vine)[["pairHigh"]],
      order = paste(rvinecopulib::get_structure(st)$order, collapse = "-"))
  }
  cat(sprintf("rep %d  loglik spread over 16 structures: %.2f .. %.2f (range %.2f)\n",
      rr, min(ll), max(ll), max(ll) - min(ll))); flush.console()
}
saveRDS(res, "out/structure.rds")
cat(sprintf("\nTRUTH pairHigh = %.4f\n", cTrue))
cat(sprintf("%-7s %9s %12s %s\n", "arm", "parErr", "pairHigh", "orders picked"))
for (a in c("fixed", "best", "worst")) {
  k <- vapply(res, function(x) x$arm == a, logical(1)); if (!any(k)) next
  cat(sprintf("%-7s %9.4f %12.4f (err %+.4f)  %s\n", a,
      mean(vapply(res[k], function(x) x$parErr, numeric(1))),
      mean(vapply(res[k], function(x) x$pairHigh, numeric(1))),
      mean(vapply(res[k], function(x) x$pairHigh, numeric(1))) - cTrue,
      paste(unique(vapply(res[k], function(x) x$order, character(1))), collapse = " ")))
}
