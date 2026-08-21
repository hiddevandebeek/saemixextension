## E2/E3: does route A recover a non-Gaussian eta dependence, and does it change
## anything anyone cares about?
##
## E2 recovery : fitted taus / families vs truth, and the EXACT marginal
##               log-likelihood (independent IS evaluator) of the Gaussian fit
##               vs the copula fit on the same data.  The Gaussian model is
##               nested, so a positive difference is a real LRT.
## E3 relevance: joint tail probabilities.  For a 2-cmt IV bolus, Cmax = Dose/V1
##               and AUC = Dose/CL, so the clinically dangerous corner is LOW V1
##               AND LOW CL simultaneously.  Gumbel puts tail dependence in the
##               upper corner (asymptotically independent in the lower), survival
##               Gumbel (rotation 180) puts it in the dangerous one.  If the
##               Gaussian fit misses that probability, the copula earns its keep;
##               if not, it does not.
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)})
setwd("C:/package/saemix-copula/copula")
source("R/simData.R"); source("R/marginalLL.R"); source("R/diagnostic.R")

args <- commandArgs(trailingOnly = TRUE)
NREP <- if (length(args) > 0) as.integer(args[1]) else 8
N    <- if (length(args) > 1) as.integer(args[2]) else 200
FAM  <- if (length(args) > 2) args[3] else "gumbel"
TAG  <- if (length(args) > 3) args[4] else sprintf("rec_N%d_%s", N, FAM)

tw <- makeTwinVines(PK_R, FAM)
ctl <- function(sd) list(seed = sd, save = FALSE, save.graphs = FALSE, print = FALSE,
                         displayProgress = FALSE, nbiter.saemix = c(150, 120),
                         nbiter.mcmc = c(2, 2, 2, 0), warnings = FALSE)

## P(all d etas below their q-quantile) -- the joint "everything low" corner.
## Both corners: Gumbel puts tail dependence in the UPPER one and is
## asymptotically independent in the lower, so report whichever moves.
corners <- function(vine, q = 0.10, nsim = 2e5) {
  u <- withSeed(31, rvinecopulib::rvinecop(nsim, vine))
  d <- ncol(u)
  c(allLow  = mean(rowSums(u < q) == d),
    allHigh = mean(rowSums(u > 1 - q) == d),
    pairLow = mean(u[, 1] < q & u[, 2] < q) / q,
    pairHigh= mean(u[, 1] > 1 - q & u[, 2] > 1 - q) / q)
}

out <- list()
for (rep in seq_len(NREP)) {
  set.seed(50000 + rep)
  s <- simPK(N, tw$alt)
  ## --- Gaussian arm (stock saemix) ---
  copulaClear()
  fG <- try(saemix::saemix(pkSaemixModel(), pkSaemixData(s$data), ctl(rep)), silent = TRUE)
  if (inherits(fG, "try-error")) { cat("rep", rep, "stock failed\n"); next }
  muG <- log(fG@results@fixed.effects); OmG <- fG@results@omega; peG <- fG@results@respar[2]
  sdG <- sqrt(diag(OmG)); vnGhat <- etaVineGaussian(cov2cor(OmG))

  ## --- Copula arm (route A, family selection) ---
  copulaSet(etaVineGaussian(cov2cor(OmG)), sdG,
            familySet = c("gaussian", "t", "clayton", "gumbel", "frank", "joe"),
            refitEvery = 2L, fitFrom = 20L)
  fC <- try(saemix::saemix(pkSaemixModel(), pkSaemixData(s$data), ctl(rep)), silent = TRUE)
  if (inherits(fC, "try-error")) { cat("rep", rep, "copula failed\n"); copulaClear(); next }
  st <- copulaGet(); vnC <- st$vine; sdC <- st$sd
  muC <- log(fC@results@fixed.effects); peC <- fC@results@respar[2]
  copulaClear()

  ## --- exact marginal likelihood of each fitted model, same data ---
  llG <- marginalLL(s$data, muG, priorMVN(OmG), peG, M = 3000, seed = 900 + rep)
  llC <- marginalLL(s$data, muC, priorVine(vnC, sdC), peC, M = 3000, seed = 900 + rep)

  famC <- paste(vapply(unlist(vnC$pair_copulas, recursive = FALSE),
                       function(b) b$family, character(1)), collapse = ",")
  tauC <- vapply(unlist(vnC$pair_copulas, recursive = FALSE),
                 function(b) if (b$family == "indep") 0 else rvinecopulib::par_to_ktau(b), numeric(1))

  out[[length(out) + 1]] <- data.frame(
    rep = rep, N = N, fam = FAM,
    llGauss = llG$ll, llCop = llC$ll, dLL = 2 * (llC$ll - llG$ll),
    essGauss = llG$essMin, essCop = llC$essMin,
    muG = paste(round(exp(muG), 3), collapse = ","), muC = paste(round(exp(muC), 3), collapse = ","),
    sdG = paste(round(sdG, 3), collapse = ","), sdC = paste(round(sdC, 3), collapse = ","),
    tauTrue = paste(round(tw$tauAlt, 3), collapse = ","),
    tauCop = paste(round(tauC, 3), collapse = ","),
    famCop = famC,
    stringsAsFactors = FALSE)
  cT <- corners(tw$alt); cG <- corners(vnGhat); cC <- corners(vnC)
  for (nm in names(cT)) {
    out[[length(out)]][[paste0(nm, "True")]]  <- cT[[nm]]
    out[[length(out)]][[paste0(nm, "Gauss")]] <- cG[[nm]]
    out[[length(out)]][[paste0(nm, "Cop")]]   <- cC[[nm]]
  }
  cat(sprintf("[%s] rep %d/%d  dLL=%.2f  ess=%.0f/%.0f  low: true=%.4f gauss=%.4f cop=%.4f\n",
              format(Sys.time(), "%H:%M:%S"), rep, NREP, out[[length(out)]]$dLL,
              llG$essMin, llC$essMin, out[[length(out)]]$lowTrue,
              out[[length(out)]]$lowGauss, out[[length(out)]]$lowCop))
  flush.console()
}
res <- do.call(rbind, out)
saveRDS(res, file.path("out", paste0(TAG, ".rds")))
cat("saved out/", TAG, ".rds  rows=", nrow(res), "\n", sep = "")
