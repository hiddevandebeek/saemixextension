## E2/E3 with a null arm, so that dLL is interpretable.
##
## The family-selected vine has at least as many parameters as the all-Gaussian
## one, so 2*(llCop - llGauss) > 0 is EXPECTED even under a Gaussian truth.
## Both arms are therefore run and compared: `gauss` calibrates the statistic,
## `alt` measures whether a genuine non-Gaussian dependence is recovered.
##
## E3: for a 2-cmt IV bolus Cmax = Dose/V1 and AUC = Dose/CL, so the corner that
## matters clinically is joint extremes of (V1, CL).  Gumbel puts tail dependence
## in the upper corner and is asymptotically independent in the lower, so both
## are reported and whichever moves is the answer.
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)})
setwd("C:/package/saemix-copula/copula")
source("R/simData.R"); source("R/marginalLL.R"); source("R/diagnostic.R")

a <- commandArgs(trailingOnly = TRUE)
NREP <- as.integer(a[1]); N <- as.integer(a[2]); FAM <- a[3]
MODE <- if (length(a) > 3) a[4] else "sa"
TAG  <- if (length(a) > 4) a[5] else sprintf("rec2_N%d_%s_%s", N, FAM, MODE)

## t and joe dropped from the search to keep the M-step affordable; gaussian is
## kept because it is the nested null and gumbel because it is the truth.
FAMSET <- c("gaussian", "clayton", "gumbel", "frank")
tw <- makeTwinVines(PK_R, FAM)
ctl <- function(sd) list(seed = sd, save = FALSE, save.graphs = FALSE, print = FALSE,
                         displayProgress = FALSE, nbiter.saemix = c(150, 100),
                         nbiter.mcmc = c(2, 2, 2, 0), warnings = FALSE)
corners <- function(vine, q = 0.10, nsim = 1e5) {
  u <- withSeed(31, rvinecopulib::rvinecop(nsim, vine)); d <- ncol(u)
  c(allLow = mean(rowSums(u < q) == d), allHigh = mean(rowSums(u > 1 - q) == d),
    pairLow = mean(u[, 1] < q & u[, 2] < q) / q,
    pairHigh = mean(u[, 1] > 1 - q & u[, 2] > 1 - q) / q)
}
out <- list()
for (arm in c("gauss", "alt")) for (rep in seq_len(NREP)) {
  set.seed(if (arm == "gauss") 60000 + rep else 50000 + rep)
  s <- simPK(N, tw[[arm]])
  dat <- pkSaemixData(s$data)
  copulaClear()
  fG <- try(saemix::saemix(pkSaemixModel(), dat, ctl(rep)), silent = TRUE)
  if (inherits(fG, "try-error")) next
  muG <- log(fG@results@fixed.effects); OmG <- fG@results@omega
  peG <- fG@results@respar[2]; sdG <- sqrt(diag(OmG)); vnGhat <- etaVineGaussian(cov2cor(OmG))
  copulaSet(etaVineGaussian(cov2cor(OmG)), sdG, familySet = FAMSET,
            mode = MODE, refitEvery = 5L, fitFrom = 30L)
  fC <- try(saemix::saemix(pkSaemixModel(), dat, ctl(rep)), silent = TRUE)
  if (inherits(fC, "try-error")) { copulaClear(); next }
  st <- copulaGet(); vnC <- st$vine; sdC <- st$sd
  muC <- log(fC@results@fixed.effects); peC <- fC@results@respar[2]
  copulaClear()
  llG <- marginalLL(s$data, muG, priorMVN(OmG), peG, M = 2000, seed = 900 + rep)
  llC <- marginalLL(s$data, muC, priorVine(vnC, sdC), peC, M = 2000, seed = 900 + rep)
  fl <- unlist(vnC$pair_copulas, recursive = FALSE)
  cT <- corners(tw[[arm]]); cG <- corners(vnGhat); cC <- corners(vnC)
  out[[length(out) + 1]] <- data.frame(rep = rep, arm = arm, N = N, mode = MODE,
    dLL = 2 * (llC$ll - llG$ll), llG = llG$ll, llC = llC$ll,
    essG = llG$essMin, essC = llC$essMin,
    famCop = paste(vapply(fl, function(b) b$family, character(1)), collapse = ","),
    tauCop = paste(round(vapply(fl, function(b)
      if (b$family == "indep") 0 else rvinecopulib::par_to_ktau(b), numeric(1)), 3), collapse = ","),
    tauTrue = paste(round(if (arm == "gauss") tw$tau else tw$tauAlt, 3), collapse = ","),
    muG = paste(round(exp(muG), 2), collapse = ","), muC = paste(round(exp(muC), 2), collapse = ","),
    sdG = paste(round(sdG, 3), collapse = ","), sdC = paste(round(sdC, 3), collapse = ","),
    allHighTrue = cT[["allHigh"]], allHighGauss = cG[["allHigh"]], allHighCop = cC[["allHigh"]],
    allLowTrue = cT[["allLow"]], allLowGauss = cG[["allLow"]], allLowCop = cC[["allLow"]],
    pairHighTrue = cT[["pairHigh"]], pairHighGauss = cG[["pairHigh"]], pairHighCop = cC[["pairHigh"]],
    stringsAsFactors = FALSE)
  cat(sprintf("[%s] %-5s rep %d/%d dLL=%7.2f ess=%.0f/%.0f pairHigh t/g/c=%.2f/%.2f/%.2f %s\n",
      format(Sys.time(), "%H:%M:%S"), arm, rep, NREP, 2 * (llC$ll - llG$ll),
      llG$essMin, llC$essMin, cT[["pairHigh"]], cG[["pairHigh"]], cC[["pairHigh"]],
      out[[length(out)]]$famCop)); flush.console()
}
res <- do.call(rbind, out); saveRDS(res, file.path("out", paste0(TAG, ".rds")))
cat("saved out/", TAG, ".rds rows=", nrow(res), "\n", sep = "")
