## Which component drives the copula path?  Freeze each in turn.  Also vary the
## burn-in length: SAEM's stability comes from the decaying-gamma phase, and a
## long gamma=1 burn-in gives the copula M-step no averaging at all.
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)})
setwd("C:/package/saemix-copula/copula"); source("R/simData.R")

set.seed(77); vnG <- etaVineGaussian(PK_R); s <- simPK(120, vnG)
ctl <- function(b, k) list(seed = 1, save = FALSE, save.graphs = FALSE, print = FALSE,
                           displayProgress = FALSE, nbiter.saemix = c(b, k),
                           nbiter.mcmc = c(2, 2, 2, 0), warnings = FALSE)
line <- function(nm, f, sdv) cat(sprintf("%-22s %-30s %-26s %.4f\n", nm,
  paste(round(f@results@fixed.effects, 3), collapse = ","),
  paste(round(sdv, 3), collapse = ","), f@results@respar[2]))

cat(sprintf("%-22s %-30s %-26s %s\n", "config", "fixed (true 10,3,2,30)", "sd (true .3,.3,.4,.4)", "pres(.1)"))
for (bk in list(c(150, 80), c(40, 190))) {
  copulaClear()
  f <- saemix::saemix(pkSaemixModel(), pkSaemixData(s$data), ctl(bk[1], bk[2]))
  line(sprintf("stock burn=%d", bk[1]), f, sqrt(diag(f@results@omega)))
  for (cfg in list(c(TRUE, TRUE), c(FALSE, TRUE), c(TRUE, FALSE), c(FALSE, FALSE))) {
    copulaSet(etaVineGaussian(PK_R), PK_SD, familySet = "gaussian",
              freezeSd = cfg[1], freezeVine = cfg[2])
    f <- try(saemix::saemix(pkSaemixModel(), pkSaemixData(s$data), ctl(bk[1], bk[2])), silent = TRUE)
    nm <- sprintf("sd%s vine%s burn=%d", if (cfg[1]) "FRZ" else "fit",
                  if (cfg[2]) "FRZ" else "fit", bk[1])
    if (inherits(f, "try-error")) cat(sprintf("%-22s ERROR\n", nm))
    else line(nm, f, copulaGet()$sd)
  }
}
copulaClear()
