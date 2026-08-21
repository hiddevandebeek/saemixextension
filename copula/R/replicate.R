## One replicate of the route-C experiment.
##
## Three levels of the SAME statistic, which together decompose exactly how much
## detectability the Gaussian-prior E-step destroys:
##   oracle : true etas                       -- upper bound, unattainable
##   draws  : conditional draws p(eta|y)      -- what route C can actually use
##   ebe    : conditional means (EBEs)        -- what people actually do

source("R/simData.R"); source("R/diagnostic.R")

runReplicate <- function(seed, arm, tw, N = 200, nsamp = 5, propErr = 0.10,
                         nbiter = c(200, 100), how = "param") {
  set.seed(seed)
  vine <- tw[[arm]]
  s <- simPK(N, vine, propErr = propErr)
  fit <- saemix::saemix(pkSaemixModel(), pkSaemixData(s$data),
                        list(seed = seed, save = FALSE, save.graphs = FALSE,
                             print = FALSE, displayProgress = FALSE,
                             nbiter.saemix = nbiter, warnings = FALSE))
  cd <- condEtaDraws(fit, nsamp = nsamp)

  ## shrinkage per eta, the standard eta-shrinkage definition
  shrink <- 1 - apply(cd$ebe, 2, stats::sd) / cd$sd

  out <- do.call(rbind, lapply(c("param", "rank"), function(hh) {
    oracle <- vineLRT(etaToUnif(s$eta, PK_SD, hh))
    draws  <- vineLRTdraws(cd$draws, cd$sd, hh)
    ebe    <- vineLRT(etaToUnif(cd$ebe, cd$sd, hh))
    data.frame(seed = seed, arm = arm, N = N, how = hh,
               statOracle = oracle$stat, statDraws = draws$stat, statEbe = ebe$stat,
               ngOracle = oracle$nNonGauss, ngDraws = draws$nNonGauss, ngEbe = ebe$nNonGauss,
               famOracle = oracle$families, famDraws = draws$families, famEbe = ebe$families,
               stringsAsFactors = FALSE)
  }))
  out$shrinkMax <- max(shrink); out$shrinkMean <- mean(shrink)
  out$sdFit <- paste(round(cd$sd, 3), collapse = ",")
  out
}
