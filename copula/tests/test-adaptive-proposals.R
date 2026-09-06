## Fast tests for the end-of-fit adaptive proposal phase.

suppressMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R")
source("R/simEta.R")
source("R/adaptiveVine.R")

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-60s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}

set.seed(52)
d <- 4L; N <- 80L; ns <- 3L; sd <- rep(.30, d)
R <- .55^abs(outer(seq_len(d), seq_len(d), "-"))
truth <- makeTwinVines(R, "gumbel")$alt
draws <- array(0, c(N, d, ns))
for (s in seq_len(ns)) draws[, , s] <- rEtaVine(N, truth, sd)

ps <- copulaAdaptiveProposalsFromSlices(
  draws, sd, truncLevels = c(1L, Inf), top = 3L,
  familySet = c("gaussian", "gumbel"))
ok("A1 adaptive phase emits a proposal set", inherits(ps, "saemixCopulaProposalSet"))
ok("A2 structure/family/truncation may vary across proposals",
   nrow(ps$proposals) == ns * 2L && length(unique(ps$proposals$truncLvl)) == 2L)
ok("A3 shortlist is bounded and every candidate has a fingerprint",
   length(ps$candidates) <= 3L &&
   all(vapply(ps$candidates, function(x) nzchar(x$key), logical(1))))
ok("A4 proposal ranking uses subject-level N in BIC",
   all(abs(ps$proposals$proposalBIC -
     (-2 * ps$proposals$proposalLogLik + log(N) * ps$proposals$npar)) < 1e-10))
ok("A5 output requires fresh fixed-model likelihood comparison",
   grepl("fresh start", ps$warning) && grepl("Gaussian baseline", ps$warning))

## A completed copula run can drive another adaptive cycle from its weighted
## subject-identified pool rather than from EBEs.
copulaSet(truth, sd, mode = "joint", guard = FALSE,
  populationAlgorithm = "common-q")
X <- matrix(1, N * ns, 1)
locMap <- matrix(0, 1, d)
Epool <- do.call(rbind, lapply(seq_len(ns), function(s) draws[, , s]))
copulaPoolUpdate(Epool, gamma = 1, nchains = ns, phiM = Epool, xM = X,
  locMap = locMap, beta = 0, betaFree = integer(0),
  subject = rep(seq_len(N), ns), iteration = 100L)
st <- copulaGet(); st$betaFitted <- 0
sl <- copulaAdaptiveSlicesFromState(st, nsamp = 2L, seed = 9L)
ok("A6 later adaptive cycles sample one particle per subject",
   identical(dim(sl), c(N, d, 2L)) && !anyNA(sl))
ps2 <- copulaAdaptiveProposalsFromState(st, nsamp = 2L, seed = 9L,
  truncLevels = 1L, familySet = c("gaussian", "gumbel"), top = 2L)
ok("A7 copula-state adaptation emits fresh candidate proposals",
   inherits(ps2, "saemixCopulaProposalSet") && length(ps2$candidates) >= 1L)
qs <- copulaAdaptiveQStep(st, ps2, criterion = "bic", includeIncumbent = TRUE,
                          maxit = 40L)
ok("A8 structural jumps are compared on one common Q measure",
   inherits(qs, "saemixCopulaQStep") && isTRUE(qs$commonMeasure) &&
     nrow(qs$ranking) >= 1L)
ok("A9 every structural candidate receives a full native parameter count",
   all(qs$ranking$npar >= d) && all(is.finite(qs$ranking$Q)) &&
     all(is.finite(qs$ranking$QBIC)))
marginAlt <- c(list(copulaMarginStudent(.30, 8)),
               lapply(rep(.30, d - 1L), copulaMarginNormal))
qm <- copulaAdaptiveMarginQStep(st, list(marginAlt), criterion = "bic",
                                includeIncumbent = TRUE, maxit = 10L)
ok("A10 margin-family jumps use the same common Q machinery",
   inherits(qm, "saemixCopulaQStep") && isTRUE(qm$commonMeasure) &&
     setequal(qm$ranking$nMarginParameters, c(d, d + 1L)))
md <- copulaMargin("binary-centred", c(prob=.5), 1e-5, 1-1e-5,
  log_density=function(x,par) dbinom(x+.5,1,par[1],log=TRUE),
  cdf=function(x,par) pbinom(x+.5,1,par[1]),
  cdf_left=function(x,par) pbinom(x-.5,1,par[1]),
  quantile=function(u,par) qbinom(u,1,par[1])-.5,
  type="discrete", centered=TRUE,
  scale=function(par) sqrt(par[1]*(1-par[1])))
typeJump <- try(copulaAdaptiveMarginQStep(st,
  list(c(list(md), lapply(rep(.30,d-1L),copulaMarginNormal)))), silent=TRUE)
ok("A11 continuous/discrete adaptive jumps are refused for this augmentation",
   inherits(typeJump,"try-error"))
copulaClear()

cat(sprintf("\n%d failure(s)\n", nFail))
if (nFail > 0) quit(status = 1)
