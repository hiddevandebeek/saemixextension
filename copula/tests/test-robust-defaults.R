## Fast contract tests for defaults that make the copula path safe by default.

suppressMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R")

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-60s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}

R4 <- .55^abs(outer(1:4, 1:4, "-"))
v4 <- etaVineGaussian(R4)

copulaSet(v4, rep(.30, 4))
s <- copulaGet()
fam <- vapply(copulaPadFlat(s$vine, 4), function(b) b$family, character(1))
ok("R1 default estimator is the joint common-Q path", identical(s$mode, "joint"), s$mode)
ok("R2 all-Gaussian default preserves the full MVN nested null", is.infinite(s$truncLvl),
   sprintf("truncLvl=%s", s$truncLvl))
ok("R2b full Gaussian conditional correlations are retained",
   all(fam == "gaussian"), paste(fam, collapse = ","))
ok("R2c canonical Gaussian vine correlation is reconstructed analytically",
   max(abs(copulaImpliedCor() - R4)) < 1e-12,
   sprintf("maxabs=%.2e", max(abs(copulaImpliedCor() - R4))))
ok("R3 omitted familySet freezes the supplied edge families",
   !is.null(s$famFixed) && identical(fam, s$famFixed))
ok("R3a fixed model stores a discrete fingerprint", isTRUE(s$modelFrozen) &&
   !is.null(s$fingerprint))

v4ng <- etaVine(rep(list(rvinecopulib::bicop_dist("gumbel", parameters = 2)), 6), 4)
copulaSet(v4ng, rep(.30, 4))
s <- copulaGet()
fam <- vapply(copulaPadFlat(s$vine, 4), function(b) b$family, character(1))
ok("R3b supplied non-Gaussian fixed model is not silently truncated",
   is.infinite(s$truncLvl) && all(fam == "gumbel"))

copulaSet(v4, rep(.30, 4), familySet = c("gaussian", "gumbel"))
s <- copulaGet()
ok("R3c family search does not silently change truncation", is.infinite(s$truncLvl))

copulaSet(v4ng, rep(.30, 4), truncLvl = 1L)
s <- copulaGet()
fam <- vapply(copulaPadFlat(s$vine, 4), function(b) b$family, character(1))
ok("R3d explicit tree-1 model is active before the first E-step",
   identical(s$truncLvl, 1L) && all(fam[4:6] == "indep"),
   paste(fam, collapse = ","))

## A finite pool is explicitly an approximation, but it must still incorporate
## new equal-weight smoothing draws rather than freezing on the earliest ties.
suppressWarnings(copulaSet(v4,rep(.30,4),poolMax=1L))
copulaPoolUpdate(matrix(1,5,4),gamma=1,nchains=1,
  subject=seq_len(5),iteration=1L)
copulaPoolUpdate(matrix(2,5,4),gamma=.5,nchains=1,
  subject=seq_len(5),iteration=2L)
s<-copulaGet()
ok("R3e finite equal-weight pool retains the newest particles",
   all(s$poolEta==2)&&all(s$poolIteration==2L))

## Non-dyadic 1/t gains must not let floating-point weight differences outrank
## recency. The cap is in complete batch identities, even if sizes differ.
suppressWarnings(copulaSet(v4,rep(.30,4),poolMax=3L))
for(t in seq_len(11L))
  copulaPoolUpdate(matrix(t,1L,4L),gamma=1/t,nchains=1L,
    subject=1L,iteration=t)
s<-copulaGet()
ok("R3f finite pool keeps newest batches under non-dyadic gains",
   identical(sort(unique(s$poolIteration)),9:11),
   paste(sort(unique(s$poolIteration)),collapse=","))

suppressWarnings(copulaSet(v4,rep(.30,4),poolMax=2L))
copulaPoolUpdate(matrix(1,2L,4L),gamma=1,nchains=1L,
  subject=1:2,iteration=1L)
copulaPoolUpdate(matrix(2,3L,4L),gamma=.5,nchains=1L,
  subject=1:3,iteration=2L)
copulaPoolUpdate(matrix(3,1L,4L),gamma=1/3,nchains=1L,
  subject=1L,iteration=3L)
s<-copulaGet()
ok("R3g finite pool retains complete batches when sizes differ",
   identical(sort(unique(s$poolIteration)),2:3)&&nrow(s$poolEta)==4L)
mixedMetadata<-try(copulaPoolUpdate(matrix(4,1L,4L),gamma=.25,nchains=1L,
  iteration=4L),silent=TRUE)
ok("R3h particle-pool metadata modes cannot change mid-run",
   inherits(mixedMetadata,"try-error"))
badDesignBundle<-try(copulaPoolUpdate(matrix(4,1L,4L),gamma=.25,nchains=1L,
  phiM=matrix(4,1L,4L),subject=1L,iteration=4L),silent=TRUE)
badMultiplicity<-try(copulaPoolUpdate(matrix(4,2L,4L),gamma=.25,nchains=2L,
  subject=1:2,iteration=4L),silent=TRUE)
ok("R3i incomplete design bundles and subject multiplicities fail closed",
   inherits(badDesignBundle,"try-error")&&inherits(badMultiplicity,"try-error"))

copulaSet(v4, rep(.30, 4), truncLvl = Inf, guard = FALSE)
s <- copulaGet()
ok("R4 explicit full-vine choice is preserved", is.infinite(s$truncLvl))
ok("R4b guard can be explicitly disabled", identical(s$guard, FALSE))

## The default guard diagnoses but never mutates the fixed model.
copulaSet(v4, rep(.30, 4), truncLvl = Inf, guardConsecutive = 3L)
g1 <- copulaGuardUpdate(10L, rep(.12, 4), v4)
g2 <- copulaGuardUpdate(11L, rep(.12, 4), v4)
g3 <- suppressWarnings(copulaGuardUpdate(12L, rep(.12, 4), v4))
s <- copulaGet()
gfam <- vapply(copulaPadFlat(g3$vine, 4), function(b) b$family, character(1))
ok("R5 guard ignores fewer than three consecutive bad updates",
   !g1$triggered && !g2$triggered)
ok("R5b default guard leaves the active fixed model unchanged",
   !g3$triggered && all(gfam == "gaussian") && identical(g3$sd, rep(.12, 4)))
ok("R5c collapse warning is recorded for audit",
   isTRUE(s$collapseDetected) && identical(s$collapseAt, 12L) &&
   !isTRUE(s$fallbackTriggered))

## A mutating fallback exists only as an explicit research/legacy action.
copulaSet(v4ng, rep(.30, 4), guardConsecutive = 1L, guardAction = "fallback")
gf <- suppressWarnings(copulaGuardUpdate(20L, rep(.12, 4), v4ng))
s <- copulaGet()
gfam <- vapply(copulaPadFlat(gf$vine, 4), function(b) b$family, character(1))
ok("R5d explicit fallback changes to a Gaussian vine", gf$triggered &&
   all(gfam %in% c("gaussian", "indep")) && isTRUE(s$fallbackTriggered))

## Direct margins diagnose collapse from their model-implied conditional
## dispersion, not from the unrelated MCMC proposal covariance.
v1<-rvinecopulib::vinecop_dist(list(),rvinecopulib::dvine_structure(1))
mDirect<-parameterMarginGamma(5,"mean")
copulaSet(v1,margins=list(mDirect),populationScale="parameter",proposalScale=.3,
  guardConsecutive=2L)
predictor<-matrix(log(8),20L,1L);weights<-rep(1/20,20)
scale0<-parameterMarginDispersion(list(mDirect),predictor,weights)
scaleCollapsed<-parameterMarginDispersion(list(parameterMarginGamma(30,"mean")),
  predictor,weights)
dg1<-copulaGuardUpdate(30L,scaleCollapsed,v1,initial=scale0)
dg2<-suppressWarnings(copulaGuardUpdate(31L,scaleCollapsed,v1,initial=scale0))
s<-copulaGet()
ok("R5e direct guard uses natural-margin dispersion and records collapse",
   !dg1$triggered&&!dg2$triggered&&isTRUE(s$collapseDetected)&&
     scaleCollapsed/scale0<.5)
directFallback<-try(copulaSet(v1,margins=list(mDirect),
  populationScale="parameter",proposalScale=.3,guardAction="fallback"),silent=TRUE)
ok("R5f direct fits reject model-changing fallback guards",
   inherits(directFallback,"try-error"))
ok("R5g a rejected declaration atomically preserves the prior state",
   identical(copulaGet()$fingerprint,
     copulaFingerprint(v1,1L,list(mDirect))))

## Preserve the established positional prefix while the expanded declaration
## interface is available by name.
copulaSet(v4,rep(.30,4),"gaussian",Inf,2L)
positional<-copulaGet()
copulaSet(v4,rep(.30,4),familySet="gaussian",poolMax=Inf,refitEvery=2L)
named<-copulaGet()
ok("R5h legacy positional and current named copulaSet calls agree",
   identical(positional$familySet,named$familySet)&&
     identical(positional$poolMax,named$poolMax)&&
     identical(positional$refitEvery,named$refitEvery)&&
     identical(positional$fingerprint,named$fingerprint))

## Candidate proposal accepts every finite-dimensional parametric family, but
## remains explicitly provisional until a fixed-model refit.
set.seed(33)
E <- rEtaVine(300, v4ng, rep(.30, 4))
cp <- copulaProposeCandidate(E, rep(.30, 4), truncLvl = 1L)
ok("R6 proposal helper returns an explicitly provisional candidate",
   inherits(cp, "saemixCopulaCandidate") && grepl("Candidate proposal only", cp$warning))
tCandidate <- try(copulaProposeCandidate(E, rep(.30, 4), familySet = "t"), silent = TRUE)
ok("R6b multi-parameter t candidates are supported",
   inherits(tCandidate, "saemixCopulaCandidate"))
badFamily <- try(copulaProposeCandidate(E, rep(.30, 4), familySet = "tll"), silent = TRUE)
ok("R6c nonparametric candidates are refused by fixed MLE path",
   inherits(badFamily, "try-error"))

## A staged run is stock Gaussian before activation and restarts its stochastic
## approximation clock when the copula model becomes active.
copulaSet(v4ng, rep(.30, 4), activeFrom = 100L, restartBurn = 20L)
ok("R7 staged copula is inactive during Gaussian warm-up",
   !copulaActiveAt(99L) && copulaActiveAt(100L))
ok("R7b staged stochastic approximation restarts at constant gain",
   identical(copulaStepSize(100L, .01), 1) &&
     identical(copulaStepSize(119L, .01), 1) &&
     abs(copulaStepSize(120L, .01) - .5) < 1e-15 &&
     abs(copulaStepSize(121L, .01) - 1/3) < 1e-15)
copulaActivateFromGaussian(100L, diag(c(.04, .09, .16, .25)), nsim = 200L)
s <- copulaGet()
ok("R7c activation inherits Gaussian marginal scales and freezes identity",
   max(abs(s$sd - c(.2, .3, .4, .5))) < 1e-14 &&
     identical(s$fingerprint, copulaFingerprint(s$vine, 4L, s$margins)))

copulaClear()
cat(sprintf("\n%d failure(s)\n", nFail))
if (nFail > 0) quit(status = 1)
