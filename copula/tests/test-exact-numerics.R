## Exact-target numerical policy: production never silently clips, floors, or
## discards positive stochastic-approximation mass.
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet=TRUE)})

nFail <- 0L
ok <- function(label, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-72s %s %s\n", label,
    if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}

vine <- rvinecopulib::vinecop_dist(
  list(list(rvinecopulib::bicop_dist("gaussian", parameters = .35))),
  rvinecopulib::dvine_structure(1:2))
margins <- list(copulaMarginNormal(.3), copulaMarginNormal(.4))
x <- matrix(c(-.2, .1, .05, -.15), ncol = 2, byrow = TRUE)
ev <- copulaMarginsEvaluate(x, margins)
manual <- rowSums(cbind(dnorm(x[, 1], 0, .3, log = TRUE),
                        dnorm(x[, 2], 0, .4, log = TRUE))) +
  log(rvinecopulib::dvinecop(ev$vine_u, vine))
exact <- copulaLogPrior(x, vine, margins = margins)
ok("E1 exact interior density equals the direct declared formula",
   max(abs(exact - manual)) < 1e-13)

boundary <- matrix(c(40, 0), nrow = 1)
failed <- try(copulaMarginsEvaluate(boundary, margins), silent = TRUE)
regularized <- copulaMarginsEvaluate(boundary, margins, eps = 1e-6,
  numericalPolicy = "clip")
ok("E2 an exact boundary PIT fails instead of changing the target",
   inherits(failed, "try-error"))
ok("E3 clipping is explicit and reports every changed coordinate",
   regularized$n_clipped == 1L && regularized$u[1, 1] == 1 - 1e-6)

pm <- list(parameterMarginNormal(.3), parameterMarginNormal(.4))
psi <- matrix(c(40, 0), nrow = 1); predictor <- matrix(0, 1, 2)
pfailed <- try(parameterMarginsEvaluate(psi, predictor, pm), silent = TRUE)
preg <- parameterMarginsEvaluate(psi, predictor, pm, eps = 1e-6,
  numericalPolicy = "clip")
ok("E4 direct-parameter exact PITs obey the same fail-fast contract",
   inherits(pfailed, "try-error") && preg$n_clipped == 1L)

copulaSet(vine, margins = margins, poolMax = Inf, poolPruneRel = 0,
  numericalPolicy = "exact", guard = FALSE,
  populationAlgorithm = "common-q")
for (k in seq_len(8L))
  copulaPoolUpdate(matrix(c(.01 * k, -.02 * k), nrow = 1),
    gamma = if (k <= 3L) 1 else 1 / (k - 2), nchains = 1L, iteration = k)
state <- copulaGet(); copulaClear()
ok("E5 proof-oriented pool retains all positive-weight particles",
   state$poolPrunedCount == 0L && state$poolPrunedMass == 0 &&
     nrow(state$poolEta) == sum(state$poolW > 0))

uniform<-copulaMarginDistribution("unif",c(min=-1,max=1),c(-1,-1),c(1,1),
  free=c(FALSE,FALSE),centered=TRUE,scale=function(par)1/sqrt(3),
  scale_is_sd=TRUE)
outside<-copulaLogPrior(matrix(c(2,0),nrow=1),vine,
  margins=list(uniform,uniform))
ok("E6 zero density outside a valid support is returned as minus infinity",
   identical(as.numeric(outside),-Inf))

badRegularized<-try(copulaSet(vine,margins=margins,
  numericalPolicy="regularized"),silent=TRUE)
badLegacy<-try(copulaSet(vine,margins=margins,mode="pool"),silent=TRUE)
copulaClear()
ok("E7 fitted states are restricted to the exact common-Q path",
   inherits(badRegularized,"try-error")&&inherits(badLegacy,"try-error"))

gLayout<-copulaParameterLayout(list(
  rvinecopulib::bicop_dist("gumbel",parameters=1)))
jLayout<-copulaParameterLayout(list(
  rvinecopulib::bicop_dist("joe",parameters=1)))
ok("E8 valid independence boundaries are not silently moved inward",
   gLayout$lower[1]==1&&gLayout$par[1]==1&&
     jLayout$lower[1]==1&&jLayout$par[1]==1)

v3<-etaVineFromFlat(rep(list(
  rvinecopulib::bicop_dist("gumbel",parameters=1.2)),3L),3L,
  rvinecopulib::dvine_structure(1:3))
extreme<-matrix(c(40,0,0),nrow=1)
cppFailure<-try(copulaMaximiseJointDirectCpp(extreme,1,c(1,1,1),v3,3L,
  maxit=2L,withMu=FALSE),silent=TRUE)
ok("E9 generic compiled optimization does not fit a saturated Gaussian start",
   inherits(cppFailure,"try-error")||identical(cppFailure$conv,-1L))

singularLayout<-copulaParameterLayout(list(
  rvinecopulib::bicop_dist("gaussian",parameters=1)))
ok("E10 singular Gaussian correlation endpoints are kept out of density fits",
   singularLayout$upper[1]<1&&singularLayout$par[1]<1)

copulaSet(vine,margins=margins,poolMax=Inf,poolPruneRel=0,
  numericalPolicy="exact",guard=FALSE,populationAlgorithm="common-q")
copulaPoolUpdate(matrix(c(0,0),nrow=1),gamma=0,nchains=1L,iteration=0L)
zeroState<-copulaGet()
badChains<-try(copulaPoolUpdate(matrix(c(0,0),nrow=1),gamma=1,
  nchains=Inf,iteration=1L),silent=TRUE)
copulaClear()
ok("E11 zero gain is a no-op and nonfinite chain counts fail cleanly",
   is.null(zeroState$poolEta)&&inherits(badChains,"try-error"))

## Both the Markov energy and the empirical-Q helper must use the same native
## multivariate-Normal specialization. The earlier implementation covered only
## copulaUeta(), so a long Gaussian-null fit could fail later in copulaLogPrior.
sd4<-c(.30,.35,.40,.45)
v4<-etaVineFromFlat(rep(list(
  rvinecopulib::bicop_dist("gaussian",parameters=.2)),6L),4L,
  rvinecopulib::dvine_structure(1:4))
m4<-lapply(sd4,copulaMarginNormal)
x4<-matrix(c(50,-50,40,-40),nrow=1)
lp4<-copulaLogPrior(x4,v4,margins=m4)
R4<-copulaGaussianDvineCor(v4,4L);O4<-diag(sd4)%*%R4%*%diag(sd4)
oracle4<--.5*(4*log(2*pi)+as.numeric(determinant(O4,logarithm=TRUE)$modulus)+
  as.numeric(x4%*%solve(O4,t(x4))))
auto4<-try(copulaMaximiseJointDirect(x4,1,sd4,v4,4L,maxit=2L,
  withMu=FALSE,backend="auto"),silent=TRUE)
ok("E12 Gaussian common-Q helper stays finite in representable extreme tails",
   is.finite(lp4)&&abs(lp4-oracle4)<1e-8&&
     identical(attr(lp4,"pit_clipped"),0L)&&
     identical(attr(lp4,"density_floored"),0L))
ok("E12b auto optimizer routes the nested Gaussian null to native arithmetic",
   !inherits(auto4,"try-error")&&identical(auto4$backend,"r"))

## L-BFGS-B needs a finite penalty for invalid trial points, but that internal
## sentinel must never become a fitted exact objective.  A direct-parameter
## Gaussian vine still uses PITs; this deliberately saturated starting measure
## must therefore fail closed rather than report sentinel equality as a
## zero-gradient optimum.
pmDirect<-list(parameterMarginNormal(1),parameterMarginNormal(1))
directBad<-parameterMaximiseJoint(
  phi=matrix(c(40,0),nrow=1),w=1,margins0=pmDirect,vine0=vine,
  transform=c(0L,0L),X=matrix(1,nrow=1,ncol=2),locMap=diag(2),
  beta0=c(0,0),betaFree=integer(),maxit=2L,
  optimizeMargins=FALSE,optimizeVine=TRUE,diagnostics=TRUE)
ok("E13 an exact direct-Q sentinel can never be certified as an optimum",
   identical(directBad$conv,-1L)&&identical(directBad$value,-Inf)&&
     identical(directBad$gradientMax,Inf))

## Adjacent representable tails must remain admissible after the fail-closed
## sentinel repair. Re-evaluate the returned model through the literal exact
## density path so this also checks that the guard is not overbroad.
directGood<-parameterMaximiseJoint(
  phi=matrix(c(4,0),nrow=1),w=1,margins0=pmDirect,vine0=vine,
  transform=c(0L,0L),X=matrix(1,nrow=1,ncol=2),locMap=diag(2),
  beta0=c(0,0),betaFree=integer(),maxit=2L,
  optimizeMargins=FALSE,optimizeVine=TRUE,diagnostics=TRUE)
goodPredictor<-matrix(0,nrow=1,ncol=2)
goodExact<-try(parameterLogPrior(matrix(c(4,0),nrow=1),goodPredictor,
  directGood$vine,directGood$margins),silent=TRUE)
ok("E13b the direct-Q repair retains representable tail objectives",
   !inherits(goodExact,"try-error")&&is.finite(goodExact)&&
     directGood$conv!=-1L&&is.finite(directGood$value)&&
     is.finite(directGood$gradientMax)&&directGood$value>-1e99)

## Exact stationarity is a valid M-step: an optimizer that returns the starting
## point must not be converted into a failure merely because Q did not rise.
stationaryR<-matrix(c(1,1,1,-1,-1,1,-1,-1),ncol=2,byrow=TRUE)%*%
  chol(diag(c(.3,.4))%*%matrix(c(1,.35,.35,1),2,2)%*%diag(c(.3,.4)))
stationaryRfit<-copulaMaximiseJointDirectR(stationaryR,rep(.25,4),
  c(.3,.4),vine,2L,maxit=20L,withMu=TRUE,diagnostics=TRUE)
stationaryCppFit<-copulaMaximiseJointDirectCpp(stationaryR,rep(.25,4),
  c(.3,.4),vine,2L,maxit=20L,withMu=TRUE,diagnostics=TRUE)
stationaryParameterFit<-parameterMaximiseJoint(stationaryR,rep(.25,4),
  list(parameterMarginNormal(.3),parameterMarginNormal(.4)),vine,
  transform=c(0L,0L),X=matrix(1,nrow=4,ncol=2),locMap=diag(2),
  beta0=c(0,0),betaFree=integer(),maxit=20L,
  optimizeMargins=TRUE,optimizeVine=TRUE,diagnostics=TRUE)
ok("E14 a stationary no-move is retained as a valid exact common-Q endpoint",
   stationaryRfit$conv==0L&&stationaryCppFit$conv==0L&&
     stationaryParameterFit$conv==0L&&
     all(is.finite(c(stationaryRfit$value,stationaryCppFit$value,
       stationaryParameterFit$value,stationaryRfit$gradientMax,
       stationaryCppFit$gradientMax,stationaryParameterFit$gradientMax))))

## A Gaussian-surrogate GLS update is not the incumbent of the copula common-Q
## block. Once a joint location has been accepted, failure must return it rather
## than the transient ordinary-saemix value.
ok("E15 failed joint updates retain the accepted common-Q location incumbent",
   identical(copulaAcceptedBeta(c(9,8),c(1,2),
     matrix(c(1,0),nrow=2),betaFree=1L),c(1,8))&&
     identical(copulaAcceptedBeta(c(9,8),c(1,NA_real_),
       matrix(c(1,0),nrow=2),betaFree=1L),c(9,8))&&
     identical(copulaHeldBetaIndices(c(1L,2L),
       matrix(c(1,0),nrow=2),2L),1L))

hybridInvalid<-copulaMaximiseMargins(x,rep(.5,2),c(Inf,.4),vine,2L,
  maxit=1L,withMu=FALSE)
qHybrid0<-sum(rep(.25,4)*copulaLogPrior(stationaryR,vine,c(.3,.4)))
hybridFit<-copulaMaximiseJointHybrid(stationaryR,rep(.25,4),c(.3,.4),
  vine,2L,maxit=2L,withMu=FALSE,cycles=1L,polishMaxit=0L)
qHybrid1<-sum(rep(.25,4)*copulaLogPrior(
  sweep(stationaryR,2,hybridFit$delta,"-"),hybridFit$vine,hybridFit$sd))
ok("E16 hybrid margin proposals fail closed and preserve literal common Q",
   identical(hybridInvalid$conv,-1L)&&identical(hybridInvalid$value,-Inf)&&
     is.finite(qHybrid1)&&qHybrid1>=qHybrid0&&
   copulaQNondecreasing(1,1)&&copulaQNondecreasing(1+1e-12,1)&&
     !copulaQNondecreasing(1-1e-12,1)&&
     !copulaQNondecreasing(NA_real_,1))

truncated4<-rvinecopulib::vinecop_dist(list(list(
  rvinecopulib::bicop_dist("gaussian",parameters=.1),
  rvinecopulib::bicop_dist("gaussian",parameters=.2),
  rvinecopulib::bicop_dist("gaussian",parameters=.3))),
  rvinecopulib::dvine_structure(1:4,trunc_lvl=1L))
truncatedFlat<-copulaPadFlat(truncated4,4L)
truncatedLayout<-copulaParameterLayout(truncatedFlat)
truncatedRebuilt<-copulaBuildNative(truncatedFlat,truncatedLayout,
  truncatedLayout$par,truncated4,4L)
truncatedU<-matrix(c(.2,.3,.4,.5,.7,.6,.5,.4),ncol=4,byrow=TRUE)
ok("E17 rebuilding a truncated vine preserves its exact fixed-model identity",
   as.integer(dim(truncatedRebuilt)["trunc_lvl"])==1L&&
     length(truncatedRebuilt$pair_copulas)==1L&&
     identical(copulaFingerprint(truncatedRebuilt,4L),
       copulaFingerprint(truncated4,4L))&&
     identical(copulaParameterLayout(copulaPadFlat(truncatedRebuilt,4L))$par,
       truncatedLayout$par)&&
     max(abs(rvinecopulib::dvinecop(truncatedU,truncatedRebuilt)-
       rvinecopulib::dvinecop(truncatedU,truncated4)))<1e-14)

cat(sprintf("\n%d failure(s)\n", nFail))
if (nFail > 0L) quit(status = 1L)
