## General marginal laws belong to the same fixed-model Q as the vine.

suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet=TRUE)})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R")

nFail <- 0L
ok <- function(lbl, pass, extra = "") {
  if (!isTRUE(pass)) nFail <<- nFail + 1L
  cat(sprintf("%-64s %s %s\n", lbl, if (isTRUE(pass)) "PASS" else "**FAIL**", extra))
}

## F1: the legacy sd shortcut and explicit Normal objects are identical.
set.seed(901)
v2 <- etaVine(list(rvinecopulib::bicop_dist("gumbel", parameters = 1.7)), 2)
sd <- c(.3, .45); E <- sweep(matrix(rnorm(800), ncol = 2), 2, sd, "*")
mn <- lapply(sd, copulaMarginNormal)
l0 <- copulaLogPrior(E, v2, sd)
l1 <- copulaLogPrior(E, v2, margins = mn)
ok("F1 explicit Normal margins equal the legacy Gaussian shortcut",
   max(abs(l0 - l1)) < 1e-13, sprintf("maxdiff=%.3g", max(abs(l0-l1))))

## F2: built-in margins are normalized and their CDF/quantile maps agree.
mt <- copulaMarginStudent(.35, 5)
ml <- copulaMarginLaplace(.42)
it <- integrate(function(x) exp(mt$log_density(x, mt$parameters)), -Inf, Inf)$value
il <- integrate(function(x) exp(ml$log_density(x, ml$parameters)), -Inf, Inf)$value
u <- seq(.01, .99, length.out = 99)
pit <- mt$cdf(mt$quantile(u, mt$parameters), mt$parameters)
pil <- ml$cdf(ml$quantile(u, ml$parameters), ml$parameters)
ok("F2 Student and Laplace margins are normalized", max(abs(c(it,il)-1)) < 1e-7)
ok("F2b Student and Laplace PIT maps are inverse", max(abs(c(pit-u,pil-u))) < 1e-10)
ma <- copulaMarginDistribution("norm", c(mean=0,sd=.4), c(0,1e-6), c(0,1e3),
  free=c(FALSE,TRUE), centered=TRUE, scale=function(par) par["sd"],
  set_scale=function(par,value) { par["sd"]<-value; par })
la <- copulaMarginLayout(list(ma))
ok("F2c generic distribution adapter separates fixed and free parameters",
   identical(names(ma$parameters),c("mean","sd")) && length(la$par)==1L &&
     names(la$par)=="sd")
maFreeLocation <- copulaMarginDistribution("norm", c(mean=0,sd=.4),
  c(-1,1e-6),c(1,1e3),free=c(TRUE,TRUE),centered=TRUE,
  scale=function(par) par["sd"])
badLocation <- try(copulaSet(v2,margins=list(maFreeLocation,maFreeLocation)),
                   silent=TRUE)
copulaClear()
ok("F2d free marginal location cannot duplicate saemix X beta",
   inherits(badLocation,"try-error"))
badWorkingScale <- try(copulaSet(v2,margins=list(ma,ma)),silent=TRUE)
maSd <- copulaMarginDistribution("norm", c(mean=0,sd=.4), c(0,1e-6),
  c(0,1e3),free=c(FALSE,TRUE),centered=TRUE,
  scale=function(par) par["sd"],scale_is_sd=TRUE)
goodActualSd <- try(copulaSet(v2,margins=list(maSd,maSd),
  populationAlgorithm="common-q"),silent=TRUE)
copulaClear()
ok("F2e end-to-end custom margins require an actual finite SD declaration",
   inherits(badWorkingScale,"try-error") && !inherits(goodActualSd,"try-error"))

## F3: native marginal shapes/scales and copula dependence are jointly recovered.
set.seed(902)
truthMargins <- list(copulaMarginStudent(.30, 5), copulaMarginLaplace(.40))
truthVine <- etaVine(list(rvinecopulib::bicop_dist("gumbel", parameters = 1.8)), 2)
eta <- copulaMarginsQuantile(rvinecopulib::rvinecop(6000, truthVine), truthMargins)
startMargins <- list(copulaMarginStudent(.39, 15), copulaMarginLaplace(.31))
startVine <- etaVine(list(rvinecopulib::bicop_dist("gumbel", parameters = 1.25)), 2)
w <- rep(1/nrow(eta), nrow(eta))
q0 <- sum(w * copulaLogPrior(eta, startVine, margins = startMargins))
jm <- copulaMaximiseJointMargins(eta, w, startMargins, startVine, 2,
                                  maxit = 150L, withMu = FALSE)
q1 <- sum(w * copulaLogPrior(eta, jm$vine, margins = jm$margins))
p1 <- jm$margins[[1]]$parameters; p2 <- jm$margins[[2]]$parameters
g <- as.numeric(copulaPadFlat(jm$vine, 2)[[1]]$parameters)
ok("F3 common Q increases with non-Gaussian marginal parameters", q1 > q0 + .02,
   sprintf("dQ=%.4f",q1-q0))
ok("F3b Student SD and df are recovered jointly",
   abs(p1["sd"]-.30) < .03 && abs(p1["df"]-5) < 2,
   sprintf("sd=%.3f df=%.2f",p1["sd"],p1["df"]))
ok("F3c Laplace SD and Gumbel parameter are recovered jointly",
   abs(p2["sd"]-.40) < .03 && abs(g-1.8) < .15,
   sprintf("sd=%.3f g=%.3f",p2["sd"],g))

## F4: a mixed discrete/continuous vine uses CDF left limits and probability
## mass, not a continuous density at a discrete PIT.
mb <- copulaMargin("bernoulli", c(prob=.3), 1e-5, 1-1e-5,
  log_density=function(x,par) dbinom(x,1,par[1],log=TRUE),
  cdf=function(x,par) pbinom(x,1,par[1]),
  cdf_left=function(x,par) pbinom(x-1,1,par[1]),
  quantile=function(u,par) qbinom(u,1,par[1]),
  random=function(n,par) rbinom(n,1,par[1]), type="discrete",
  scale=function(par) sqrt(par[1]*(1-par[1])), centered=FALSE)
mi <- list(mb, copulaMarginNormal(.5))
vi <- rvinecopulib::vinecop_dist(
  list(list(rvinecopulib::bicop_dist("indep", var_types=c("d","c")))),
  rvinecopulib::dvine_structure(1:2), var_types=c("d","c"))
set.seed(903); X <- cbind(rbinom(1000,1,.3), rnorm(1000,sd=.5))
lmix <- copulaLogPrior(X, vi, margins=mi)
loracle <- dbinom(X[,1],1,.3,log=TRUE)+dnorm(X[,2],sd=.5,log=TRUE)
ok("F4 mixed discrete/continuous independence likelihood is exact",
   max(abs(lmix-loracle)) < 1e-12,
   sprintf("maxdiff=%.3g",max(abs(lmix-loracle))))

rho <- .45
vg <- rvinecopulib::vinecop_dist(
  list(list(rvinecopulib::bicop_dist("gaussian",parameters=rho,
                                     var_types=c("d","c")))),
  rvinecopulib::dvine_structure(1:2),var_types=c("d","c"))
z <- X[,2]/.5; cut <- qnorm(.7)
p0givenz <- pnorm((cut-rho*z)/sqrt(1-rho^2))
ldependent <- dnorm(X[,2],sd=.5,log=TRUE) +
  ifelse(X[,1]==0,log(p0givenz),log1p(-p0givenz))
lmixDependent <- copulaLogPrior(X,vg,margins=mi)
ok("F4c dependent Bernoulli/Normal mixed likelihood matches rectangle oracle",
   max(abs(lmixDependent-ldependent))<5e-7,
   sprintf("maxdiff=%.3g",max(abs(lmixDependent-ldependent))))

badCenter <- try(copulaSet(vi, margins=mi), silent=TRUE)
copulaClear()
ok("F4b uncentred eta margins are rejected by the identified model API",
   inherits(badCenter, "try-error"))

## F5: a frozen block must be held fixed inside the common-Q optimization,
## rather than optimized provisionally and discarded afterward.
fm <- copulaMaximiseJointMargins(eta, w, startMargins, startVine, 2,
  maxit=40L, withMu=FALSE, optimizeMargins=FALSE, optimizeVine=TRUE)
fv <- copulaMaximiseJointMargins(eta, w, startMargins, startVine, 2,
  maxit=40L, withMu=FALSE, optimizeMargins=TRUE, optimizeVine=FALSE)
sameMargins <- all(vapply(seq_along(startMargins), function(j)
  identical(fm$margins[[j]]$parameters, startMargins[[j]]$parameters), logical(1)))
sameVine <- identical(copulaFingerprint(fv$vine, 2),
                      copulaFingerprint(startVine, 2))
ok("F5 frozen margins are held fixed inside the optimizer", sameMargins)
ok("F5b frozen vine parameters are held fixed inside the optimizer", sameVine)

fixedNormal <- copulaMarginNormal(.30); fixedNormal$free[] <- FALSE
fn <- copulaMaximiseJoint(E, rep(1/nrow(E),nrow(E)), c(.30,.30), v2, 2,
  maxit=20L, withMu=FALSE,
  margins=list(fixedNormal, fixedNormal), optimizeVine=FALSE)
ok("F5c fixed Normal masks bypass the all-free Gaussian fast path",
   identical(unname(fn$sd), c(.30,.30)) && fn$nMarginParameters==0L)

copulaSet(v2, margins=startMargins, activeFrom=5L, freezeSd=TRUE,
          freezeVine=TRUE, guard=FALSE, populationAlgorithm="common-q")
before <- copulaGet()
copulaActivateFromGaussian(5L, diag(c(.8,.9)^2))
after <- copulaGet(); copulaClear()
sameActivatedVine <- identical(
  lapply(copulaPadFlat(before$vine,2), `[[`, "parameters"),
  lapply(copulaPadFlat(after$vine,2), `[[`, "parameters"))
sameActivatedMargins <- all(vapply(seq_along(before$margins), function(j)
  identical(before$margins[[j]]$parameters, after$margins[[j]]$parameters), logical(1)))
ok("F6 staged activation honors genuinely frozen margin and vine blocks",
   sameActivatedVine && sameActivatedMargins)

copulaSet(v2, sd=c(.3,.4), mode="joint", guard=FALSE,
  populationAlgorithm="common-q")
Erank <- matrix(rnorm(80),ncol=2); Xrank <- cbind(1,1)
Xrank <- Xrank[rep(1,nrow(Erank)),,drop=FALSE]
badRank <- try(copulaPoolUpdate(Erank, gamma=1, nchains=1,
  phiM=Erank, xM=Xrank, locMap=rbind(c(1,0),c(1,0)), beta=c(0,0),
  betaFree=1:2), silent=TRUE)
copulaClear()
ok("F7 rank-deficient beta-to-eta location designs are rejected",
   inherits(badRank,"try-error"))

cat(sprintf("\n%d failure(s)\n", nFail))
if(nFail>0) quit(status=1)
