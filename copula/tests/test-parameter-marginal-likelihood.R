## Mathematical/numerical oracles for the observed marginal likelihood.
suppressMessages({library(devtools); load_all("C:/package/saemix-copula",quiet=TRUE)})

nFail<-0L
ok<-function(label,pass,extra="") {
  if(!isTRUE(pass)) nFail<<-nFail+1L
  cat(sprintf("%-72s %s %s\n",label,if(isTRUE(pass))"PASS" else "**FAIL**",extra))
}
v1<-rvinecopulib::vinecop_dist(list(),rvinecopulib::dvine_structure(1))

## Analytic Gaussian convolution: this is the marginal likelihood, not merely
## a normalization check.
y<-1.35; mu<-.4; tau<-.65; sigma<-.28
mN<-parameterMarginNormal(tau,"identity")
llQ<-log(integrate(function(psi)
  exp(stats::dnorm(y,psi,sigma,log=TRUE)+
      parameterLogPrior(matrix(psi,ncol=1),matrix(mu,nrow=length(psi)),v1,list(mN))),
  -Inf,Inf,rel.tol=1e-10)$value)
llA<-stats::dnorm(y,mu,sqrt(tau^2+sigma^2),log=TRUE)
ok("ML1 direct Normal observed likelihood equals analytic convolution",
   abs(llQ-llA)<2e-8,sprintf("error=%.3g",llQ-llA))

## A nonlinear transformation must give the same integrated likelihood in psi
## and phi only when the working-coordinate Jacobian is present.
y<-8.2; pred<-log(7.5); sigma<-.7
mL<-parameterMarginLognormal(.35,"median")
llPsi<-log(integrate(function(psi)
  exp(stats::dnorm(y,psi,sigma,log=TRUE)+
      parameterLogPrior(matrix(psi,ncol=1),matrix(pred,nrow=length(psi)),v1,list(mL))),
  0,Inf,rel.tol=2e-9)$value)
state<-list(populationScale="parameter",d=1L,vine=v1,margins=list(mL),cores=1L)
llPhi<-log(integrate(function(phi)
  exp(stats::dnorm(y,exp(phi),sigma,log=TRUE)-
      parameterUphi(matrix(phi,ncol=1),matrix(pred,nrow=length(phi)),1L,state)),
  -Inf,Inf,rel.tol=2e-9)$value)
ok("ML2 natural- and working-scale marginal likelihoods agree",
   abs(llPsi-llPhi)<2e-8,sprintf("error=%.3g",llPsi-llPhi))

## The median-anchored lognormal is exactly the stock additive Normal model on
## log(phi), point by point, including the forward Jacobian.
grid<-seq(-1.5,2.0,length.out=101); sdlog<-.42; pred<-.3
mLN<-parameterMarginLognormal(sdlog,"median")
stLN<-list(populationScale="parameter",d=1L,vine=v1,margins=list(mLN),cores=1L)
lpPhi<--parameterUphi(matrix(grid,ncol=1),matrix(pred,nrow=length(grid)),1L,stLN)
lpStock<-stats::dnorm(grid,pred,sdlog,log=TRUE)
ok("ML3 transformed Gaussian nested null agrees pointwise",
   max(abs(lpPhi-lpStock))<2e-12,
   sprintf("max error=%.3g",max(abs(lpPhi-lpStock))))

## Fisher identity for a non-Gaussian direct Gamma population.  The observed
## score is differentiated outside the integral; the complete score is
## averaged under the exact quadrature posterior.
y<-9.4; sigma<-1.1; beta<-log(8.7); shape<-7.3
obsLL<-function(beta,shape) {
  mg<-parameterMarginGamma(shape,"mean")
  log(integrate(function(psi)
    exp(stats::dnorm(y,psi,sigma,log=TRUE)+
        parameterLogPrior(matrix(psi,ncol=1),matrix(beta,nrow=length(psi)),v1,list(mg))),
    0,Inf,rel.tol=2e-9)$value)
}
h<-1e-4
scoreObs<-c(beta=(obsLL(beta+h,shape)-obsLL(beta-h,shape))/(2*h),
            shape=(obsLL(beta,shape+h)-obsLL(beta,shape-h))/(2*h))
mg<-parameterMarginGamma(shape,"mean")
Li<-exp(obsLL(beta,shape))
posteriorScore<-function(which) integrate(function(psi) {
  lpPlus<-if(which=="beta")
    parameterLogPrior(matrix(psi,ncol=1),matrix(beta+h,nrow=length(psi)),v1,list(mg)) else
    parameterLogPrior(matrix(psi,ncol=1),matrix(beta,nrow=length(psi)),v1,
      list(parameterMarginGamma(shape+h,"mean")))
  lpMinus<-if(which=="beta")
    parameterLogPrior(matrix(psi,ncol=1),matrix(beta-h,nrow=length(psi)),v1,list(mg)) else
    parameterLogPrior(matrix(psi,ncol=1),matrix(beta,nrow=length(psi)),v1,
      list(parameterMarginGamma(shape-h,"mean")))
  sc<-(lpPlus-lpMinus)/(2*h)
  exp(stats::dnorm(y,psi,sigma,log=TRUE)+
      parameterLogPrior(matrix(psi,ncol=1),matrix(beta,nrow=length(psi)),v1,list(mg)))*sc/Li
},0,Inf,rel.tol=2e-8)$value
scoreFI<-c(beta=posteriorScore("beta"),shape=posteriorScore("shape"))
ok("ML4 Fisher identity holds for beta and Gamma shape",
   max(abs(scoreObs-scoreFI))<2e-5,
   sprintf("max error=%.3g",max(abs(scoreObs-scoreFI))))

## Common-particle Q recovers direct population parameters without using EBEs.
set.seed(812); N<-2500L; betaTruth<-log(10); shapeTruth<-9
psi<-stats::rgamma(N,shapeTruth,scale=10/shapeTruth); phi<-matrix(log(psi),ncol=1)
X<-matrix(1,N,1); locMap<-matrix(1,1,1); w<-rep(1/N,N)
fitQ<-parameterMaximiseJoint(phi,w,list(parameterMarginGamma(4,"mean")),v1,1L,
  X,locMap,beta0=log(8),betaFree=1L,maxit=200L,optimizeVine=FALSE)
gradQ<-function(p) {
  mm<-parameterMarginGamma(unname(p[2]),"mean")
  pred<-matrix(p[1],N,1)
  -mean(parameterLogPrior(exp(phi),pred,v1,list(mm))+parameterLogJacobian(phi,1L))
}
g<-c((gradQ(c(fitQ$beta[1]+h,fitQ$margins[[1]]$parameters[1]))-
      gradQ(c(fitQ$beta[1]-h,fitQ$margins[[1]]$parameters[1])))/(2*h),
     (gradQ(c(fitQ$beta[1],fitQ$margins[[1]]$parameters[1]+h))-
      gradQ(c(fitQ$beta[1],fitQ$margins[[1]]$parameters[1]-h)))/(2*h))
ok("ML5 joint direct-Q optimizer is stationary and recovers its truth",
   max(abs(g))<2e-4 && abs(exp(fitQ$beta[1])-10)<.2 &&
     abs(fitQ$margins[[1]]$parameters[1]-shapeTruth)<1.2,
   sprintf("mean=%.3f shape=%.3f |grad|=%.3g",exp(fitQ$beta[1]),
     fitQ$margins[[1]]$parameters[1],max(abs(g))))

## Native pair-copula parameters are part of the same direct Q, rather than a
## post-hoc fit to EBEs or marginal ranks.
set.seed(813); N2<-3000L; gTruth<-1.85
makeV2<-function(g) rvinecopulib::vinecop_dist(
  list(list(rvinecopulib::bicop_dist("gumbel",parameters=g))),
  rvinecopulib::dvine_structure(1:2))
vTruth<-makeV2(gTruth); u<-rvinecopulib::rvinecop(N2,vTruth)
m2<-list(parameterMarginGamma(8,"mean"),parameterMarginLognormal(.32,"median"))
for(j in 1:2)m2[[j]]$free[]<-FALSE
beta2<-log(c(10,4)); pred2<-matrix(rep(beta2,each=N2),N2,2)
psi2<-parameterMarginsQuantile(u,pred2,m2); phi2<-log(psi2)
fitC<-parameterMaximiseJoint(phi2,rep(1/N2,N2),m2,makeV2(1.2),c(1L,1L),
  matrix(1,N2,2),diag(2),beta2,betaFree=integer(),maxit=160L,
  optimizeMargins=FALSE,optimizeVine=TRUE)
gEst<-as.numeric(copulaPadFlat(fitC$vine,2)[[1]]$parameters)
ok("ML6 native copula parameter is optimized inside the common direct Q",
   fitC$conv==0L && abs(gEst-gTruth)<.12,
   sprintf("truth=%.3f estimate=%.3f",gTruth,gEst))

## Fisher identity also covers a native copula parameter.  Common uniform
## quadrature nodes make the observed-likelihood derivative and posterior
## expected complete score independently reproducible.
set.seed(814); nq<-120000L; uq<-matrix(runif(2*nq,1e-5,1-1e-5),nq,2)
predq<-matrix(rep(beta2,each=nq),nq,2)
psiq<-parameterMarginsQuantile(uq,predq,m2)
logResponse<-stats::dnorm(12.5,psiq[,1]+psiq[,2],1.4,log=TRUE)
mcLL<-function(g) {
  lc<-log(rvinecopulib::dvinecop(uq,makeV2(g)))
  mx<-max(logResponse+lc); log(mean(exp(logResponse+lc-mx)))+mx
}
hg<-2e-4
scoreCopObs<-(mcLL(gTruth+hg)-mcLL(gTruth-hg))/(2*hg)
lc0<-log(rvinecopulib::dvinecop(uq,makeV2(gTruth)))
sc<-(log(rvinecopulib::dvinecop(uq,makeV2(gTruth+hg)))-
     log(rvinecopulib::dvinecop(uq,makeV2(gTruth-hg))))/(2*hg)
lw<-logResponse+lc0;ww<-exp(lw-max(lw));scoreCopFI<-sum(ww*sc)/sum(ww)
ok("ML7 Fisher identity holds for the native copula parameter",
   abs(scoreCopObs-scoreCopFI)<2e-6,
   sprintf("error=%.3g",scoreCopObs-scoreCopFI))

for(tr in c(2L,3L)) {
  link<-if(tr==2L)"probit" else "logit";mb<-parameterMarginBeta(14,link)
  stb<-list(populationScale="parameter",d=1L,vine=v1,margins=list(mb),cores=1L)
  normPhi<-integrate(function(phi)exp(-parameterUphi(
    matrix(phi,ncol=1),matrix(.25,nrow=length(phi)),tr,stb)),
    -Inf,Inf,rel.tol=2e-9)$value
  ok(paste0("ML8 bounded transform ",tr," includes the exact Jacobian"),
     abs(normPhi-1)<2e-8,sprintf("integral=%.10f",normPhi))
}
## The stock probit helper historically floored the lower tail at 1e-30.
## A direct-natural density with the exact Normal Jacobian cannot use that
## many-to-one numerical map: the density would no longer normalize exactly.
## Test a representable tail below the former floor pointwise, because ordinary
## quadrature cannot resolve a discrepancy of order 1e-30 in the total mass.
phiTail<--12;psiTail<-stats::pnorm(phiTail)
mbTail<-parameterMarginBeta(14,"probit")
stTail<-list(populationScale="parameter",d=1L,vine=v1,
  margins=list(mbTail),cores=1L,numericalPolicy="exact")
lpTailWorking<--parameterUphi(matrix(phiTail,ncol=1),matrix(.25,nrow=1),
  2L,stTail)
lpTailNatural<-parameterLogPrior(matrix(psiTail,ncol=1),matrix(.25,nrow=1),
  v1,list(mbTail),numericalPolicy="exact")+stats::dnorm(phiTail,log=TRUE)
ok("ML8c probit direct-margin map has no target-changing lower-tail floor",
   abs(as.numeric(transphi(matrix(phiTail,ncol=1),2L))/psiTail-1)<1e-14&&
     abs(as.numeric(lpTailWorking-lpTailNatural))<1e-12,
   sprintf("psi=%.3g log-density error=%.3g",psiTail,
     as.numeric(lpTailWorking-lpTailNatural)))
mgSaturation<-parameterMarginGamma(8,"mean")
copulaSet(v1,margins=list(mgSaturation),populationScale="parameter",
  proposalScale=.3,guard=FALSE)
saturated<-parameterUphi(matrix(-1000,ncol=1),matrix(log(8),ncol=1),1L)
saturationState<-copulaGet();copulaClear()
ok("ML8d finite-precision support-transform saturation is rejected and counted",
   is.infinite(saturated)&&saturated>0&&
     identical(attr(saturated,"transform_invalid"),1L)&&
     identical(saturationState$transformInvalid,1L))
probitGrid<-matrix(c(-2,0,2),ncol=1)
ok("ML8e probit delta-method derivative is the derivative of Phi",
   max(abs(dtransphi(probitGrid,2L)-stats::dnorm(probitGrid)))<1e-15&&
     max(abs(derivphi(probitGrid,2L)-1/stats::dnorm(probitGrid)))<1e-12)
lpGuard<-parameterLogPrior(matrix(c(7,9,11),ncol=1),
  matrix(log(9),nrow=3),v1,list(parameterMarginGamma(8,"mean")))
ok("ML9 clipping and density floors are inactive in likelihood oracles",
   identical(attr(lpGuard,"pit_clipped"),0L) &&
     identical(attr(lpGuard,"density_floored"),0L))

## Every built-in continuous direct margin must define the same observed
## marginal likelihood in its declared natural coordinate and in saemix's
## internal Markov coordinate.  This checks normalization, support, link, and
## Jacobian together rather than checking the CDF callbacks in isolation.
families<-list(
  normal=list(m=parameterMarginNormal(.7),pred=.4,tr=0L,lo=-Inf,hi=Inf,y=.8,s=.5),
  student=list(m=parameterMarginStudent(.7,5),pred=.4,tr=0L,lo=-Inf,hi=Inf,y=.8,s=.5),
  laplace=list(m=parameterMarginLaplace(.7),pred=.4,tr=0L,lo=-Inf,hi=Inf,y=.8,s=.5),
  lognormal=list(m=parameterMarginLognormal(.35,"median"),pred=log(2),tr=1L,
    lo=0,hi=Inf,y=2.3,s=.5),
  gamma=list(m=parameterMarginGamma(7,"mean"),pred=log(2),tr=1L,
    lo=0,hi=Inf,y=2.3,s=.5),
  weibull=list(m=parameterMarginWeibull(3,"mean"),pred=log(2),tr=1L,
    lo=0,hi=Inf,y=2.3,s=.5),
  beta=list(m=parameterMarginBeta(15,"logit"),pred=stats::qlogis(.4),tr=3L,
    lo=0,hi=1,y=.55,s=.15))
familyErrors<-vapply(families,function(z) {
  st<-list(populationScale="parameter",d=1L,vine=v1,margins=list(z$m),cores=1L)
  natural<-integrate(function(psi)exp(stats::dnorm(z$y,psi,z$s,log=TRUE)+
    parameterLogPrior(matrix(psi,ncol=1),matrix(z$pred,nrow=length(psi)),
      v1,list(z$m))),z$lo,z$hi,rel.tol=2e-8)$value
  internal<-integrate(function(phi) {
    psi<-transphi(matrix(phi,ncol=1),z$tr)[,1]
    exp(stats::dnorm(z$y,psi,z$s,log=TRUE)-parameterUphi(
      matrix(phi,ncol=1),matrix(z$pred,nrow=length(phi)),z$tr,st))
  },-Inf,Inf,rel.tol=2e-8)$value
  abs(log(natural)-log(internal))
},numeric(1))
ok("ML10 every built-in margin preserves the observed marginal likelihood",
   max(familyErrors)<2e-7,
   sprintf("max error=%.3g (%s)",max(familyErrors),names(which.max(familyErrors))))

## Fisher's identity covers response parameters as well as the population
## block.  This is the score for the residual standard deviation in ML4.
## Re-evaluate because obsLL closes over sigma.
obsSigma<-function(s) {
  mg<-parameterMarginGamma(shape,"mean")
  log(integrate(function(psi)exp(stats::dnorm(y,psi,s,log=TRUE)+
    parameterLogPrior(matrix(psi,ncol=1),matrix(beta,nrow=length(psi)),v1,list(mg))),
    0,Inf,rel.tol=2e-9)$value)
}
scoreSigmaObs<-(obsSigma(sigma+h)-obsSigma(sigma-h))/(2*h)
LiSigma<-exp(obsSigma(sigma))
scoreSigmaFI<-integrate(function(psi) {
  completeScore<--1/sigma+(y-psi)^2/sigma^3
  exp(stats::dnorm(y,psi,sigma,log=TRUE)+parameterLogPrior(
    matrix(psi,ncol=1),matrix(beta,nrow=length(psi)),v1,
    list(parameterMarginGamma(shape,"mean"))))*completeScore/LiSigma
},0,Inf,rel.tol=2e-8)$value
ok("ML11 Fisher identity holds for a response/residual parameter",
   abs(scoreSigmaObs-scoreSigmaFI)<2e-5,
   sprintf("error=%.3g",scoreSigmaObs-scoreSigmaFI))

cat(sprintf("\n%d failure(s)\n",nFail)); if(nFail>0)quit(status=1)
