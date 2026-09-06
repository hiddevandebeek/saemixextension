## End-to-end direct lognormal/Gaussian-vine null versus ordinary saemix.
suppressMessages({library(devtools);load_all("C:/package/saemix-copula",quiet=TRUE)})
nFail<-0L
ok<-function(label,pass,extra="") {
  if(!isTRUE(pass))nFail<<-nFail+1L
  cat(sprintf("%-70s %s %s\n",label,if(isTRUE(pass))"PASS" else "**FAIL**",extra))
}
vGauss<-function(rho) rvinecopulib::vinecop_dist(
  list(list(rvinecopulib::bicop_dist("gaussian",parameters=rho))),
  rvinecopulib::dvine_structure(1:2))
set.seed(451);N<-75L;truth<-c(V=10,CL=3);sdTruth<-c(.26,.31);rhoTruth<-.55
pred<-matrix(rep(log(truth),each=N),N,2)
margins<-list(parameterMarginLognormal(sdTruth[1],"median"),
              parameterMarginLognormal(sdTruth[2],"median"))
psi<-parameterMarginsQuantile(rvinecopulib::rvinecop(N,vGauss(rhoTruth)),pred,margins)
times<-c(.25,.5,1,2,4,6,8,12);dose<-100
dd<-data.frame(id=rep(seq_len(N),each=length(times)),dose=dose,time=rep(times,N))
modelFun<-function(psi,id,xidep) {
  V<-psi[id,1];CL<-psi[id,2];xidep[,1]/V*exp(-CL/V*xidep[,2])
}
mu<-modelFun(psi,dd$id,cbind(dd$dose,dd$time));dd$y<-pmax(mu*(1+.06*rnorm(nrow(dd))),1e-7)
dat<-saemix::saemixData(name.data=dd,header=TRUE,name.group="id",
  name.predictors=c("dose","time"),name.response="y",verbose=FALSE)
makeModel<-function()saemix::saemixModel(model=modelFun,modeltype="structural",description="",
  psi0=matrix(truth,nrow=1,dimnames=list(NULL,names(truth))),transform.par=c(1,1),
  covariance.model=matrix(1,2,2),omega.init=diag(c(.25^2,.25^2)),
  error.model="proportional",verbose=FALSE)
ctl<-list(seed=451,save=FALSE,save.graphs=FALSE,print=FALSE,displayProgress=FALSE,
  nbiter.saemix=c(100,70),nbiter.mcmc=c(2,2,2,0),warnings=FALSE,ll.is=FALSE,fim=FALSE)
copulaClear();stock<-saemix::saemix(makeModel(),dat,ctl)
sd0<-sqrt(diag(stock@results@omega));rho0<-cov2cor(stock@results@omega)[1,2]
pop<-copulaPopulation(vGauss(rho0),
  margins=list(parameterMarginLognormal(sd0[1],"median"),
               parameterMarginLognormal(sd0[2],"median")),
  scale="parameter",proposalScale=sd0,mode="joint",activeFrom=81L,
  restartBurn=30L,refitEvery=5L,jointMaxit=15L,jointFinalMaxit=120L,
  guard=FALSE,populationAlgorithm="common-q")
direct<-saemix::saemix(makeModel(),dat,ctl,population=pop)
omDirect<-populationCovariance(direct,predictor=log(direct@results@fixed.psi),
                              scale="phi",nsim=80000L)
fixRel<-max(abs(direct@results@fixed.psi-stock@results@fixed.psi)/stock@results@fixed.psi)
omRel<-max(abs(omDirect-stock@results@omega)/pmax(abs(stock@results@omega),.02))
ok("G1 direct Gaussian null fixed effects match ordinary saemix",
   fixRel<.04,sprintf("max relative difference=%.3f",fixRel))
ok("G2 direct Gaussian null working-scale covariance matches stock",
   omRel<.16,sprintf("max scaled difference=%.3f",omRel))
ok("G3 Gaussian direct model has the same number of population parameters",
   direct@results@npar.est==stock@results@npar.est,
   sprintf("stock=%d direct=%d",stock@results@npar.est,direct@results@npar.est))
## The two declarations below induce the same density on phi.  With a common
## seed, their defensive-IS likelihood calculations must therefore agree draw
## for draw, not merely asymptotically.
directState<-list(vine=vGauss(rho0),
  margins=list(parameterMarginLognormal(sd0[1],"median"),
               parameterMarginLognormal(sd0[2],"median")),d=2L,
  populationScale="parameter",transform=c(1L,1L),proposalOmega=stock@results@omega)
etaState<-list(vine=vGauss(rho0),
  margins=list(copulaMarginNormal(sd0[1]),copulaMarginNormal(sd0[2])),
  sd=sd0,d=2L,populationScale="transformed-additive",
  transform=c(1L,1L),proposalOmega=stock@results@omega)
scoreAs<-function(fit,state) {
  attr(fit,"saemix.copula")<-copulaSnapshot(state,fit@model@indx.omega,
    fit@model@name.modpar[fit@model@indx.omega])
  slot(fit,"options")[["nmc.is"]]<-400L
  suppressWarnings(llisCopula.saemix(fit,defensive=.25,batch=100L,seed=902L))
}
llDirectDeclaration<-scoreAs(stock,directState)
llEtaDeclaration<-scoreAs(stock,etaState)
ok("G4 observed marginal likelihood is invariant to the equivalent declaration",
   abs(llDirectDeclaration@results@ll.is-llEtaDeclaration@results@ll.is)<1e-10,
   sprintf("absolute difference=%.3g",
     abs(llDirectDeclaration@results@ll.is-llEtaDeclaration@results@ll.is)))
## A low-level state used by objective-unit tests must not contaminate the
## documented population=NULL entry point.
copulaSet(vGauss(rho0),sd0,familySet=NULL,guard=FALSE)
ambientFingerprint<-copulaGet()$fingerprint
isolated<-saemix::saemix(makeModel(),dat,ctl,population=NULL)
ambientAfter<-copulaGet();ambientRestored<-copulaActive()&&
  identical(ambientAfter$fingerprint,ambientFingerprint)
copulaClear()
ok("G5 population=NULL is isolated from ambient low-level copula state",
   is.null(attr(isolated,"saemix.copula",exact=TRUE))&&ambientRestored&&
     identical(isolated@results@allpar,stock@results@allpar))
cat(sprintf("\n%d failure(s)\n",nFail));if(nFail>0)quit(status=1)
