## End-to-end direct heterogeneous margins plus non-Gaussian dependence.
suppressMessages({library(devtools);load_all("C:/package/saemix-copula",quiet=TRUE)})
nFail<-0L
ok<-function(label,pass,extra="") {
  if(!isTRUE(pass))nFail<<-nFail+1L
  cat(sprintf("%-72s %s %s\n",label,if(isTRUE(pass))"PASS" else "**FAIL**",extra))
}
makeV2<-function(g) rvinecopulib::vinecop_dist(
  list(list(rvinecopulib::bicop_dist("gumbel",parameters=g))),
  rvinecopulib::dvine_structure(1:2))

set.seed(141);N<-90L;truth<-c(V=10,CL=3);gTruth<-1.75
truthMargins<-list(parameterMarginGamma(9,"mean"),
                   parameterMarginLognormal(.30,"median"))
pred<-matrix(rep(log(truth),each=N),N,2)
psi<-parameterMarginsQuantile(rvinecopulib::rvinecop(N,makeV2(gTruth)),
                              pred,truthMargins)
colnames(psi)<-names(truth)
times<-c(.25,.5,1,2,4,6,8,12);dose<-100
dd<-data.frame(id=rep(seq_len(N),each=length(times)),dose=dose,
               time=rep(times,N))
modelFun<-function(psi,id,xidep) {
  V<-psi[id,1];CL<-psi[id,2]
  xidep[,1]/V*exp(-CL/V*xidep[,2])
}
mu<-modelFun(psi,dd$id,cbind(dd$dose,dd$time))
dd$y<-pmax(mu*(1+.06*rnorm(nrow(dd))),1e-7)
dat<-saemix::saemixData(name.data=dd,header=TRUE,name.group="id",
  name.predictors=c("dose","time"),name.response="y",verbose=FALSE)
mod<-saemix::saemixModel(model=modelFun,modeltype="structural",description="",
  psi0=matrix(truth,nrow=1,dimnames=list(NULL,names(truth))),
  transform.par=c(1,1),covariance.model=matrix(1,2,2),
  omega.init=diag(c(.25^2,.25^2)),error.model="proportional",verbose=FALSE)
ctl<-list(seed=141,save=FALSE,save.graphs=FALSE,print=FALSE,
  displayProgress=FALSE,nbiter.saemix=c(90,60),nbiter.mcmc=c(2,2,2,0),
  warnings=FALSE,ll.is=FALSE,fim=FALSE)
pop<-copulaPopulation(makeV2(1.25),
  margins=list(parameterMarginGamma(5,"mean"),
               parameterMarginLognormal(.45,"median")),
  scale="parameter",proposalScale=c(.25,.25),mode="joint",refitEvery=5L,
  jointMaxit=18L,jointFinalMaxit=140L,guard=FALSE)
fit<-saemix::saemix(mod,dat,ctl,population=pop);st<-copulaGet(fit)
gEst<-as.numeric(copulaPadFlat(st$vine,2)[[1]]$parameters)
shapeEst<-st$margins[[1]]$parameters["shape"]
sdlogEst<-st$margins[[2]]$parameters["sdlog"]
ok("C1 heterogeneous direct margins and Gumbel vine fit jointly",
  st$lastJoint$final && st$lastJoint$conv==0L &&
    st$lastJoint$nMarginParameters==2L && st$lastJoint$nCopulaParameters==1L)
ok("C2 natural population anchors remain close to truth",
  max(abs(fit@results@fixed.psi-truth)/truth)<.12,
  paste(sprintf("%s=%.3f",names(truth),fit@results@fixed.psi),collapse=" "))
ok("C3 direct margin shapes and native copula parameter recover",
  abs(shapeEst-9)<4 && abs(sdlogEst-.30)<.12 && abs(gEst-gTruth)<.40,
  sprintf("shape=%.2f sdlog=%.3f gumbel=%.3f",shapeEst,sdlogEst,gEst))
slot(fit,"options")[["nmc.is"]]<-600L
scored<-suppressWarnings(llisCopula.saemix(fit,defensive=.25,batch=100L,seed=19L))
ok("C4 fitted heterogeneous model has a finite observed marginal likelihood",
  is.finite(scored@results@ll.is) &&
    all(is.finite(attr(scored,"saemix.copula.likelihood",exact=TRUE)$ess)),
  sprintf("logLik=%.2f ESSmin=%.1f",scored@results@ll.is,
          attr(scored,"saemix.copula.likelihood",exact=TRUE)$ess_min))

cat(sprintf("\n%d failure(s)\n",nFail));if(nFail>0)quit(status=1)
