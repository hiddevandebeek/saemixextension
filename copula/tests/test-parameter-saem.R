## End-to-end saemix fit with a direct natural-parameter population margin.
suppressMessages({library(devtools);load_all("C:/package/saemix-copula",quiet=TRUE)})
nFail<-0L
ok<-function(label,pass,extra="") {
  if(!isTRUE(pass))nFail<<-nFail+1L
  cat(sprintf("%-70s %s %s\n",label,if(isTRUE(pass))"PASS" else "**FAIL**",extra))
}

set.seed(919); N<-70L; shapeTruth<-10; meanTruth<-8; sigmaTruth<-.7
psiI<-stats::rgamma(N,shapeTruth,scale=meanTruth/shapeTruth)
dd<-data.frame(id=rep(seq_len(N),each=4L),time=rep(seq_len(4L),N))
dd$y<-psiI[dd$id]+stats::rnorm(nrow(dd),sd=sigmaTruth)
f<-function(psi,id,xidep) psi[id,1]
dat<-saemix::saemixData(name.data=dd,header=TRUE,name.group="id",
  name.predictors="time",name.response="y",verbose=FALSE)
mod<-saemix::saemixModel(model=f,modeltype="structural",description="",
  psi0=matrix(meanTruth,nrow=1,dimnames=list(NULL,"level")),transform.par=1,
  covariance.model=matrix(1,1,1),omega.init=matrix(.20^2,1,1),
  error.model="constant",verbose=FALSE)
ctl<-list(seed=919,save=FALSE,save.graphs=FALSE,print=FALSE,
  displayProgress=FALSE,nbiter.saemix=c(70,50),nbiter.mcmc=c(2,2,1,0),
  warnings=FALSE,ll.is=FALSE,fim=FALSE)
v1<-rvinecopulib::vinecop_dist(list(),rvinecopulib::dvine_structure(1))
copulaClear()
population<-copulaPopulation(v1,margins=list(parameterMarginGamma(5,"mean")),
  scale="parameter",proposalScale=.25,mode="joint",
  refitEvery=4L,jointMaxit=25L,jointFinalMaxit=120L,guard=FALSE)
fit<-saemix::saemix(mod,dat,ctl,population=population)
st<-copulaGet(fit)

estMean<-as.numeric(fit@results@fixed.psi)
estShape<-st$margins[[1]]$parameters["shape"]
ok("S1 direct parameter population survives the standard saemix path",
  identical(st$populationScale,"parameter") && st$lastJoint$final &&
    st$lastJoint$nMarginParameters==1L)
ok("S2 direct Gamma mean and shape are recovered",
  abs(estMean-meanTruth)<.45 && abs(estShape-shapeTruth)<3.5,
  sprintf("mean=%.3f shape=%.3f",estMean,estShape))
ok("S3 fitted state stores proposal covariance separately",
  inherits(attr(fit,"saemix.copula",exact=TRUE),"saemixCopulaSnapshot") &&
    all(is.finite(st$proposalOmega)))
ok("S3b direct fits do not mislabel proposal covariance as population omega",
  all(is.na(fit@results@omega[fit@model@indx.omega,fit@model@indx.omega])) &&
    all(is.finite(fit@results@cond.mean.psi)))
pcov<-populationCovariance(fit,predictor=log(meanTruth),scale="psi",nsim=20000L)
ok("S3c conditional natural-scale population covariance is available explicitly",
  is.finite(pcov[1,1]) && abs(pcov[1,1]-meanTruth^2/estShape)<.8,
  sprintf("variance=%.3f",pcov[1,1]))
ok("S3d fitted state retains native direct-population convergence traces",
  length(st$trace)>20L && isTRUE(tail(st$trace,1)[[1]]$final) &&
    all(c("beta","margin","copula","tau","gamma","guardTriggered") %in%
      names(tail(st$trace,1)[[1]])) &&
    length(attr(fit,"saemix.copula",exact=TRUE)$trace)==length(st$trace),
  sprintf("updates=%d",length(st$trace)))
ok("S4 package-native population argument restores global state after fitting",
  !copulaActive() && identical(
    copulaFingerprint(v1,1L,population$arguments$margins),st$fingerprint))

slot(fit,"options")[["nmc.is"]]<-500L
scored<-suppressWarnings(llisCopula.saemix(fit,defensive=.25,batch=100L,seed=31L))
isd<-attr(scored,"saemix.copula.likelihood",exact=TRUE)
ok("S5 observed marginal likelihood works after copulaClear",
  is.finite(scored@results@ll.is) && isd$draws_used==500L &&
    identical(isd$engine,"parameter-copula-defensive-is") &&
    isd$pit_clipped==0L && isd$density_floored==0L,
  sprintf("logLik=%.3f ESSmin=%.1f",scored@results@ll.is,isd$ess_min))
tmp<-tempfile(fileext=".rds");saveRDS(scored,tmp);restored<-readRDS(tmp);unlink(tmp)
rescored<-suppressWarnings(llisCopula.saemix(restored,defensive=.25,batch=100L,seed=31L))
ok("S6 direct population and marginal likelihood survive serialization",
  identical(scored@results@ll.is,rescored@results@ll.is))
simPsi<-simulateIndividualParameters(restored,nsim=300L,seed=77L)
ok("S7 standard saemix simulation draws from the persisted direct population",
  all(is.finite(simPsi)) && abs(mean(simPsi[,1])-estMean)<.12 &&
    abs(stats::var(simPsi[,1])-estMean^2/estShape)<.35,
  sprintf("mean=%.3f variance=%.3f",mean(simPsi[,1]),stats::var(simPsi[,1])))
mapped<-map.saemix(restored)
ok("S8 MAP estimation and individual predictions use the persisted direct prior",
  nrow(mapped@results@map.psi)==N && all(is.finite(mapped@results@map.psi[,1])) &&
    all(is.na(mapped@results@map.shrinkage)))
conditional<-conddist.saemix(restored,nsamp=3L,max.iter=50L,plot=FALSE)
ok("S9 conditional-distribution sampling uses the persisted direct posterior",
  identical(dim(conditional@results@psi.samp),c(N,1L,3L)) &&
    all(is.finite(conditional@results@cond.mean.psi)) &&
    all(conditional@results@cond.var.phi[,1]>0))

cat(sprintf("\n%d failure(s)\n",nFail));if(nFail>0)quit(status=1)
