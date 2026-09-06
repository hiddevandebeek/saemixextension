## Basic production workflow: Gaussian chain, shared-incumbent-Q R-vine
## learning, then a fresh fixed-model SAEM-MCMC run.
suppressMessages({library(devtools);load_all("C:/package/saemix-copula",quiet=TRUE)})
nFail<-0L
ok<-function(label,pass,extra=""){
  if(!isTRUE(pass))nFail<<-nFail+1L
  cat(sprintf("%-72s %s %s\n",label,if(isTRUE(pass))"PASS" else "**FAIL**",extra))
}

set.seed(451);N<-35L;times<-c(0,1,2,3)
b<-rvinecopulib::bicop_dist("gumbel",parameters=1.7)
truthVine<-rvinecopulib::vinecop_dist(list(list(b)),
  rvinecopulib::dvine_structure(1:2))
u<-rvinecopulib::rvinecop(N,truthVine)
eta<-sweep(qnorm(u),2,c(.30,.20),"*")
theta<-c(2,.5);psiI<-sweep(eta,2,theta,"+")
dd<-data.frame(id=rep(seq_len(N),each=length(times)),
  time=rep(times,N))
dd$y<-psiI[dd$id,1]+psiI[dd$id,2]*dd$time+rnorm(nrow(dd),sd=.12)
f<-function(psi,id,xidep)psi[id,1]+psi[id,2]*xidep[,1]
dat<-saemix::saemixData(name.data=dd,header=TRUE,name.group="id",
  name.predictors="time",name.response="y",verbose=FALSE)
mod<-saemix::saemixModel(model=f,modeltype="structural",description="",
  psi0=matrix(theta,nrow=1,dimnames=list(NULL,c("intercept","slope"))),
  transform.par=c(0,0),covariance.model=matrix(1,2,2),
  omega.init=diag(c(.25^2,.15^2)),error.model="constant",verbose=FALSE)
ctl<-list(seed=451,save=FALSE,save.graphs=FALSE,print=FALSE,
  displayProgress=FALSE,nbiter.saemix=c(35,20),nbiter.mcmc=c(2,2,1,0),
  map=FALSE,warnings=FALSE,ll.is=FALSE,fim=FALSE)
fit<-suppressWarnings(saemixCopula(mod,dat,ctl,nsamp=3L,max.iter=40L,
  familySet=c("gaussian","gumbel","clayton","frank"),
  learningSeed=9173L,learningValidationSeed=27183L,
  populationControl=list(jointMaxit=15L,jointFinalMaxit=80L,guard=FALSE)))
learning<-attr(fit,"saemix.copula.learning",exact=TRUE)
gaussian<-attr(fit,"saemix.gaussian.fit",exact=TRUE)
state<-copulaGet(fit)
ok("W1 workflow returns the ordinary fitted saemix object API",
  inherits(fit,"SaemixObject")&&identical(fit@results@status,"fitted"))
ok("W2 ordinary Gaussian chain is preserved for comparison",
  inherits(gaussian,"SaemixObject")&&
    is.null(attr(gaussian,"saemix.copula",exact=TRUE)))
ok("W3 learning uses several one-draw-per-subject posterior slices",
  inherits(learning,"saemixCopulaLearning")&&
    identical(dim(learning$posterior$eta),c(N,2L,3L)))
ok("W4 automatic and forced-family proposals receive a common-Q ranking",
  nrow(learning$proposals)>3L&&nrow(learning$ranking)>=1L&&
    all(c("Q","QAIC_step","QBIC_step","trainingDelta","stationary",
      "selected","modelChanged")%in%
      names(learning$ranking))&&
    any(learning$proposals$source=="forced-family")&&
    any(learning$proposals$source=="incumbent")&&
    is.list(learning$commonQ$validation)&&
    identical(learning$commonQ$trainingSeed,9173L)&&
    identical(learning$commonQ$validationSeed,27183L)&&
    inherits(learning$vine,"vinecop_dist"))
ok("W5 legacy workflow uses an explicit fixed-fingerprint common-Q refit",
  isTRUE(state$modelFrozen)&&isTRUE(state$lastJoint$final)&&
    state$lastJoint$conv==0L&&
    identical(state$populationAlgorithm,"common-q")&&
    identical(state$fingerprint,copulaFingerprint(state$vine,2L,state$margins)),
  sprintf("algorithm=%s",state$populationAlgorithm))

gCov<-gaussian;mCov<-mod
mCov@betaest.model<-rbind(c(1,1),c(1,1))
mCov@psi0<-rbind(mCov@psi0,c(0,0))
gCov@results@betaC<-c(.11,-.07)
warmed<-copulaWarmStartModel(mCov,gCov)
ok("W6 warm start transfers fitted covariate coefficients as well as intercepts",
  max(abs(warmed@psi0[1,]-gCov@results@fixed.psi))<1e-12&&
    max(abs(warmed@psi0[2,]-c(.11,-.07)))<1e-12)

badMod<-mod;badMod@transform.par<-c(1,0)
badExternal<-try(saemixCopula(badMod,dat,ctl,gaussianFit=gaussian),silent=TRUE)
ok("W7 incompatible externally supplied Gaussian fits are rejected early",
  inherits(badExternal,"try-error"))

badDat<-dat;badDat@data[1,badDat@name.response]<-
  badDat@data[1,badDat@name.response]+1
badProvenance<-try(saemixCopula(mod,badDat,ctl,gaussianFit=gaussian),silent=TRUE)
ok("W8 same-size Gaussian fits from a different dataset are rejected",
  inherits(badProvenance,"try-error"))

refitError<-tryCatch(suppressWarnings(saemixCopula(mod,dat,ctl,
  gaussianFit=gaussian,nsamp=3L,max.iter=40L,
  familySet=c("gaussian","gumbel","clayton","frank"),
  populationControl=list(jointBackend="not-a-backend",guard=FALSE))),
  error=function(e)e)
ok("W9 a thrown fresh refit preserves the completed learning audit",
  inherits(refitError,"saemixCopulaRefitError")&&
    identical(refitError$stage,"fixed-refit")&&
    inherits(refitError$learning,"saemixCopulaLearning")&&
    is.data.frame(refitError$learning$ranking)&&
    is.list(refitError$learning$commonQ$validation)&&
    all(c("learning","fixed_fit","copula_workflow")%in%names(refitError$timing)))
ok("W10 refit failure retains the original error message and compact cause",
  inherits(refitError,"error")&&is.character(refitError$message)&&
    identical(refitError$message,refitError$cause$message)&&
    is.character(refitError$cause$class))

badContract<-try(saemixCopula(mod,dat,ctl,
  populationControl=list(familySet="gumbel")),silent=TRUE)
ok("W11 high-level overrides cannot undo the frozen discrete-model contract",
  inherits(badContract,"try-error")&&
    grepl("fixed-model contract",as.character(badContract),fixed=TRUE))

cat(sprintf("\n%d failure(s)\n",nFail));if(nFail>0)quit(status=1L)
