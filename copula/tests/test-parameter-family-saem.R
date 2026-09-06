## Every built-in direct continuous margin through the actual saemix path.
suppressMessages({library(devtools);load_all("C:/package/saemix-copula",quiet=TRUE)})
nFail<-0L
ok<-function(label,pass,extra=""){
  if(!isTRUE(pass))nFail<<-nFail+1L
  cat(sprintf("%-68s %s %s\n",label,if(isTRUE(pass))"PASS" else "**FAIL**",extra))
}
v1<-rvinecopulib::vinecop_dist(list(),rvinecopulib::dvine_structure(1))
specs<-list(
  normal=list(tr=0L,anchor=2,truth=parameterMarginNormal(.55),
    start=parameterMarginNormal(.8),noise=.25,anchorTol=.25,nativeTol=.30),
  student=list(tr=0L,anchor=2,truth=parameterMarginStudent(.55,6),
    start=parameterMarginStudent(.8,12),noise=.25,anchorTol=.30,nativeTol=8),
  laplace=list(tr=0L,anchor=2,truth=parameterMarginLaplace(.55),
    start=parameterMarginLaplace(.8),noise=.25,anchorTol=.30,nativeTol=.35),
  lognormal=list(tr=1L,anchor=2,truth=parameterMarginLognormal(.28,"median"),
    start=parameterMarginLognormal(.45,"median"),noise=.20,anchorTol=.20,nativeTol=.18),
  gamma=list(tr=1L,anchor=2,truth=parameterMarginGamma(8,"mean"),
    start=parameterMarginGamma(4,"mean"),noise=.20,anchorTol=.20,nativeTol=4),
  weibull=list(tr=1L,anchor=2,truth=parameterMarginWeibull(3,"mean"),
    start=parameterMarginWeibull(1.8,"mean"),noise=.20,anchorTol=.22,nativeTol=1.4),
  beta=list(tr=3L,anchor=.4,truth=parameterMarginBeta(20),
    start=parameterMarginBeta(10),noise=.04,anchorTol=.06,nativeTol=12))

fit_one<-function(z,index){
  set.seed(3300L+index);N<-80L
  predictor<-matrix(transpsi(matrix(z$anchor,nrow=1),z$tr)[1],N,1)
  psiI<-parameterMarginsQuantile(matrix(runif(N),ncol=1),predictor,list(z$truth))[,1]
  dd<-data.frame(id=rep(seq_len(N),each=4L),time=rep(seq_len(4L),N))
  dd$y<-psiI[dd$id]+stats::rnorm(nrow(dd),sd=z$noise)
  f<-function(psi,id,xidep)psi[id,1]
  dat<-saemix::saemixData(name.data=dd,header=TRUE,name.group="id",
    name.predictors="time",name.response="y",verbose=FALSE)
  mod<-saemix::saemixModel(model=f,modeltype="structural",description="",
    psi0=matrix(z$anchor,nrow=1,dimnames=list(NULL,"level")),transform.par=z$tr,
    covariance.model=matrix(1,1,1),omega.init=matrix(.35^2,1,1),
    error.model="constant",verbose=FALSE)
  ctl<-list(seed=3300L+index,save=FALSE,save.graphs=FALSE,print=FALSE,
    displayProgress=FALSE,nbiter.saemix=c(80,60),nbiter.mcmc=c(2,2,1,0),
    warnings=FALSE,ll.is=FALSE,fim=FALSE)
  pop<-copulaPopulation(v1,margins=list(z$start),scale="parameter",
    proposalScale=.35,mode="joint",refitEvery=1L,jointMaxit=20L,
    jointFinalMaxit=120L,guard=FALSE)
  fit<-saemix::saemix(mod,dat,ctl,population=pop);state<-copulaGet(fit)
  anchorError<-abs(as.numeric(fit@results@fixed.psi)-z$anchor)
  nativeError<-max(abs(state$margins[[1]]$parameters-z$truth$parameters))
  list(pass=state$lastJoint$final&&state$lastJoint$conv==0L&&
      anchorError<z$anchorTol&&nativeError<z$nativeTol,
    anchor=as.numeric(fit@results@fixed.psi),native=state$margins[[1]]$parameters,
    anchorError=anchorError,nativeError=nativeError)
}
for(i in seq_along(specs)){
  ans<-fit_one(specs[[i]],i)
  ok(paste0("F",i," ",names(specs)[i]," margin fits through direct SAEM-MCMC"),
    ans$pass,sprintf("anchor=%.3f native=%s",ans$anchor,
      paste(sprintf("%.3f",ans$native),collapse=",")))
}
cat(sprintf("\n%d failure(s)\n",nFail));if(nFail>0)quit(status=1)
