## Direct-parameter finite adaptive proposal/Q layer.
suppressMessages({library(devtools);load_all("C:/package/saemix-copula",quiet=TRUE)})
nFail<-0L
ok<-function(label,pass,extra=""){
 if(!isTRUE(pass))nFail<<-nFail+1L
 cat(sprintf("%-70s %s %s\n",label,if(isTRUE(pass))"PASS" else "**FAIL**",extra))
}
makeV<-function(family,par)rvinecopulib::vinecop_dist(
 list(list(rvinecopulib::bicop_dist(family,parameters=par))),
 rvinecopulib::dvine_structure(1:2))
set.seed(991);N<-260L;S<-3L;gTruth<-1.9;beta<-log(c(10,4))
margins<-list(parameterMarginGamma(9,"mean"),parameterMarginLognormal(.3,"median"))
for(j in 1:2)margins[[j]]$free[]<-FALSE
u<-rvinecopulib::rvinecop(N*S,makeV("gumbel",gTruth))
pred<-matrix(rep(beta,each=N*S),N*S,2)
phi<-log(parameterMarginsQuantile(u,pred,margins))
state<-list(populationScale="parameter",d=2L,margins=margins,
 vine=makeV("gaussian",.45),transform=c(1L,1L),
 poolPhi=phi,poolX=matrix(1,N*S,2),poolW=rep(1/(N*S),N*S),
 poolSubject=rep(seq_len(N),S),locMap=diag(2),betaCurrent=beta,
 betaFitted=beta,betaFree=integer(),cores=1L)
sl<-parameterAdaptiveSlices(state,nsamp=4L,seed=7L)
ok("A1 adaptive slices contain one posterior particle per subject",
   identical(dim(sl$u),c(N,2L,4L)) && all(sl$u>0&sl$u<1))
ps<-parameterAdaptiveVineCandidates(state,nsamp=4L,seed=7L,
 familySet=c("gaussian","gumbel","clayton","frank"),top=3L)
ok("A2 standard rvine selection proposes structure/family candidates",
   inherits(ps,"saemixParameterVineCandidates") && length(ps$candidates)>=1L)
qs<-parameterAdaptiveQStep(state,ps,criterion="bic",includeIncumbent=TRUE,
 maxit=120L,minImprovement=0)
incQ<-qs$ranking$Q[qs$ranking$incumbent][1]
ok("A3 candidate acceptance uses a common absolute-phi particle Q",
   isTRUE(qs$commonMeasure) && qs$winner$conv==0L && qs$winner$value*N>=incQ-1e-6)
pop<-copulaPopulation(state$vine,margins=margins,scale="parameter",
 proposalScale=c(.25,.25),mode="joint")
fixed<-parameterAdaptivePopulation(pop,qs)
ok("A4 learned model becomes a fixed-fingerprint refit specification",
   inherits(fixed,"saemixPopulation") && is.null(fixed$arguments$familySet) &&
   inherits(fixed$arguments$vine,"vinecop_dist"))

candidateMargins<-list(
 list(parameterMarginGamma(5,"mean"),parameterMarginLognormal(.45,"median")),
 list(parameterMarginLognormal(.45,"mean"),parameterMarginLognormal(.45,"median")))
pc<-parameterAdaptiveCandidates(state,marginSets=candidateMargins,nsamp=3L,seed=8L,
 familySet=c("gaussian","gumbel","clayton","frank"),top=2L)
candidateNames<-vapply(pc$candidates,function(z)z$margins[[1]]$name,character(1))
ok("A5 proposal generation retains every declared compatible margin family",
   inherits(pc,"saemixParameterPopulationCandidates") &&
     all(c("gamma","lognormal")%in%candidateNames))
mq<-parameterAdaptiveQStep(state,list(
 list(vine=makeV("gumbel",1.5),margins=candidateMargins[[1]]),
 list(vine=makeV("gumbel",1.5),margins=candidateMargins[[2]])),
 criterion="bic",includeIncumbent=FALSE,maxit=160L,minImprovement=0)
ok("A6 margin-family jumps use the same absolute-phi auxiliary function",
   isTRUE(mq$commonMeasure) && identical(mq$winner$margins[[1]]$name,"gamma"))

## Smoothing tracking: fully solve the first decreasing-gain auxiliary
## function, then take one warm-started contraction step per later iteration.
v1<-rvinecopulib::vinecop_dist(list(),rvinecopulib::dvine_structure(1))
set.seed(992);p1<-matrix(log(rgamma(120,8,scale=10/8)),ncol=1)
copulaSet(v1,margins=list(parameterMarginGamma(5,"mean")),
  populationScale="parameter",proposalScale=.3,jointMaxit=7L,
  jointSmoothMaxit=1L,jointFinalMaxit=35L,guard=FALSE)
.cop$curEta<-p1-log(10);.cop$poolPhi<-p1;.cop$poolEta<-.cop$curEta
.cop$poolX<-matrix(1,nrow(p1),1);.cop$poolW<-rep(1/nrow(p1),nrow(p1))
.cop$locMap<-matrix(1,1,1);.cop$betaCurrent<-log(10);.cop$betaFree<-integer()
.cop$transform<-1L
copulaMstep(10L,final=FALSE,nbiterExploration=10L);explore<-copulaGet()$lastJoint
copulaMstep(11L,final=FALSE,nbiterExploration=10L);transition<-copulaGet()$lastJoint
copulaMstep(12L,final=FALSE,nbiterExploration=10L);tracking<-copulaGet()$lastJoint
ok("A7 smoothing starts with a full solve then warm-started tracking",
   explore$maxit==7L&&!explore$trackingSmoothing&&
     transition$maxit==35L&&transition$transitionSmoothing&&
     tracking$maxit==1L&&tracking$trackingSmoothing)
copulaClear()
cat(sprintf("\n%d failure(s)\n",nFail));if(nFail>0)quit(status=1)
