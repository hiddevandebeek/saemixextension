## Exact Normal/lognormal score collapse into the mutable C++ vine context.
suppressPackageStartupMessages({library(devtools); load_all("C:/package/saemix-copula", quiet=TRUE)})
n_fail <- 0L
ok <- function(label, pass, detail="") {
  if(!isTRUE(pass)) n_fail <<- n_fail + 1L
  cat(sprintf("%-70s %s %s\n", label, if(isTRUE(pass))"PASS" else "**FAIL**", detail))
}

set.seed(82501)
d <- 4L; n <- 240L
vine <- rvinecopulib::vinecop_dist(list(
  list(rvinecopulib::bicop_dist("frank", parameters=2.2),
       rvinecopulib::bicop_dist("clayton", rotation=180, parameters=1.6),
       rvinecopulib::bicop_dist("gaussian", parameters=.55)),
  list(rvinecopulib::bicop_dist("joe", rotation=270, parameters=1.4),
       rvinecopulib::bicop_dist("frank", parameters=-2.1)),
  list(rvinecopulib::bicop_dist("gaussian", parameters=.25))),
  rvinecopulib::rvine_structure_sim(d))
margins <- list(copulaMarginNormal(.28), copulaMarginNormal(.36),
  copulaMarginCovariateLognormal(log(72), .20),
  copulaMarginCovariateLognormal(log(91), .25))
u <- rvinecopulib::rvinecop(n, vine)
joint <- copulaMarginsQuantile(u, margins)
X <- matrix(1, n, 2L); locMap <- cbind(diag(2L), matrix(0,2L,2L))
beta <- c(.12,-.08); joint[,1:2] <- sweep(joint[,1:2],2,beta,"+")
w <- rexp(n); w <- w/sum(w)

transformed <- joint; transformed[,3:4] <- log(transformed[,3:4])
mu <- c(beta, margins[[3]]$parameters[["meanlog"]],
             margins[[4]]$parameters[["meanlog"]])
z <- sweep(transformed, 2, mu, "-")
normalMargins <- lapply(c(.28,.36,.20,.25), copulaMarginNormal)
qOriginal <- sum(w*copulaLogPrior(joint-copulaLocation(X,beta,locMap),
  vine,margins=margins,numericalPolicy="exact"))
qNormal <- sum(w*copulaLogPrior(z,vine,margins=normalMargins,
  numericalPolicy="exact"))
jacobian <- sum(w*rowSums(log(joint[,3:4,drop=FALSE])))
ok("GM1 transformed Normal Q differs only by the lognormal Jacobian",
   abs(qOriginal-(qNormal-jacobian))<2e-10,
   sprintf("error=%.3g",abs(qOriginal-(qNormal-jacobian))))

startMargins <- margins
startMargins[[1]]$parameters[] <- .34
startMargins[[2]]$parameters[] <- .31
startMargins[[3]]$parameters[] <- c(log(69),.24)
startMargins[[4]]$parameters[] <- c(log(96),.21)
fit <- copulaMaximiseJointMarginsGaussianized(joint,w,startMargins,vine,d,
  maxit=35L,X=X,locMap=locMap,beta0=c(.08,-.03),betaFree=1:2,
  cores=1L,diagnostics=TRUE)
qStart <- sum(w*copulaLogPrior(joint-copulaLocation(X,c(.08,-.03),locMap),
  vine,margins=startMargins,numericalPolicy="exact"))
qFit <- sum(w*copulaLogPrior(joint-copulaLocation(X,fit$beta,locMap),fit$vine,
  margins=fit$margins,numericalPolicy="exact"))
ok("GM2 collapsed C++ M-step increases the literal arbitrary-margin Q",
   is.list(fit)&&qFit+1e-10>=qStart,
   sprintf("Q %.6f -> %.6f",qStart,qFit))
ok("GM3 collapsed backend is reported explicitly",
   identical(fit$backend,"cpp-gaussianized-margins-pure"), fit$backend)

oldState<-as.list(.cop);on.exit(copulaRestoreState(oldState),add=TRUE)
copulaSet(vine,margins=margins,conditioning=list(values=joint[1:20,3:4],
  variableName=c("WT","eGFR")),proposalCorrelation="gaussian-score",
  warmStartOnActivate=FALSE)
proposalR<-copulaImpliedCor()
expectedR<-copulaGaussianRvineCor(copulaAsGaussian(vine,d),d)
ok("GM4 proposal correlation collapses analytically without changing the vine",
   max(abs(proposalR-expectedR))<1e-12 && identical(.cop$vine,vine),
   sprintf("maxabs=%.3g",max(abs(proposalR-expectedR))))

repeated<-joint[rep(1:20,each=8),3:4,drop=FALSE]
compressedEval<-copulaMarginsEvaluate(repeated,margins[3:4],
  numericalPolicy="exact")
preindexed<-copulaMarginsCompression(repeated,1:2)
preindexedEval<-copulaMarginsEvaluate(repeated,margins[3:4],
  numericalPolicy="exact",compression=preindexed)
literalU<-cbind(margins[[3]]$cdf(repeated[,1],margins[[3]]$parameters),
                margins[[4]]$cdf(repeated[,2],margins[[4]]$parameters))
literalLog<-margins[[3]]$log_density(repeated[,1],margins[[3]]$parameters)+
  margins[[4]]$log_density(repeated[,2],margins[[4]]$parameters)
ok("GM5 repeated covariate margin evaluation is an exact compression",
   identical(compressedEval$u,literalU)&&
     identical(compressedEval$log_margin,as.numeric(literalLog))&&
     identical(preindexedEval,compressedEval),"")

changedMargins<-margins[3:4]
changedMargins[[1]]<-copulaMarginWithParameters(changedMargins[[1]],c(log(75),.23))
changedMargins[[2]]<-copulaMarginWithParameters(changedMargins[[2]],c(log(88),.19))
changedCached<-copulaMarginsEvaluate(repeated,changedMargins,
  numericalPolicy="exact",compression=preindexed)
changedLiteral<-copulaMarginsEvaluate(repeated,changedMargins,
  numericalPolicy="exact",compression=vector("list",2L))
ok("GM6 pre-indexing remains exact while margin parameters change",
   identical(changedCached,changedLiteral),"")

cppWeighted<-copula_weighted_log_pdf(vine,u,w,threads=2L)
referenceWeighted<-sum(w*log(rvinecopulib::dvinecop(u,vine)))
ok("GM7 direct pair-log recursion equals the library R-vine density",
   abs(cppWeighted-referenceWeighted)<2e-12,
   sprintf("error=%.3g",abs(cppWeighted-referenceWeighted)))

repIndex<-rep(1:40,each=4L)
repJoint<-joint[repIndex,,drop=FALSE]
repJoint[,1:2]<-repJoint[,1:2]+matrix(rnorm(nrow(repJoint)*2,sd=.01),ncol=2)
repTransformed<-repJoint;repTransformed[,3:4]<-log(repTransformed[,3:4])
repX<-matrix(1,nrow(repJoint),2L)
repXaug<-cbind(repX,1,1)
repLaug<-rbind(locMap,c(0,0,1,0),c(0,0,0,1))
repBeta<-c(beta,margins[[3]]$parameters[["meanlog"]],
                margins[[4]]$parameters[["meanlog"]])
repSd<-c(margins[[1]]$parameters[["sd"]],margins[[2]]$parameters[["sd"]],
         margins[[3]]$parameters[["sdlog"]],margins[[4]]$parameters[["sdlog"]])
repLayout<-copulaParameterLayout(copulaPadFlat(vine,d))
repPar<-c(repBeta,log(repSd),repLayout$par)
repW<-rep(1/nrow(repJoint),nrow(repJoint))
repCtx<-copula_q_context_create(vine,repTransformed,repXaug,repLaug,
  repBeta,1:4,repW)
repCpp<-copula_q_context_eval_const(repCtx,repPar,2L)
repResidual<-repTransformed-copulaLocation(repXaug,repBeta,repLaug)
repLiteral<--sum(repW*copulaLogPrior(repResidual,vine,
  margins=lapply(repSd,copulaMarginNormal),numericalPolicy="exact"))
ok("GM8 compiled repeated-row algebra preserves the literal Q",
   abs(repCpp-repLiteral)<2e-11,
   sprintf("error=%.3g",abs(repCpp-repLiteral)))

rawMarginLayout<-copulaMarginLayout(margins)
rawPar<-c(beta,rawMarginLayout$par,repLayout$par)
rawLower<-c(rep(-Inf,2),rawMarginLayout$lower,repLayout$lower)
rawUpper<-c(rep( Inf,2),rawMarginLayout$upper,repLayout$upper)
rawToCpp<-c(1L,2L,5L,6L,3L,7L,4L,8L,9:14)
logTransform<-seq_along(rawPar)%in%c(3L,4L,6L,8L)
toCpp<-function(raw) {
  ans<-numeric(length(raw));ans[rawToCpp]<-raw
  ans[5:8]<-log(ans[5:8]);ans
}
graphGradient<-copula_q_context_gradient_mapped(repCtx,rawPar,toCpp(rawPar),
  rawLower,rawUpper,rawToCpp,logTransform,threads=2L,step_size=1e-3)
literalGradient<-vapply(seq_along(rawPar),function(j){
  lo<-hi<-rawPar;lo[j]<-lo[j]-1e-3;hi[j]<-hi[j]+1e-3
  (copula_q_context_eval_const(repCtx,toCpp(hi),1L)-
   copula_q_context_eval_const(repCtx,toCpp(lo),1L))/2e-3
},numeric(1))
ok("GM9 dependency-cached gradient equals full central differences",
   max(abs(graphGradient-literalGradient))<2e-10,
   sprintf("max error=%.3g",max(abs(graphGradient-literalGradient))))

populationResidual<-joint-copulaLocation(X,beta,locMap)
populationLocation<-c(0,0,margins[[3]]$parameters[["meanlog"]],
                           margins[[4]]$parameters[["meanlog"]])
populationScale<-c(margins[[1]]$parameters[["sd"]],
  margins[[2]]$parameters[["sd"]],margins[[3]]$parameters[["sdlog"]],
  margins[[4]]$parameters[["sdlog"]])
cppEnergy<-copula_ueta_gaussianized(vine,populationResidual,
  populationLocation,populationScale,c(FALSE,FALSE,TRUE,TRUE),threads=2L)
literalEnergy<--copulaLogPrior(populationResidual,vine,margins=margins,
  numericalPolicy="exact")
ok("GM10 fused MCMC energy equals the arbitrary-margin R likelihood",
   max(abs(cppEnergy-literalEnergy))<2e-12,
   sprintf("max error=%.3g",max(abs(cppEnergy-literalEnergy))))

cat(sprintf("\n%d failure(s)\n",n_fail)); if(n_fail)quit(status=1L)
