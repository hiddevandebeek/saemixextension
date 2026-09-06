## End-to-end SAEM smoke test with genuinely non-Gaussian eta margins.
suppressMessages({library(devtools); load_all("C:/package/saemix-copula", quiet=TRUE)})
setwd("C:/package/saemix-copula/copula")
source("R/etaCopula.R"); source("R/simEta.R"); source("R/simpleModels.R")

nFail <- 0L
ok <- function(lbl, pass, extra="") {
  if(!isTRUE(pass)) nFail <<- nFail+1L
  cat(sprintf("%-64s %s %s\n",lbl,if(isTRUE(pass)) "PASS" else "**FAIL**",extra))
}

set.seed(931)
m <- MODELS$iv1; N <- 55L
truthMargins <- list(copulaMarginStudent(.30, 5), copulaMarginLaplace(.30))
truthVine <- makeTwinVines(etaR(2), "gumbel")$alt
eta <- copulaMarginsQuantile(rvinecopulib::rvinecop(N, truthVine), truthMargins)
psi <- sweep(exp(eta), 2, m$true, "*"); colnames(psi) <- m$names
nt <- length(m$times)
dd <- data.frame(id=rep(seq_len(N),each=nt),dose=m$dose,time=rep(m$times,N))
f <- m$f(psi,dd$id,cbind(dd$dose,dd$time)); dd$y <- pmax(f*(1+.08*rnorm(nrow(dd))),1e-6)
dat <- saemixDataFor(dd); mod <- saemixModelFor(m)
ctl <- list(seed=931,save=FALSE,save.graphs=FALSE,print=FALSE,
  displayProgress=FALSE,nbiter.saemix=c(80,50),nbiter.mcmc=c(2,2,2,0),warnings=FALSE)

copulaClear(); stock <- saemix::saemix(mod,dat,ctl)
sd0 <- sqrt(diag(stock@results@omega)); tau0 <- 2/pi*asin(cov2cor(stock@results@omega)[1,2])
v0 <- etaVine(list(copulaBicopFromTau(
  rvinecopulib::bicop_dist("gumbel",parameters=1.2),tau0)),2)
startMargins <- list(copulaMarginStudent(sd0[1], 12),copulaMarginLaplace(sd0[2]))
copulaSet(v0,margins=startMargins,mode="joint",activeFrom=61L,restartBurn=30L,
  refitEvery=5L,jointMaxit=5L,jointFinalMaxit=60L,guardAction="warn")
fit <- saemix::saemix(mod,dat,ctl); st <- copulaGet(); stObject <- copulaGet(fit)
copulaClear()

ok("S1 non-Gaussian margins survive the complete SAEM path",
   identical(vapply(st$margins,`[[`,character(1),"name"),c("student","laplace")))
ok("S1b every native marginal parameter reaches the joint M-step",
   st$lastJoint$nMarginParameters==3L && abs(st$margins[[1]]$parameters["df"]-12)>.2,
   sprintf("student=(%.3f,%.2f) laplaceSD=%.3f",
     st$margins[[1]]$parameters["sd"],st$margins[[1]]$parameters["df"],
     st$margins[[2]]$parameters["sd"]))
ok("S1c final flexible-margin optimizer converges",st$lastJoint$conv==0L,
   sprintf("backend=%s",st$lastJoint$backend))
ok("S1d fitted structural parameters and scales remain finite",
   all(is.finite(fit@results@fixed.effects)) && all(is.finite(st$sd)) && all(st$sd>0))
ok("S1e fitted copula state persists on the returned saemix object",
   identical(copulaFingerprint(stObject$vine, stObject$d, stObject$margins),
             copulaFingerprint(st$vine, st$d, st$margins)))
ok("S1f AIC/BIC parameter count uses margins plus native vine parameters",
   fit@results@npar.est == stock@results@npar.est + 1L,
   sprintf("stock=%d flexible=%d", stock@results@npar.est, fit@results@npar.est))

fit@options$nmc.is <- 200L
scored <- suppressWarnings(llis.saemix(fit, defensive=.25,
                                       batch=50L, seed=77L))
isd <- attr(scored,"saemix.copula.likelihood",exact=TRUE)
ok("S2 copula-aware observed likelihood works after global state is cleared",
   is.finite(scored@results@ll.is) && isd$draws_used==200L &&
     length(isd$ess)==N && all(is.finite(isd$ess)))
tmp <- tempfile(fileext=".rds"); saveRDS(scored,tmp); restored <- readRDS(tmp); unlink(tmp)
rescored <- suppressWarnings(llisCopula.saemix(restored, defensive=.25,
                                                batch=50L, seed=77L))
ok("S2b persisted margin/vine snapshot survives serialization",
   identical(scored@results@ll.is,rescored@results@ll.is))

## Under the full Gaussian-vine null the new scorer must be draw-for-draw
## identical to stock importance sampling when its defensive component is off.
stock@options$nmc.is <- 200L
sdg <- sqrt(diag(stock@results@omega))
gst <- list(vine=etaVineGaussian(cov2cor(stock@results@omega)),
            margins=lapply(sdg,copulaMarginNormal),sd=sdg,d=2L)
gcopy <- stock
attr(gcopy,"saemix.copula") <- copulaSnapshot(gst,
  etaIndex=stock@model@indx.omega,
  variableName=stock@model@name.modpar[stock@model@indx.omega])
legacyLL <- withSeed(78L, llis.saemix(stock))@results@ll.is
copulaLL <- suppressWarnings(llisCopula.saemix(gcopy, defensive=0,
                                                batch=100L, seed=78L))@results@ll.is
ok("S3 Gaussian-vine likelihood scorer equals stock IS draw for draw",
   abs(legacyLL-copulaLL)<1e-8,
   sprintf("difference=%.3g",copulaLL-legacyLL))

cat(sprintf("\n%d failure(s)\n",nFail)); if(nFail>0) quit(status=1)
