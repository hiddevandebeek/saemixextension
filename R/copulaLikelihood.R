## Copula-aware observed likelihood by defensive importance sampling.
## The public llis.saemix() dispatches here whenever the fitted SaemixObject
## carries an immutable copula snapshot.  No mutable global .cop state is read.

copulaLogAddExp <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  out[is.infinite(m) & m < 0] <- -Inf
  out
}

copulaRowLogSumExp <- function(x) {
  x <- as.matrix(x)
  m <- apply(x, 1L, max)
  out <- m + log(rowSums(exp(x - m)))
  out[is.infinite(m) & m < 0] <- -Inf
  out
}

copulaResponseLogLikBatch <- function(saemixObject, phi, B) {
  dat <- saemixObject["data"]
  model <- saemixObject["model"]
  res <- saemixObject["results"]
  N <- dat["N"]; ntot <- dat["ntot.obs"]
  yobs <- dat["data"][, dat["name.response"]]
  xnames <- c(dat["name.predictors"], dat["name.cens"],
              dat["name.mdv"], dat["name.ytype"])
  xind <- dat["data"][, xnames, drop=FALSE]
  IdB <- rep(0:(B - 1L), each=ntot) * N +
    rep(dat["data"][, "index"], B)
  XB <- do.call(rbind, rep(list(xind), B))
  yB <- rep(yobs, B)
  psi <- transphi(phi, model["transform.par"])
  f <- model["model"](psi, IdB, XB)
  idxExp <- which(model["error.model"] == "exponential")
  if (length(idxExp)) for (j in idxExp) {
    take <- XB$ytype == j
    f[take] <- log(cutoff(f[take]))
  }
  if (model["modeltype"] == "structural") {
    g <- error(f, res["respar"], XB$ytype)
    obs <- -0.5 * ((yB - f) / g)^2 - log(g) - 0.5 * log(2*pi)
    obs[!is.finite(f) | !is.finite(g) | g <= 0 | !is.finite(obs)] <- -Inf
  } else {
    obs <- as.numeric(f)
    obs[!is.finite(obs)] <- -Inf
  }
  io <- matrix(0, nrow=N, ncol=max(dat["nind.obs"]))
  for (i in seq_len(N)) io[i, seq_len(dat["nind.obs"][i])] <- 1
  ioB <- matrix(rep(t(io), B), ncol=ncol(io), byrow=TRUE)
  DYF <- matrix(0, nrow=ncol(io), ncol=N * B)
  DYF[which(t(ioB) != 0)] <- obs
  matrix(colSums(DYF), nrow=N, ncol=B)
}

llisCopula.saemix <- function(saemixObject, defensive=0.25, batch=100L,
                              seed=NULL) {
  state <- copulaGet(saemixObject)
  if (!inherits(state, "saemixCopulaSnapshot"))
    stop("fitted object has no immutable copula snapshot")
  if (!identical(copulaFingerprint(state$vine, state$d, state$margins),
                 state$fingerprint))
    stop("stored copula fingerprint is inconsistent")
  if (any(vapply(state$margins[seq_len(state$dEta %||% state$d)],
      `[[`, character(1), "type") != "continuous"))
    stop("copula importance sampling currently requires continuous eta margins")
  model <- saemixObject["model"]; dat <- saemixObject["data"]
  res <- saemixObject["results"]
  etaIndex <- as.integer(model["indx.omega"])
  dEta <- as.integer(state$dEta %||% state$d)
  dConditioning <- as.integer(state$dConditioning %||% 0L)
  hasConditioning <- dConditioning > 0L
  if (!identical(etaIndex, as.integer(state$etaIndex)) ||
      length(etaIndex) != dEta || state$d != dEta + dConditioning)
    stop("stored copula/eta ordering does not match the fitted model")
  if (hasConditioning && (is.null(state$conditioning) ||
      ncol(as.matrix(state$conditioning)) != dConditioning ||
      nrow(as.matrix(state$conditioning)) != dat["N"]))
    stop("stored conditioning values are incompatible with the fitted model")
  exactPosterior <- copulaLinearGaussianPosterior(saemixObject, state,
    etaIndex, res["mean.phi"][, etaIndex, drop = FALSE])
  if (!is.null(exactPosterior)) {
    conditioningLogLik <- if (hasConditioning &&
        identical(state$likelihoodTarget %||% "joint", "joint"))
      copulaGaussianFremConditioningLogDensity(as.matrix(state$conditioning),
        state$vine, state$margins, dEta) else rep(0, dat["N"])
    perSubject <- exactPosterior$responseLogLikelihood + conditioningLogLik
    observed <- sum(perSubject)
    saemixObject["results"]["LL"] <- observed
    saemixObject["results"]["ll.is"] <- observed
    np <- saemixObject["results"]["npar.est"]
    saemixObject["results"]["aic.is"] <- -2 * observed + 2 * np
    saemixObject["results"]["bic.is"] <- -2 * observed + log(dat["N"]) * np
    saemixObject["results"]["bic.covariate.is"] <- -2 * observed +
      log(dat["N"]) * saemixObject["results"]["nbeta.random"] +
      log(sum(dat["nind.obs"])) * saemixObject["results"]["nbeta.fixed"]
    attr(saemixObject, "saemix.copula.likelihood") <- list(
      method = "exact-linear-gaussian", exact = TRUE,
      per_subject_loglik = perSubject, se_loglik_subject = rep(0, dat["N"]),
      se_loglik_total = 0, se_deviance = 0,
      ess = rep(Inf, dat["N"]), ess_fraction = rep(1, dat["N"]),
      ess_min = Inf, ess_mean = Inf, max_weight_fraction = rep(0, dat["N"]),
      invalid_draws = integer(dat["N"]), pit_clipped = 0L,
      density_floored = 0L, transform_invalid = 0L,
      likelihood_target = state$likelihoodTarget %||% "joint",
      state_fingerprint = state$fingerprint)
    return(saemixObject)
  }
  M <- as.integer(saemixObject["options"]$nmc.is)
  if (length(M) != 1L || is.na(M) || M < 1L) stop("nmc.is must be positive")
  batch <- max(1L, as.integer(batch)); defensive <- as.numeric(defensive)
  if (!is.finite(defensive) || defensive < 0 || defensive >= 1)
    stop("defensive must be in [0,1)")
  if (hasConditioning && defensive > 0) {
    warning(paste0("defensive prior-mixture sampling is not yet available for ",
      "joint parameter-covariate vines; using the exact fitted density in the ",
      "importance weights with defensive=0"), call. = FALSE)
    defensive <- 0
  }

  run <- function() {
    N <- dat["N"]; d <- dEta; nu <- saemixObject["options"]$nu.is
    meanPhi <- res["mean.phi"][, etaIndex, drop=FALSE]
    condMean <- res["cond.mean.phi"][, etaIndex, drop=FALSE]
    condVar <- res["cond.var.phi"][, etaIndex, drop=FALSE]
    condVar[!is.finite(condVar) | condVar <= 1e-12] <- 1e-12
    condSd <- sqrt(condVar)
    logSumW <- logSumW2 <- rep(-Inf, N)
    maxLogW <- rep(-Inf, N); invalid <- integer(N)
    pitClipped<-0L;densityFloored<-0L;transformInvalid<-0L
    LL <- numeric(); cumulative <- integer(); used <- 0L; priorBlocks <- 0L; k <- 0L
    logConst <- 0
    idxExp <- which(model["error.model"] == "exponential")
    if (length(idxExp)) {
      yobs <- dat["data"][, dat["name.response"]]
      logConst <- -sum(yobs[dat["data"][, "ytype"] %in% idxExp])
    }
    while (used < M) {
      B <- min(batch, M - used); k <- k + 1L
      nPrior <- round(B * defensive); eps <- nPrior / B
      priorBlocks <- priorBlocks + nPrior
      nT <- B - nPrior
      subjectRow <- rep(seq_len(N), B)
      pop <- meanPhi[subjectRow, , drop=FALSE]
      centre <- condMean[subjectRow, , drop=FALSE]
      sdev <- condSd[subjectRow, , drop=FALSE]
      z <- trnd.mlx(nu, N * B, d)
      z <- matrix(z, nrow=N * B, ncol=d)
      phiEta <- centre + sdev * z
      if (nPrior) {
        take <- rep(c(rep(FALSE, nT), rep(TRUE, nPrior)), each=N)
        uprior <- rvinecopulib::rvinecop(sum(take), state$vine)
        phiEta[take, ] <- pop[take, , drop=FALSE] +
          copulaMarginsQuantile(uprior, state$margins)
      }
      eta <- phiEta - pop
      logPrior <- if (hasConditioning) {
        conditioningRows <- as.matrix(state$conditioning)[subjectRow, , drop=FALSE]
        jointState <- cbind(eta, conditioningRows)
        copulaGaussianFremLogPrior(
          jointState, state$vine, state$margins, dEta,
          likelihoodTarget = state$likelihoodTarget)
      } else {
        copulaLogPrior(eta, state$vine, margins=state$margins)
      }
      zEval <- (phiEta - centre) / sdev
      logQT <- rowSums(log(tpdf.mlx(zEval, nu))) - rowSums(log(sdev))
      logQ <- if (eps == 0) logQT else if (eps == 1) logPrior else
        copulaLogAddExp(log1p(-eps) + logQT, log(eps) + logPrior)
      phi <- res["cond.mean.phi"][subjectRow, , drop=FALSE]
      phi[, etaIndex] <- phiEta
      logObs <- copulaResponseLogLikBatch(saemixObject, phi, B)
      invalid <- invalid + rowSums(!is.finite(logObs))
      logW <- logObs + matrix(logPrior - logQ, nrow=N, ncol=B)
      ## A heavy-tailed proposal may produce a few non-finite structural-model
      ## predictions. Such draws have zero importance weight; they must not
      ## poison the row maximum and discard every other valid draw for that
      ## subject. The subject is rejected below only when its entire cumulative
      ## sample has zero finite weight.
      logW[!is.finite(logW)] <- -Inf
      batchSum <- copulaRowLogSumExp(logW)
      batchSum2 <- copulaRowLogSumExp(2 * logW)
      logSumW <- copulaLogAddExp(logSumW, batchSum)
      logSumW2 <- copulaLogAddExp(logSumW2, batchSum2)
      maxLogW <- pmax(maxLogW, apply(logW, 1L, max))
      used <- used + B
      cumulative[k] <- used
      if (any(!is.finite(logSumW)))
        stop("all importance draws are invalid for subject(s): ",
             paste(which(!is.finite(logSumW)), collapse=", "))
      LL[k] <- sum(logSumW - log(used)) + logConst
    }
    logLi <- logSumW - log(M)
    ess <- exp(2 * logSumW - logSumW2)
    seSubject <- sqrt(pmax(1 / ess - 1 / M, 0))
    diag <- list(engine="copula-defensive-is", draws_requested=M,
      draws_used=M, batch_size=batch, cumulative_draws=cumulative,
      defensive_requested=defensive, defensive=priorBlocks/M, nu=nu,
      per_subject_loglik=logLi,
      ess=ess, ess_fraction=ess/M, ess_min=min(ess), ess_mean=mean(ess),
      max_weight_fraction=exp(maxLogW-logSumW),
      se_loglik_subject=seSubject,
      se_loglik_total=sqrt(sum(seSubject^2)),
      se_deviance=2*sqrt(sum(seSubject^2)), invalid_draws=invalid,
      pit_clipped=pitClipped,density_floored=densityFloored,
      transform_invalid=transformInvalid,
      likelihood_target=state$likelihoodTarget %||% "joint",
      state_fingerprint=state$fingerprint)
    saemixObject["results"]["LL"] <- LL
    saemixObject["results"]["ll.is"] <- tail(LL, 1L)
    np <- saemixObject["results"]["npar.est"]
    saemixObject["results"]["aic.is"] <- -2*tail(LL,1L) + 2*np
    saemixObject["results"]["bic.is"] <- -2*tail(LL,1L) + log(N)*np
    saemixObject["results"]["bic.covariate.is"] <-
      -2*tail(LL,1L) + log(N)*saemixObject["results"]["nbeta.random"] +
      log(sum(dat["nind.obs"]))*saemixObject["results"]["nbeta.fixed"]
    attr(saemixObject, "saemix.copula.likelihood") <- diag
    if (any(diag$ess_fraction < .01) || diag$se_deviance > 1 ||
        any(diag$max_weight_fraction > .1))
      warning("copula likelihood importance sampling has low precision; inspect ",
              "attr(fit, 'saemix.copula.likelihood')", call.=FALSE)
    saemixObject
  }
  if (is.null(seed)) run() else withSeed(seed, run())
}

copulaSubjectResponseLogLik <- function(object,i,phiRow) {
  dat<-object["data"];model<-object["model"];res<-object["results"]
  take<-dat["data"][,"index"]==i
  xnames<-c(dat["name.predictors"],dat["name.cens"],dat["name.mdv"],dat["name.ytype"])
  xi<-dat["data"][take,xnames,drop=FALSE]
  yi<-dat["data"][take,dat["name.response"]]
  psi<-transphi(matrix(phiRow,nrow=1),model["transform.par"])
  val<-model["model"](psi,rep(1L,sum(take)),xi)
  if(model["modeltype"]=="likelihood")return(sum(val))
  idxExp<-which(model["error.model"]=="exponential")
  if(length(idxExp))for(j in idxExp) {
    z<-xi$ytype==j;val[z]<-log(cutoff(val[z]))
  }
  g<-error(val,res["respar"],xi$ytype)
  if(any(!is.finite(val))||any(!is.finite(g))||any(g<=0))return(-Inf)
  sum(-.5*((yi-val)/g)^2-log(g)-.5*log(2*pi))
}

mapCopula.saemix <- function(object) {
  state<-copulaGet(object);idx<-as.integer(state$etaIndex)
  if(!identical(idx,as.integer(object["model"]["indx.omega"])))
    stop("stored population ordering does not match the fitted model")
  phiMap<-object["results"]["phi"]
  predictor<-object["results"]["mean.phi"][,idx,drop=FALSE]
  transform<-object["model"]["transform.par"][idx]
  hasConditioning <- (state$dConditioning %||% 0L) > 0L
  for(i in seq_len(object["data"]["N"])) {
    phii<-phiMap[i,]
    objective<-function(p) {
      row<-phii;row[idx]<-p
      llY<-copulaSubjectResponseLogLik(object,i,row)
      eta<-matrix(p-predictor[i,],nrow=1)
      priorState <- if (hasConditioning)
        cbind(eta, as.matrix(state$conditioning)[i, , drop=FALSE]) else eta
      lp<-copulaGaussianFremLogPrior(priorState, state$vine,
        state$margins, state$dEta, "joint")
      out<--(llY+lp)
      if(is.finite(out))out else 1e100
    }
    op<-suppressWarnings(try(stats::optim(phii[idx],objective,method="BFGS",
      control=list(maxit=200L)),silent=TRUE))
    if(!inherits(op,"try-error")&&is.finite(op$value))phiMap[i,idx]<-op$par
  }
  mapPsi<-data.frame(transphi(phiMap,object["model"]["transform.par"]))
  colnames(mapPsi)<-object["model"]["name.modpar"]
  object["results"]["map.psi"]<-mapPsi
  object["results"]["map.phi"]<-data.frame(phiMap)
  eta<-phiMap-object["results"]["mean.phi"]
  colnames(eta)<-paste0("eta.",object["model"]["name.modpar"])
  object["results"]["map.eta"]<-eta
  shrink<-rep(NA_real_,ncol(phiMap));names(shrink)<-paste0("Sh.",object["model"]["name.modpar"],".%")
  object["results"]["map.shrinkage"]<-shrink
  object
}

## Exact individual posterior for the important affine-response submodel.
## Gaussian-copula conditioning gives a Normal prior for eta | c_obs; an
## affine structural model with constant Gaussian error then gives the usual
## conjugate Normal posterior.  Numerical probes only detect applicability;
## they do not approximate the posterior used after the identity is verified.
copulaLinearGaussianPosterior <- function(object, state, idx, predictor) {
  model <- object["model"]; dat <- object["data"]; res <- object["results"]
  dEta <- length(idx); N <- dat["N"]
  if (!identical(model["modeltype"], "structural") ||
      !all(model["error.model"] == "constant") ||
      !identical(state$populationScale, "transformed-additive") ||
      !copulaIsFullGaussianVine(state$vine, state$d) ||
      any(vapply(state$margins[dEta + seq_len(state$d - dEta)], function(m)
        identical(m$type, "discrete"), logical(1))) ||
      any(model["transform.par"][idx] != 0) ||
      !all(vapply(state$margins[seq_len(dEta)], function(m)
        identical(m$name, "normal") && isTRUE(m$centered), logical(1))))
    return(NULL)
  ytypes <- sort(unique(dat["data"][, "ytype"]))
  if (any(2L * ytypes > length(res["respar"])) ||
      any(abs(res["respar"][2L * ytypes]) > 1e-12)) return(NULL)
  conditioning <- if ((state$dConditioning %||% 0L) > 0L)
    as.matrix(state$conditioning) else matrix(numeric(), N, 0L)
  conditional <- try(copulaGaussianFremConditional(
    conditioning, state$vine, state$margins, dEta), silent = TRUE)
  if (inherits(conditional, "try-error")) return(NULL)
  scale <- diag(vapply(state$margins[seq_len(dEta)], function(m)
    m$scale(m$parameters), numeric(1L)), dEta)
  posteriorMean <- matrix(NA_real_, N, dEta)
  posteriorCovariance <- vector("list", N)
  responseLogLikelihood <- rep(NA_real_, N)
  evaluate <- function(phiRow, xi) {
    psi <- transphi(matrix(phiRow, nrow = 1L), model["transform.par"])
    as.numeric(model["model"](psi, rep(1L, nrow(xi)), xi))
  }
  for (i in seq_len(N)) {
    take <- dat["data"][, "index"] == i
    xnames <- c(dat["name.predictors"], dat["name.cens"],
      dat["name.mdv"], dat["name.ytype"])
    xi <- dat["data"][take, xnames, drop = FALSE]
    yi <- dat["data"][take, dat["name.response"]]
    priorMean <- as.numeric(predictor[i, ] +
      scale %*% conditional$mean[i, ])
    priorCov <- scale %*% conditional$covariance[[i]] %*% scale
    baseRow <- object["results"]["mean.phi"][i, ]
    baseRow[idx] <- priorMean
    f0 <- try(evaluate(baseRow, xi), silent = TRUE)
    if (inherits(f0, "try-error") || any(!is.finite(f0))) return(NULL)
    step <- pmax(sqrt(diag(priorCov)), 1e-4)
    J <- matrix(NA_real_, length(f0), dEta)
    for (j in seq_len(dEta)) {
      plus <- minus <- baseRow
      plus[idx[j]] <- plus[idx[j]] + step[j]
      minus[idx[j]] <- minus[idx[j]] - step[j]
      fp <- try(evaluate(plus, xi), silent = TRUE)
      fm <- try(evaluate(minus, xi), silent = TRUE)
      if (inherits(fp, "try-error") || inherits(fm, "try-error")) return(NULL)
      J[, j] <- (fp - fm) / (2 * step[j])
    }
    ## Detect interactions/curvature with two multivariate probes.
    for (coefficient in list(seq(.23, .41, length.out = dEta),
                             seq(-.37, .19, length.out = dEta))) {
      probe <- baseRow; displacement <- step * coefficient
      probe[idx] <- probe[idx] + displacement
      actual <- try(evaluate(probe, xi), silent = TRUE)
      expected <- f0 + as.numeric(J %*% displacement)
      if (inherits(actual, "try-error") || any(!is.finite(actual)) ||
          max(abs(actual - expected)) > 1e-8 * (1 + max(abs(expected))))
        return(NULL)
    }
    g <- error(f0, res["respar"], xi$ytype)
    if (any(!is.finite(g)) || any(g <= 0)) return(NULL)
    weightedJ <- J / g
    postCov <- try(solve(solve(priorCov) + crossprod(weightedJ)), silent = TRUE)
    if (inherits(postCov, "try-error")) return(NULL)
    postMean <- priorMean + postCov %*%
      crossprod(J, (yi - f0) / g^2)
    posteriorMean[i, ] <- postMean
    posteriorCovariance[[i]] <- (postCov + t(postCov)) / 2
    responseCov <- J %*% priorCov %*% t(J) + diag(g^2, length(g))
    responseLogLikelihood[i] <- copulaGaussianLogDensity(
      matrix(yi - f0, nrow = 1L), responseCov)
  }
  list(mean = posteriorMean, covariance = posteriorCovariance,
    responseLogLikelihood = responseLogLikelihood)
}

conddistCopula.saemix <- function(object,nsamp=1,max.iter=NULL,plot=FALSE,...) {
  state<-copulaGet(object);idx<-as.integer(state$etaIndex);N<-object["data"]["N"]
  nsamp<-max(1L,as.integer(nsamp));burn<-if(is.null(max.iter))
    max(100L,sum(object["options"]$nbiter.saemix)) else max(50L,as.integer(max.iter))
  transform<-object["model"]["transform.par"][idx]
  predictor<-object["results"]["mean.phi"][,idx,drop=FALSE]
  start<-object["results"]["phi"]
  exactPosterior <- copulaLinearGaussianPosterior(object,state,idx,predictor)
  if (!is.null(exactPosterior)) {
    phiSamples <- array(NA_real_, c(N, ncol(start), nsamp))
    phiMean <- start; phiVar <- matrix(0, N, ncol(start))
    phiMean[, idx] <- exactPosterior$mean
    for (i in seq_len(N)) {
      phiVar[i, idx] <- diag(exactPosterior$covariance[[i]])
      draws <- matrix(rnorm(nsamp * length(idx)), ncol = length(idx)) %*%
        chol(exactPosterior$covariance[[i]])
      draws <- sweep(draws, 2L, exactPosterior$mean[i, ], "+")
      for (s in seq_len(nsamp)) {
        phiSamples[i, , s] <- start[i, ]
        phiSamples[i, idx, s] <- draws[s, ]
      }
    }
    psiSamples <- phiSamples
    for (s in seq_len(nsamp))
      psiSamples[, , s] <- transphi(matrix(phiSamples[, , s], nrow = N),
        object["model"]["transform.par"])
    object["results"]["phi.samp"] <- phiSamples
    object["results"]["psi.samp"] <- psiSamples
    object["results"]["phi.samp.var"] <- array(rep(phiVar, nsamp),
      dim = dim(phiSamples))
    object["results"]["cond.mean.phi"] <- phiMean
    object["results"]["cond.var.phi"] <- phiVar
    object["results"]["cond.mean.psi"] <- transphi(phiMean,
      object["model"]["transform.par"])
    object["results"]["cond.mean.eta"] <- phiMean -
      object["results"]["mean.phi"]
    shrink <- rep(NA_real_, ncol(phiMean))
    names(shrink) <- paste0("Sh.", object["model"]["name.modpar"], ".%")
    object["results"]["cond.shrinkage"] <- shrink
    attr(object, "saemix.copula.conditional") <-
      list(backend = "exact-linear-gaussian")
    return(object)
  }
  proposal<-as.matrix(state$proposalOmega)
  hasConditioning <- (state$dConditioning %||% 0L) > 0L
  gaussianConditionalProposal <- hasConditioning
  if(any(!is.finite(proposal))||inherits(try(chol(proposal),silent=TRUE),"try-error"))
    proposal<-diag(rep(.1^2,length(idx)),length(idx))
  cholProposal<-chol(proposal)*.45
  phiSamples<-array(NA_real_,c(N,ncol(start),nsamp))
  phiMean<-phiVar<-matrix(NA_real_,N,ncol(start));phiMean[]<-start;phiVar[]<-0
  logPrior<-function(p,i) {
    value <- try({
      eta <- matrix(p-predictor[i,],nrow=1)
      priorState <- if (hasConditioning)
        cbind(eta, as.matrix(state$conditioning)[i, , drop=FALSE]) else eta
      copulaGaussianFremLogPrior(priorState, state$vine,
        state$margins, state$dEta, "joint")
    }, silent=TRUE)
    if(inherits(value,"try-error") || length(value)!=1L || !is.finite(value))
      -Inf else as.numeric(value)
  }
  priorDraw<-function(i) {
    if (gaussianConditionalProposal) {
      eta <- copulaGaussianFremRandEta(
        as.matrix(state$conditioning)[i, , drop=FALSE], state$vine,
        state$margins, state$dEta)
      return(predictor[i,]+as.numeric(eta))
    }
    predictor[i,]+as.numeric(copulaMarginsQuantile(
      rvinecopulib::rvinecop(1L,state$vine),state$margins))
  }
  nkeep<-max(100L,nsamp);thin<-2L
  for(i in seq_len(N)) {
    current<-start[i,idx];row<-start[i,]
    row[idx]<-current;llY<-copulaSubjectResponseLogLik(object,i,row)
    lp<-logPrior(current,i);kept<-matrix(NA_real_,nkeep,length(idx));kk<-0L
    total<-burn+nkeep*thin
    for(iter in seq_len(total)) {
      ## For Gaussian-copula FREM this is an exact draw from
      ## p(phi_i | c_i,obs).  Its density cancels the fitted conditional prior
      ## in the independence-MH ratio, just as the unconditional prior does.
      independence<-(!hasConditioning || gaussianConditionalProposal) &&
        iter%%5L==0L
      candidate<-if(independence)priorDraw(i) else
        as.numeric(current+rnorm(length(idx))%*%cholProposal)
      row[idx]<-candidate;llC<-copulaSubjectResponseLogLik(object,i,row)
      lpC<-logPrior(candidate,i)
      logRatio<-if(independence)llC-llY else llC+lpC-llY-lp
      if(is.finite(logRatio)&&log(runif(1))<min(0,logRatio)) {
        current<-candidate;llY<-llC;lp<-lpC
      }
      if(iter>burn&&(iter-burn)%%thin==0L) {
        kk<-kk+1L;kept[kk,]<-current
      }
    }
    phiMean[i,idx]<-colMeans(kept);phiVar[i,idx]<-apply(kept,2,stats::var)
    choose<-round(seq(1,nkeep,length.out=nsamp))
    for(s in seq_len(nsamp)) {
      phiSamples[i,,s]<-start[i,];phiSamples[i,idx,s]<-kept[choose[s],]
    }
  }
  psiSamples<-phiSamples
  for(s in seq_len(nsamp)) {
    phs<-matrix(phiSamples[,,s],nrow=N,ncol=ncol(start))
    psiSamples[,,s]<-transphi(phs,object["model"]["transform.par"])
  }
  object["results"]["phi.samp"]<-phiSamples
  object["results"]["psi.samp"]<-psiSamples
  object["results"]["phi.samp.var"]<-array(rep(phiVar,nsamp),dim=dim(phiSamples))
  object["results"]["cond.mean.phi"]<-phiMean
  object["results"]["cond.var.phi"]<-phiVar
  meanPsi<-apply(psiSamples,c(1,2),mean)
  object["results"]["cond.mean.psi"]<-matrix(meanPsi,nrow=N,ncol=ncol(start))
  object["results"]["cond.mean.eta"]<-phiMean-object["results"]["mean.phi"]
  shrink<-rep(NA_real_,ncol(phiMean));names(shrink)<-paste0("Sh.",object["model"]["name.modpar"],".%")
  object["results"]["cond.shrinkage"]<-shrink
  if(isTRUE(plot))warning("copula conditional sampler does not yet draw convergence plots",call.=FALSE)
  object
}
