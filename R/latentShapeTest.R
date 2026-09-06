## Score test for the shape of a latent individual-parameter margin.
##
## The incumbent fit specifies each eta margin; on its Gaussian-score scale
## xi = Phi^-1{F_0(eta)} the margin is standard Normal. Embed it in the tilted
## family
##
##   f_c(x) = f_0(x) exp{c He_k(xi(x))} / Z(c),
##
## whose score at c = 0 is He_k(xi): the third Hermite polynomial is the
## skewness direction (Gamma and Weibull margins are left-skewed on the log
## scale), the fourth is the tail-weight direction (Student and Laplace
## positive, mixtures negative). This is Neyman's smooth test applied to a
## latent variable: by Fisher's identity the observed-data score for c is the
## posterior expectation of the complete-data score, which the incumbent's
## posterior draws supply, and no candidate family is fitted.
##
## Two details make it exact for the copula model. First, the coordinate's
## Gaussian score z_j moves with c, dz_j/dc = -He_{k-1}(z_j), so the copula
## quadratic form contributes -[(Omega - I) z]_j He_{k-1}(z_j); with the code's
## row convention zInfluence = z (I - Omega) that term is
## -zInfluence_j He_{k-1}(z_j). Second, the incumbent parameters are estimated,
## so the statistic uses the efficient score: the tilt scores are projected off
## the population-block scores of the incumbent (locations, scales,
## correlations) with the empirical outer-product information, which is the
## Delattre-Kuhn information formed from the same posterior means.
##
## Shrinkage is handled by construction. Thin individual data shrink each
## subject's posterior mean of He_3 by (1 - s)^3 and its variance across
## subjects by (1 - s)^6, so the statistic loses power exactly as the
## manuscript's attenuation law says, without any correction being applied.
##
## The residual-error parameters are projected off too, from centered
## differences of the response log likelihood at each draw; on data whose
## generator floored small predictions their score was the largest nuisance
## term, because the floored points inflate the fitted residual SD.

copulaHermite <- function(k, x) switch(as.character(as.integer(k)),
  `0` = rep(1, length(x)), `1` = x, `2` = x^2 - 1, `3` = x^3 - 3 * x,
  `4` = x^4 - 6 * x^2 + 3, `5` = x^5 - 10 * x^3 + 15 * x,
  `6` = x^6 - 15 * x^4 + 45 * x^2 - 15,
  `7` = x^7 - 21 * x^5 + 105 * x^3 - 105 * x,
  `8` = x^8 - 28 * x^6 + 210 * x^4 - 420 * x^2 + 105,
  stop("Hermite order must be between 0 and 8"))

## Tilt scores for every requested (coordinate, order) pair, per row.
copulaTiltRows <- function(z, Omega, coordinates, orders) {
  d <- ncol(z); zInfluence <- z %*% (diag(d) - Omega)
  answer <- NULL
  for (j in coordinates) for (k in orders)
    answer <- cbind(answer, copulaHermite(k, z[, j]) -
      zInfluence[, j] * copulaHermite(k - 1L, z[, j]))
  answer
}

## Minus the second derivative of the tilted complete log density in the tilt
## coefficients, per row. For one coordinate j and orders k, k',
##
##   -d2 l / dc_k dc_k' = (Omega - I)_jj He_{k-1} He_{k'-1}
##                        + [(Omega - I) z]_j z''_{kk'} + 1{k = k'} k!,
##
## where z''_{kk'} = -sum_r r! C(k,r) C(k',r) He_{k+k'-2r-1}(z_j)
##                   + z_j He_{k-1}(z_j) He_{k'-1}(z_j)
## comes from the second derivative of the tilted CDF through the Hermite
## product expansion. Across coordinates j != j' only the copula couples them:
##   -d2 l / dc_{jk} dc_{j'k'} = (Omega - I)_{jj'} He_{k-1}(z_j) He_{k'-1}(z_j').
## Both were checked against a numerically integrated tilted family.
copulaTiltHessianRows <- function(z, Omega, coordinates, orders) {
  d <- ncol(z); Om <- Omega - diag(d); Omz <- z %*% Om
  labels <- expand.grid(k = orders, j = coordinates)[, c("j", "k")]
  q <- nrow(labels); answer <- matrix(0, nrow(z), q * q)
  for (a in seq_len(q)) for (b in seq_len(q)) {
    j <- labels$j[a]; k <- labels$k[a]; jp <- labels$j[b]; kp <- labels$k[b]
    zj <- z[, j]
    value <- if (j == jp) {
      zpp <- zj * copulaHermite(k - 1L, zj) * copulaHermite(kp - 1L, zj)
      for (r in 0:min(k, kp)) {
        m <- k + kp - 2L * r
        if (m >= 1L) zpp <- zpp - factorial(r) * choose(k, r) * choose(kp, r) *
          copulaHermite(m - 1L, zj)
      }
      Om[j, j] * copulaHermite(k - 1L, zj) * copulaHermite(kp - 1L, zj) +
        Omz[, j] * zpp + (k == kp) * factorial(k)
    } else Om[j, jp] * copulaHermite(k - 1L, zj) *
      copulaHermite(kp - 1L, z[, jp])
    answer[, (a - 1L) * q + b] <- value
  }
  answer
}

copulaLatentShapeTest <- function(object, draws = 500L, max.iter = NULL,
    seed = NULL, thin = 2L, orders = c(3L, 4L), coordinates = NULL,
    posterior = NULL, bootstrap = 0L, bootstrapDraws = 200L,
    bootstrapBurn = 400L, bootstrapThin = 2L) {
  state <- copulaGet(object)
  n <- object["data"]["N"]
  if (is.null(posterior))
    posterior <- copulaPosteriorEtaDraws(object, draws, max.iter, seed, thin)
  if (any(dim(posterior)[1:2] != c(n, length(state$etaIndex))))
    stop("posterior draws do not match the fitted model")
  conditioning <- if ((state$dConditioning %||% 0L) > 0L)
    as.matrix(state$conditioning) else NULL
  ## Residual-error scores per subject and draw, by centered differences of
  ## the response log likelihood in each free residual parameter. The
  ## population-block scores are formed inside the core; these are the only
  ## nuisance scores that need the structural model.
  residual <- copulaLatentResidualScores(object, posterior)
  answer <- copulaLatentShapeTestCore(state$margins, state$vine, posterior,
    conditioning, orders = orders, coordinates = coordinates,
    populationScale = state$populationScale,
    diagnostics = attr(posterior, "diagnostics"),
    nuisanceExtra = residual)
  bootstrap <- as.integer(bootstrap)
  if (bootstrap > 0L) {
    ## Parametric bootstrap of the outer-product statistic under the fitted
    ## incumbent: data are simulated from the fitted model, the same posterior
    ## sampler and the same statistic are applied at the fitted parameters,
    ## and the observed statistic is referred to that distribution. This is
    ## the reference the asymptotic chi-square approximates, and at a few
    ## hundred subjects the approximation is loose for the outer-product form.
    reference <- matrix(NA_real_, bootstrap, nrow(answer$table))
    for (b in seq_len(bootstrap)) {
      simulated <- copulaLatentShapeSimulate(object, seed = if (is.null(seed))
        NULL else seed + 1000L * b)
      posteriorB <- copulaPosteriorEtaDraws(simulated, bootstrapDraws,
        bootstrapBurn, if (is.null(seed)) NULL else seed + 1000L * b + 1L,
        bootstrapThin)
      residualB <- copulaLatentResidualScores(simulated, posteriorB)
      tableB <- copulaLatentShapeTestCore(state$margins, state$vine,
        posteriorB, conditioning, orders = orders, coordinates = coordinates,
        populationScale = state$populationScale, nuisanceExtra = residualB,
        louis = FALSE)$table
      reference[b, ] <- tableB$statistic.opg
    }
    observed <- answer$table$statistic.opg
    answer$table$p.value.boot <- vapply(seq_along(observed), function(i)
      (1 + sum(reference[, i] >= observed[i])) / (bootstrap + 1), numeric(1))
    answer$bootstrap <- list(replicates = bootstrap, statistics = reference)
  }
  answer
}

## One data set simulated from the fitted incumbent: latent parameters from
## the fitted conditional population law, responses from the structural model
## and the fitted residual-error model. Returned as a copy of the fitted object
## carrying the simulated responses, so every downstream tool sees the same
## design, covariates and fitted parameters.
copulaLatentShapeSimulate <- function(object, seed = NULL) {
  state <- copulaGet(object); index <- as.integer(state$etaIndex)
  dEta <- length(index); n <- object["data"]["N"]
  layout <- copulaResponseLayout(object, 1L)
  withSeed(seed, {
    eta <- if ((state$dConditioning %||% 0L) > 0L)
      copulaGaussianFremConditionalKernel(as.matrix(state$conditioning),
        state$vine, state$margins, dEta)$random() else
      copulaMarginsQuantile(rvinecopulib::rvinecop(n, state$vine),
        state$margins)[, seq_len(dEta), drop = FALSE]
    phi <- object["results"]["mean.phi"]
    phi[, index] <- phi[, index, drop = FALSE] + eta
    psi <- transphi(phi, layout$transform)
    f <- layout$model["model"](psi, layout$IdB, layout$XB)
    if (!identical(layout$modeltype, "structural"))
      stop("the latent shape bootstrap requires a structural response model")
    g <- error(f, layout$respar, layout$XB$ytype)
    y <- f + g * stats::rnorm(length(f))
  })
  simulated <- object
  data <- simulated["data"]
  frame <- data["data"]; frame[, data["name.response"]] <- y
  simulated@data@data <- frame
  if (methods::.hasSlot(simulated@data, "yorig") &&
      length(simulated@data@yorig) == length(y)) simulated@data@yorig <- y
  simulated
}

copulaLatentResidualScores <- function(object, posterior, step = 1e-4) {
  state <- copulaGet(object); index <- as.integer(state$etaIndex)
  n <- dim(posterior)[1L]; samples <- dim(posterior)[3L]
  layout <- copulaResponseLayout(object, 1L)
  free <- copulaScoreResidualIndices(layout$errorModel)
  free <- free[free <= length(layout$respar)]
  if (!length(free)) return(NULL)
  centre <- object["results"]["mean.phi"][, index, drop = FALSE]
  phiTemplate <- object["results"]["mean.phi"]
  score <- hessian <- lapply(free, function(k) matrix(NA_real_, n, samples))
  names(score) <- names(hessian) <- paste0("residual.", free)
  for (k in seq_along(free)) {
    plus <- minus <- layout
    h <- step * max(1, abs(layout$respar[free[k]]))
    plus$respar[free[k]] <- plus$respar[free[k]] + h
    minus$respar[free[k]] <- minus$respar[free[k]] - h
    for (m in seq_len(samples)) {
      phi <- phiTemplate; phi[, index] <- centre + posterior[, , m]
      llPlus <- as.numeric(copulaResponseLogLikBatch(object, phi, 1L, plus))
      llMinus <- as.numeric(copulaResponseLogLikBatch(object, phi, 1L, minus))
      llCentre <- as.numeric(copulaResponseLogLikBatch(object, phi, 1L, layout))
      score[[k]][, m] <- (llPlus - llMinus) / (2 * h)
      hessian[[k]][, m] <- -(llPlus - 2 * llCentre + llMinus) / h^2
    }
  }
  list(score = score, hessian = hessian)
}

## The statistic from posterior draws of eta (subjects x coordinates x
## samples), the incumbent margins and vine, and the conditioning covariates
## (one row per subject). Separated from the fitted object so that it can be
## checked against exact posterior draws where those exist.
copulaLatentShapeTestCore <- function(margins, vine, posterior,
    conditioning = NULL, orders = c(3L, 4L), coordinates = NULL,
    populationScale = "transformed-additive", diagnostics = NULL,
    nuisanceExtra = NULL, louis = TRUE) {
  if (!identical(populationScale, "transformed-additive"))
    stop("the latent shape test is defined for a transformed-additive incumbent")
  orders <- sort(unique(as.integer(orders)))
  if (any(orders < 2L | orders > 6L))
    stop("orders must lie between 2 and 6")
  n <- dim(posterior)[1L]; dEta <- dim(posterior)[2L]
  samples <- dim(posterior)[3L]; d <- length(margins)
  coordinates <- if (is.null(coordinates)) seq_len(dEta) else
    as.integer(coordinates)
  if (any(coordinates < 1L | coordinates > dEta))
    stop("coordinates must index the eta margins")
  eta <- copulaFlattenEtaDraws(posterior)
  subject <- rep(seq_len(n), samples)
  if (!is.null(conditioning)) {
    conditioning <- as.matrix(conditioning)
    if (nrow(conditioning) != n || ncol(conditioning) != d - dEta)
      stop("conditioning covariates do not match the margins")
    if (anyNA(conditioning))
      stop("the latent shape test requires fully observed conditioning covariates")
    conditioning <- conditioning[subject, , drop = FALSE]
  } else if (d != dEta) stop("conditioning covariates are required")
  E <- if (is.null(conditioning)) eta else cbind(eta, conditioning)
  evaluated <- copulaGaussianFremEvaluateMargins(E, margins)
  if (!all(evaluated$valid))
    stop("a posterior draw fell outside the incumbent margins")
  R <- copulaGaussianRvineCor(vine, d); Omega <- solve(R)
  z <- evaluated$z

  ## Tilt scores per row.
  tilt <- copulaTiltRows(z, Omega, coordinates, orders)
  tiltLabel <- character(); tiltCoordinate <- integer(); tiltOrder <- integer()
  for (j in coordinates) for (k in orders) {
    tiltLabel <- c(tiltLabel, sprintf("eta%d.He%d", j, k))
    tiltCoordinate <- c(tiltCoordinate, j); tiltOrder <- c(tiltOrder, k)
  }
  colnames(tilt) <- tiltLabel

  ## Population-block scores of the incumbent per row, in native coordinates.
  ## Locations enter as a common shift of every eta coordinate (withMu), which
  ## is the intercept score whatever the structural design.
  layout <- copulaGaussianFremScoreLayout(margins, vine, d, dEta,
    withMu = TRUE)
  current <- copulaScoreToInternal(layout$native, layout$lower, layout$upper)
  w <- rep(1 / nrow(E), nrow(E))
  complete <- copulaGaussianFremCompleteScoreInternal(current, layout, E, w,
    1e-5, perRow = TRUE)
  if (is.null(complete) || is.null(complete$nativeRows))
    stop("the incumbent complete-data score could not be evaluated on the draws")
  ## Internal coordinates throughout: the projection is invariant to the
  ## nuisance parameterization and the Louis Hessian below is formed by
  ## perturbing the internal vector.
  internalScale <- copulaScoreInternalDerivative(current, layout$lower,
    layout$upper)
  nuisanceRows <- sweep(complete$nativeRows * nrow(E), 2L, internalScale, "*")
  nuisanceLabel <- c(paste0("location.", seq_len(layout$nLocation)),
    if (layout$nMargin) unlist(lapply(seq_len(d), function(j)
      if (!length(layout$marginLayout$index[[j]])) character() else
        paste0("margin.", j, ".", seq_along(layout$marginLayout$index[[j]]))))
    else character(),
    if (layout$nEdge) paste0("correlation.", seq_len(layout$nEdge)) else
      character())

  ## Posterior means per subject: the observed-data score contributions.
  perSubject <- function(rows) rowsum(rows, subject, reorder = TRUE) / samples
  tiltSubject <- perSubject(tilt)
  nuisanceSubject <- perSubject(nuisanceRows)
  colnames(nuisanceSubject) <- nuisanceLabel
  extraRows <- NULL; extraHessian <- NULL
  if (!is.null(nuisanceExtra)) {
    scores <- if (!is.null(nuisanceExtra$score)) nuisanceExtra$score else
      nuisanceExtra
    extraRows <- do.call(cbind, lapply(scores, function(m) {
      m <- as.matrix(m)
      if (any(dim(m) != c(n, samples)))
        stop("extra nuisance scores must be subjects x draws")
      as.vector(m)
    }))
    colnames(extraRows) <- names(scores)
    nuisanceSubject <- cbind(nuisanceSubject, perSubject(extraRows))
    if (!is.null(nuisanceExtra$hessian))
      extraHessian <- vapply(nuisanceExtra$hessian, function(m)
        sum(as.matrix(m)) / samples, numeric(1))
  }
  S <- cbind(nuisanceSubject, tiltSubject)
  information <- crossprod(S)
  q <- ncol(tiltSubject); p <- ncol(nuisanceSubject)
  nuisanceIndex <- seq_len(p); tiltIndex <- p + seq_len(q)
  Itt <- information[tiltIndex, tiltIndex, drop = FALSE]
  Itn <- information[tiltIndex, nuisanceIndex, drop = FALSE]
  Inn <- information[nuisanceIndex, nuisanceIndex, drop = FALSE]
  InnInverse <- tryCatch(chol2inv(chol(Inn)), error = function(e)
    MASS::ginv(Inn))
  efficient <- Itt - Itn %*% InnInverse %*% t(Itn)
  efficient <- (efficient + t(efficient)) / 2
  scoreTotal <- colSums(tiltSubject)
  ## Efficient score: the tilt score with its regression on the nuisance
  ## scores removed. At an exact maximum the nuisance scores sum to zero and
  ## this is the plain tilt score; a stochastic-approximation estimate is not
  ## exactly at the maximum, and the projection removes that first-order
  ## effect from the numerator as well as from the variance.
  nuisanceTotal <- colSums(nuisanceSubject)
  scoreEfficient <- as.numeric(scoreTotal - Itn %*% InnInverse %*% nuisanceTotal)
  names(scoreEfficient) <- tiltLabel

  ## Monte Carlo error of the score from two halves of the draws.
  half <- samples %/% 2L
  firstHalf <- subject <= n & rep(seq_len(samples), each = n) <= half
  tiltFirst <- rowsum(tilt[firstHalf, , drop = FALSE],
    subject[firstHalf]) / half
  tiltSecond <- rowsum(tilt[!firstHalf, , drop = FALSE],
    subject[!firstHalf]) / (samples - half)
  mcse <- abs(colSums(tiltFirst) - colSums(tiltSecond)) / 2

  ## Louis observed information from the same draws,
  ##   I = sum_i { E[-d2 l | y_i] - Var(d l | y_i) },
  ## with the expected complete Hessian formed analytically in the tilt block
  ## and by centered differences of the per-row scores in the population
  ## block, and the posterior covariance from the rows. The outer-product
  ## form above is consistent for the same matrix but over-rejects in samples
  ## of a few hundred subjects (Davidson and MacKinnon 1983); the Louis form
  ## is what the reported p-values use.
  pTotal <- ncol(nuisanceSubject); pPop <- ncol(nuisanceRows)
  allRows <- cbind(nuisanceRows, if (is.null(extraRows)) NULL else extraRows,
    tilt)
  Sall <- cbind(nuisanceSubject, tiltSubject)
  posteriorCovariance <- crossprod(allRows) / samples - crossprod(Sall)
  hessian <- if (isTRUE(louis)) matrix(0, pTotal + q, pTotal + q) else NULL
  popIndex <- seq_len(pPop); allTilt <- pTotal + seq_len(q)
  stepInternal <- 1e-4
  for (b in if (isTRUE(louis)) popIndex else integer()) {
    plus <- minus <- current
    plus[b] <- plus[b] + stepInternal; minus[b] <- minus[b] - stepInternal
    gp <- copulaGaussianFremCompleteScoreInternal(plus, layout, E, w, 1e-5,
      perRow = TRUE)
    gm <- copulaGaussianFremCompleteScoreInternal(minus, layout, E, w, 1e-5,
      perRow = TRUE)
    if (is.null(gp) || is.null(gm)) { hessian <- NULL; break }
    popPlus <- colSums(gp$nativeRows) * nrow(E) *
      copulaScoreInternalDerivative(plus, layout$lower, layout$upper)
    popMinus <- colSums(gm$nativeRows) * nrow(E) *
      copulaScoreInternalDerivative(minus, layout$lower, layout$upper)
    tiltPlus <- colSums(copulaTiltRows(gp$evaluated$z,
      solve(gp$candidate$correlation), coordinates, orders))
    tiltMinus <- colSums(copulaTiltRows(gm$evaluated$z,
      solve(gm$candidate$correlation), coordinates, orders))
    column <- -c(popPlus - popMinus, rep(0, pTotal - pPop),
      tiltPlus - tiltMinus) / (2 * stepInternal * samples)
    hessian[, b] <- column; hessian[b, ] <- column
  }
  if (!is.null(hessian)) {
    tiltHessian <- colSums(copulaTiltHessianRows(z, Omega, coordinates,
      orders)) / samples
    hessian[allTilt, allTilt] <- matrix(tiltHessian, q, q, byrow = TRUE)
    if (pTotal > pPop && !is.null(extraHessian)) {
      extraIndex <- pPop + seq_len(pTotal - pPop)
      ## The response and population parts of the complete log density are
      ## additively separable, so their complete cross derivatives vanish.
      hessian[cbind(extraIndex, extraIndex)] <- extraHessian
    }
    louis <- hessian - posteriorCovariance
    louis <- (louis + t(louis)) / 2
    InnLouis <- louis[seq_len(pTotal), seq_len(pTotal), drop = FALSE]
    InnLouisInverse <- tryCatch(chol2inv(chol(InnLouis)),
      error = function(e) NULL)
    efficientLouis <- if (is.null(InnLouisInverse)) NULL else
      louis[allTilt, allTilt] - louis[allTilt, seq_len(pTotal)] %*%
        InnLouisInverse %*% louis[seq_len(pTotal), allTilt]
    if (!is.null(efficientLouis))
      efficientLouis <- (efficientLouis + t(efficientLouis)) / 2
    scoreLouis <- if (is.null(InnLouisInverse)) NULL else as.numeric(
      scoreTotal - louis[allTilt, seq_len(pTotal)] %*% InnLouisInverse %*%
        nuisanceTotal)
  } else { louis <- NULL; efficientLouis <- NULL; scoreLouis <- NULL }
  louisUsable <- !is.null(efficientLouis) &&
    all(eigen(efficientLouis, symmetric = TRUE, only.values = TRUE)$values > 0)

  ## The Louis block is used when it is positive definite for the directions
  ## tested. Away from the null the observed information need not be, and
  ## the outer-product form, which is always positive semi-definite, is used
  ## instead; the `information` column of the table says which.
  louisBlockUsable <- function(which) {
    if (is.null(efficientLouis)) return(FALSE)
    block <- efficientLouis[which, which, drop = FALSE]
    all(eigen(block, symmetric = TRUE, only.values = TRUE)$values > 0)
  }
  statistic <- function(which, form = c("auto", "louis", "opg")) {
    form <- match.arg(form)
    if (form == "auto") form <- if (louisBlockUsable(which)) "louis" else "opg"
    if (form == "louis" && !louisBlockUsable(which)) return(c(statistic =
      NA_real_, df = length(which), p.value = NA_real_, louis = 0))
    Ieff <- if (form == "louis") efficientLouis[which, which, drop = FALSE] else
      efficient[which, which, drop = FALSE]
    s <- if (form == "louis") scoreLouis[which] else scoreEfficient[which]
    inverse <- tryCatch(chol2inv(chol(Ieff)), error = function(e)
      MASS::ginv(Ieff))
    value <- as.numeric(t(s) %*% inverse %*% s)
    c(statistic = value, df = length(which),
      p.value = stats::pchisq(value, length(which), lower.tail = FALSE),
      louis = as.numeric(form == "louis"))
  }
  rows <- list()
  for (j in coordinates) {
    for (k in orders) {
      which <- which(tiltCoordinate == j & tiltOrder == k)
      st <- statistic(which); so <- statistic(which, "opg")
      rows[[length(rows) + 1L]] <- data.frame(coordinate = j,
        test = sprintf("He%d", k), statistic = st[["statistic"]],
        df = st[["df"]], p.value = st[["p.value"]],
        information = if (st[["louis"]] > 0) "louis" else "opg",
        statistic.opg = so[["statistic"]], p.value.opg = so[["p.value"]],
        score = scoreEfficient[which], score.raw = scoreTotal[which],
        score.se = sqrt(efficient[which, which]),
        score.se.louis = if (louisBlockUsable(which))
          sqrt(efficientLouis[which, which]) else NA_real_,
        score.mcse = mcse[which],
        attenuation = if (k == 3L) stats::var(tiltSubject[, which]) /
          factorial(3L) else NA_real_)
    }
    if (length(orders) > 1L) {
      which <- which(tiltCoordinate == j)
      st <- statistic(which); so <- statistic(which, "opg")
      rows[[length(rows) + 1L]] <- data.frame(coordinate = j,
        test = paste0("joint(", paste0("He", orders, collapse = ","), ")"),
        statistic = st[["statistic"]], df = st[["df"]],
        p.value = st[["p.value"]],
        information = if (st[["louis"]] > 0) "louis" else "opg",
        statistic.opg = so[["statistic"]],
        p.value.opg = so[["p.value"]], score = NA_real_,
        score.raw = NA_real_, score.se = NA_real_, score.se.louis = NA_real_,
        score.mcse = NA_real_, attenuation = NA_real_)
    }
  }
  table <- do.call(rbind, rows)
  rownames(table) <- NULL
  unadjusted <- vapply(seq_len(q), function(i)
    scoreTotal[i]^2 / Itt[i, i], numeric(1))
  structure(list(table = table, score = scoreEfficient, scoreRaw = scoreTotal,
    nuisanceScore = nuisanceTotal, scoreLouis = scoreLouis,
    efficientInformation = efficient, information = information,
    louisInformation = louis, efficientLouisInformation = efficientLouis,
    louisUsable = louisUsable,
    unadjustedStatistic = stats::setNames(unadjusted, tiltLabel),
    subjectScores = tiltSubject, nuisanceScores = nuisanceSubject,
    draws = samples, subjects = n, orders = orders,
    coordinates = coordinates,
    diagnostics = diagnostics,
    limitations = paste0("nuisance scores are projected off with the ",
      "outer-product information; Monte Carlo error of the score is ",
      "reported as score.mcse")),
    class = "saemixLatentShapeTest")
}

print.saemixLatentShapeTest <- function(x, digits = 4L, ...) {
  cat("Latent margin shape test (Neyman smooth test on the incumbent posterior)\n")
  cat(sprintf("%d subjects, %d posterior draws per subject; orders %s\n",
    x$subjects, x$draws, paste(x$orders, collapse = ", ")))
  print(format(x$table[, c("coordinate", "test", "statistic", "df", "p.value",
    "statistic.opg", "score", "score.se.louis", "score.mcse", "attenuation")],
    digits = digits), row.names = FALSE)
  if (!isTRUE(x$louisUsable))
    cat("Louis information not positive definite here; reported statistics are NA and the outer-product form is in statistic.opg.\n")
  cat("Negative He3 score: left-skewed on the incumbent's Gaussian-score scale ",
    "(Gamma, Weibull on the log scale); positive He4: heavier tails; ",
    "negative He4: lighter shoulders or a mixture.\n", sep = "")
  cat("attenuation = Var(posterior mean of He3) / 6, about (1 - s)^6 for ",
    "eta-shrinkage s.\n", sep = "")
  invisible(x)
}
