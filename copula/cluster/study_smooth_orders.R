## Preliminary: does a data-driven smooth test beat the fixed skewness and
## tail-weight pair?
##
## The latent shape test currently fixes two directions, He3 and He4. Expanding
## the average posterior of the Gaussian scores in the Hermite basis gives
## Delta(b) - 1 = sum_k c_k He_k(b) / k!, and the Cramer-von Mises distance of
## the gradient function of Verbeke and Molenberghs (2013) from one is the sum
## of the squared normalised coefficients. So the omnibus gradient statistic and
## Neyman's smooth statistic are the same object, and the only question is where
## to truncate. Ledwina's data-driven rule chooses the truncation by a Schwarz
## penalty and is consistent against essentially all smooth alternatives.
##
## This chunk driver computes, per coordinate and from one set of posterior
## draws and one parametric bootstrap:
##
##   S4      the present statistic, orders 3 and 4
##   S6      a fixed omnibus statistic, orders 3 to 6
##   Sdd     the data-driven statistic, orders 3..Khat with
##           Khat = argmax_K { S_K - (K - 2) log N }
##
## and bootstrap p-values for all three, using the same replicates so that the
## comparison is paired. Under the null the Hermite directions are orthogonal,
## so S_K is taken as the sum of the per-order efficient statistics; the
## bootstrap carries whatever correlation remains.
##
## Three scenarios:
##   alt         Gamma clearance, lognormal volume  (power on coordinate 2,
##               size on coordinate 1)
##   null        both parameter margins lognormal, incumbent keeps the Normal
##               CRP margin      (size under a misspecified covariate margin)
##   nullexact   both lognormal and the incumbent given the correct Gamma CRP
##               margin          (size under a correctly specified joint model)
##
## Environment: REPS SEED0 TAG OUTDIR PKG DEV FUNCTIONS SCENARIO ITERS BOOTSTRAP

suppressPackageStartupMessages({
  dev <- Sys.getenv("DEV", ""); pkg <- Sys.getenv("PKG", "")
  if (nzchar(dev)) devtools::load_all(dev, quiet = TRUE) else
    if (nzchar(pkg)) library(saemix, lib.loc = pkg) else library(saemix)
  library(rvinecopulib)
})
`%||%` <- function(x, y) if (is.null(x)) y else x

REPS  <- as.integer(Sys.getenv("REPS", "10"))
SEED0 <- as.integer(Sys.getenv("SEED0", "0"))
TAG   <- Sys.getenv("TAG", "s01")
OUT   <- Sys.getenv("OUTDIR", "out_smooth_orders")
FUN   <- Sys.getenv("FUNCTIONS", "functions.R")
SCEN  <- Sys.getenv("SCENARIO", "alt")
ITERS <- as.integer(Sys.getenv("ITERS", "3000"))
BOOT  <- as.integer(Sys.getenv("BOOTSTRAP", "59"))
ORDERS <- 3:6
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(FUN)) stop("cannot find functions.R at ", FUN)
source(FUN)
stopifnot(SCEN %in% c("alt", "null", "nullexact"))

## The null generator: the same population with the Gamma clearance replaced by
## a lognormal of equal log-scale variance, so only the shape differs.
simulate_null <- function(seed, n) {
  set.seed(seed)
  truth <- combined_truth(); R <- combined_correlation(truth)
  z <- matrix(rnorm(n * 3L), n, 3L) %*% chol(R)
  margins <- list(copulaNaturalMarginLognormal(truth$sdlogV),
    copulaNaturalMarginLognormal(sqrt(trigamma(truth$gammaCL))),
    copulaMarginCovariateGamma(truth$gammaCRP, truth$scaleCRP))
  typical <- matrix(rep(c(truth$V, truth$CL), each = n), nrow = n)
  psi <- copulaNaturalMarginsQuantile(pnorm(z[, 1:2, drop = FALSE]), typical,
    margins[1:2])
  crp <- margins[[3L]]$quantile(pnorm(z[, 3L]), margins[[3L]]$parameters)
  data <- data.frame(id = rep(seq_len(n), each = length(truth$times)),
    dose = truth$dose, time = rep(truth$times, n))
  data$y <- combined_pk(psi, data$id, cbind(data$dose, data$time)) *
    (1 + truth$residual * rnorm(nrow(data)))
  list(data = data, crp = crp, truth = truth)
}

## Smooth statistics of increasing order and the Schwarz choice among them.
smooth_statistics <- function(perOrder, n) {
  cumulative <- cumsum(perOrder)
  penalised <- cumulative - (seq_along(perOrder) - 1L) * log(n)
  khat <- which.max(penalised)
  c(S4 = cumulative[2L], S6 = cumulative[length(cumulative)],
    Sdd = cumulative[khat], khat = ORDERS[khat])
}

for (k in seq_len(REPS)) {
  i <- SEED0 + k
  path <- file.path(OUT, sprintf("%s_%03d.rds", SCEN, i))
  if (file.exists(path)) { cat(sprintf("[%s] %d done\n", TAG, i)); next }
  started <- Sys.time()
  answer <- tryCatch({
    seedBase <- 2400000L + 1000L * i
    simulation <- if (identical(SCEN, "alt"))
      combined_simulate(seedBase + 1L, 250L) else
      simulate_null(seedBase + 1L, 250L)
    sxdata <- combined_data(simulation$data)
    ## The incumbent: Normal CRP margin, except in the exact-null scenario
    ## where it is given the correct Gamma covariate margin.
    margins <- if (identical(SCEN, "nullexact"))
      list(copulaMarginNormal(.25), copulaMarginNormal(.40),
        copulaFitCovariateMargin(simulation$crp, "gamma")) else
      combined_standard_margins(simulation$crp)
    population <- combined_population(combined_initial_vine(), margins,
      simulation$crp, "transformed-additive")
    fit <- saemix(combined_model(), sxdata,
      combined_control(seedBase + 11L, ITERS), population = population)
    ## The exported entry point also forms Louis' observed information, whose
    ## cross terms need Hermite polynomials of order 2K; beyond K = 4 that
    ## exceeds the implemented range. The statistic itself does not need it, so
    ## the core is called directly with louis = FALSE, exactly as the package's
    ## own bootstrap path does.
    state <- copulaGet(fit)
    conditioning <- if ((state$dConditioning %||% 0L) > 0L)
      as.matrix(state$conditioning) else NULL
    core <- function(object, draws, burn, seed, thin) {
      posterior <- saemix:::copulaPosteriorEtaDraws(object, draws, burn, seed, thin)
      residual <- saemix:::copulaLatentResidualScores(object, posterior)
      saemix:::copulaLatentShapeTestCore(state$margins, state$vine, posterior,
        conditioning, orders = ORDERS, coordinates = NULL,
        populationScale = state$populationScale, nuisanceExtra = residual,
        louis = FALSE)
    }
    observedCore <- core(fit, 400L, 800L, seedBase + 77L, 4L)
    table <- observedCore$table
    reference <- matrix(NA_real_, BOOT, nrow(table))
    for (b in seq_len(BOOT)) {
      simulated <- saemix:::copulaLatentShapeSimulate(fit,
        seed = seedBase + 77L + 1000L * b)
      reference[b, ] <- core(simulated, 200L, 400L,
        seedBase + 77L + 1000L * b + 1L, 2L)$table$statistic.opg
    }
    n <- fit["data"]["N"]
    rows <- lapply(1:2, function(j) {
      keep <- which(table$coordinate == j & table$test != "joint" &
        !grepl("joint", table$test))
      keep <- keep[order(as.integer(sub("He", "", table$test[keep])))]
      stopifnot(length(keep) == length(ORDERS))
      observed <- smooth_statistics(table$statistic.opg[keep], n)
      boot <- t(apply(reference[, keep, drop = FALSE], 1L,
        function(x) smooth_statistics(x, n)))
      pvalue <- function(name) (1 + sum(boot[, name] >= observed[[name]])) /
        (nrow(boot) + 1)
      data.frame(scenario = SCEN, replicate = i, coordinate = j,
        S4 = observed[["S4"]], S6 = observed[["S6"]], Sdd = observed[["Sdd"]],
        khat = observed[["khat"]],
        p_S4 = pvalue("S4"), p_S6 = pvalue("S6"), p_Sdd = pvalue("Sdd"),
        p_He3 = (1 + sum(reference[, keep[1L]] >= table$statistic.opg[keep[1L]])) /
          (BOOT + 1),
        p_He4 = (1 + sum(reference[, keep[2L]] >= table$statistic.opg[keep[2L]])) /
          (BOOT + 1),
        khat_boot_median = stats::median(boot[, "khat"]))
    })
    result <- list(scenario = SCEN, replicate = i, iterations = ITERS,
      bootstrap = BOOT, orders = ORDERS, subjects = n,
      summary = do.call(rbind, rows), table = table,
      elapsedSeconds = as.numeric(difftime(Sys.time(), started, units = "secs")))
    temporary <- paste0(path, ".tmp-", Sys.getpid()); saveRDS(result, temporary)
    if (!file.rename(temporary, path)) stop("atomic rename failed")
    copulaClear()
    "ok"
  }, error = function(e) paste("FAILED:", conditionMessage(e)))
  cat(sprintf("[%s] %s %d %s (%.1f min)\n", TAG, SCEN, i, answer,
    as.numeric(difftime(Sys.time(), started, units = "mins"))))
  flush.console()
}
cat(sprintf("[%s] chunk finished\n", TAG))
