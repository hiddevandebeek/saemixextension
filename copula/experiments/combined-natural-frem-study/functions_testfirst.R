## Test-first variant of the combined natural-parameter and covariate-margin
## study: the workflow the manuscript proposes, evaluated end to end.
##
## The frozen 500-replicate study ranks candidate clearance families directly
## from the incumbent's posterior draws and holds volume lognormal. That is the
## ranking-and-refit route. Here the score test of the latent margins runs
## first on the standard fit, and a coordinate is offered to the ranking only
## where a direction rejects; a coordinate the test clears keeps its lognormal
## margin, and if the test clears both, no ranking is run at all.
##
## Everything else -- the generator, the seeds, the standard arm, the covariate
## screening, the refit and every reported endpoint -- is that of the frozen
## study, so replicate i here and replicate i there are the same dataset and
## the same Gaussian comparator.
##
## Source functions.R first; this file adds the test-first entry points and
## leaves the originals untouched.

## Which coordinates the ranking may change: those where a direction of the
## score test rejects at `level`.
combined_test_first_candidates <- function(test, dEta = 2L, level = .05,
                                           registry = c("lognormal", "gamma",
                                             "weibull")) {
  table <- test$table
  if (!is.data.frame(table) || !"p.value.boot" %in% names(table))
    stop("the score test must be bootstrap-calibrated for the test-first workflow")
  ## The directional tests only. The core also reports a joint row, but the
  ## workflow acts on a direction, and the reported calibration is directional.
  rejects <- vapply(seq_len(dEta), function(j) {
    rows <- table$coordinate == j & table$test %in% c("He3", "He4")
    any(rows) && any(table$p.value.boot[rows] < level)
  }, logical(1))
  list(rejects = rejects, level = level,
    candidates = lapply(rejects, function(r) if (r) registry else "lognormal"))
}

combined_fit_pair_test_first <- function(simulation, seedBase,
                                         iterations = 1500L,
                                         screenDraws = 700L,
                                         screenRetry = 2000L,
                                         likelihoodDraws = 2500L,
                                         testDraws = 400L, testThin = 4L,
                                         testMaxIter = 800L,
                                         bootstrap = 59L,
                                         bootstrapDraws = 200L,
                                         level = .05) {
  sxdata <- combined_data(simulation$data)
  standardPopulation <- combined_population(combined_initial_vine(),
    combined_standard_margins(simulation$crp), simulation$crp,
    "transformed-additive")
  standardTime <- system.time(standard <- saemix(combined_model(), sxdata,
    combined_control(seedBase + 11L, iterations),
    population = standardPopulation))["elapsed"]
  start <- combined_fit_start(standard)

  ## The workflow's first question, asked of the standard fit and nothing else:
  ## should either margin bend, and in which direction.
  testTime <- system.time(test <- copulaLatentShapeTest(standard,
    draws = testDraws, seed = seedBase + 77L, thin = testThin,
    max.iter = testMaxIter, bootstrap = bootstrap,
    bootstrapDraws = bootstrapDraws, bootstrapBurn = 400L,
    bootstrapThin = 2L))["elapsed"]
  decision <- combined_test_first_candidates(test, dEta = 2L, level = level)

  ## The covariate is observed, so its margin is screened from its values as in
  ## the frozen study; the score test governs the latent parameter margins only.
  flexibleCovariate <- combined_flexible_covariate_margin(simulation$crp)
  standardState <- copulaGet(standard)
  incumbentMargins <- c(standardState$margins[1:2], list(flexibleCovariate))
  incumbentPopulation <- combined_population(standardState$vine,
    incumbentMargins, simulation$crp, "transformed-additive")
  incumbentTime <- system.time(incumbent <- saemix(
    combined_model(start$fixed, start$residual, start$omega), sxdata,
    combined_control(seedBase + 21L, iterations),
    population = incumbentPopulation))["elapsed"]

  anyReject <- any(decision$rejects)
  selection <- NULL; selectionTime <- 0
  if (anyReject) {
    select <- function(draws, maxIter, seed) suppressWarnings(
      copulaSelectParameterMargins(incumbent,
        supports = c("positive", "positive"),
        candidates = decision$candidates,
        posteriorDraws = draws, max.iter = maxIter, seed = seed,
        optimizerMaxit = 200L, minimumEssFraction = .005,
        minimumPosteriorEssFraction = .01, maximumMcse = 1))
    selectionTime <- system.time(selection <- select(screenDraws, 450L,
      seedBase + 31L))["elapsed"]
    if (!isTRUE(selection$diagnostics$selectionResolved)) {
      retry <- system.time(selection <- select(screenRetry, 700L,
        seedBase + 32L))["elapsed"]
      selectionTime <- selectionTime + retry
    }
  }
  ## No rejection, or a rejection the ranking could not resolve: the standard
  ## families are kept and only the covariate margin has moved.
  retainStandard <- !anyReject ||
    !isTRUE(selection$diagnostics$selectionResolved)
  selectedFamilies <- if (retainStandard) c("lognormal", "lognormal") else
    selection$families
  finalPopulation <- if (retainStandard)
    combined_natural_lognormal_population(copulaGet(incumbent), simulation$crp) else
    selection$population
  finalStart <- combined_fit_start(incumbent)
  flexibleTime <- system.time(flexible <- saemix(
    combined_model(finalStart$fixed, finalStart$residual, finalStart$omega),
    sxdata, combined_control(seedBase + 41L, iterations),
    population = finalPopulation))["elapsed"]

  for (arm in c("standard", "flexible")) {
    fit <- get(arm); fit@options$nmc.is <- likelihoodDraws
    fit <- suppressWarnings(llisCopula.saemix(fit, defensive = 0,
      batch = 100L, seed = seedBase + if (arm == "standard") 51L else 52L))
    assign(arm, fit)
  }
  if (is.null(selection)) selection <- list(table = NULL, diagnostics = NULL,
    families = c("lognormal", "lognormal"))
  selection$retainedStandard <- retainStandard
  selection$rankingRun <- anyReject
  list(standard = standard, incumbent = incumbent, flexible = flexible,
    test = test, decision = decision, selection = selection,
    selectedFamilies = selectedFamilies,
    covariateFamily = flexibleCovariate$name,
    elapsed = c(standard = unname(standardTime), test = unname(testTime),
      incumbent = unname(incumbentTime), selection = unname(selectionTime),
      flexible = unname(flexibleTime)))
}

## Side channel: the frozen replicate driver builds its own record and does not
## know about the score test, so the fit stage leaves the test here and the
## driver's output is amended with it afterwards.
.combinedTestFirst <- new.env(parent = emptyenv())

combined_run_replicate_test_first <- function(replicate, resultRoot,
                                              force = FALSE, nSubjects = 250L,
                                              nVpc = 50000L,
                                              iterations = 1500L,
                                              screenDraws = 700L,
                                              screenRetry = 2000L,
                                              likelihoodDraws = 2500L,
                                              bootstrap = 59L) {
  replicate <- as.integer(replicate)
  path <- file.path(resultRoot, sprintf("replicate_%03d.rds", replicate))
  if (file.exists(path) && !force) return(path)
  ## combined_run_replicate resolves combined_fit_pair lexically in the
  ## environment functions.R was sourced into; swap it there for the duration.
  env <- environment(combined_run_replicate)
  original <- get("combined_fit_pair", envir = env)
  on.exit(assign("combined_fit_pair", original, envir = env), add = TRUE)
  assign("combined_fit_pair", function(simulation, seedBase, iterations,
      screenDraws, screenRetry, likelihoodDraws) {
    fits <- combined_fit_pair_test_first(simulation, seedBase, iterations,
      screenDraws, screenRetry, likelihoodDraws, bootstrap = bootstrap)
    .combinedTestFirst$test <- fits$test
    .combinedTestFirst$decision <- fits$decision
    .combinedTestFirst$elapsed <- fits$elapsed
    .combinedTestFirst$rankingRun <- isTRUE(fits$selection$rankingRun)
    fits
  }, envir = env)
  combined_run_replicate(replicate = replicate, resultRoot = resultRoot,
    force = force, nSubjects = nSubjects, nVpc = nVpc,
    iterations = iterations, screenDraws = screenDraws,
    screenRetry = screenRetry, likelihoodDraws = likelihoodDraws)

  result <- readRDS(path)
  result$schema <- 2L
  result$workflow <- "test-first"
  result$shapeTest <- .combinedTestFirst$test$table
  result$shapeTestDiagnostics <- .combinedTestFirst$test$diagnostics
  result$shapeTestRejects <- .combinedTestFirst$decision$rejects
  result$shapeTestLevel <- .combinedTestFirst$decision$level
  result$rankingRun <- .combinedTestFirst$rankingRun
  result$stageSeconds <- .combinedTestFirst$elapsed
  temporary <- paste0(path, ".amend-", Sys.getpid())
  saveRDS(result, temporary)
  if (!file.rename(temporary, path))
    stop("could not atomically amend ", path)
  path
}
