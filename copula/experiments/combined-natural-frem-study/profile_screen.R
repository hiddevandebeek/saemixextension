suppressPackageStartupMessages({
  library(devtools); devtools::load_all("C:/package/saemix-copula", quiet = TRUE)
  library(rvinecopulib)
})
root <- "C:/package/saemix-copula/copula/experiments/combined-natural-frem-study"
source(file.path(root, "functions.R"))
simulation <- combined_simulate(1401001L, 100L)
sxdata <- combined_data(simulation$data)
standard <- saemix(combined_model(), sxdata, combined_control(1401011L, 500L),
  population = combined_population(combined_initial_vine(),
    combined_standard_margins(simulation$crp), simulation$crp,
    "transformed-additive"))
start <- combined_fit_start(standard); state <- copulaGet(standard)
incumbent <- saemix(combined_model(start$fixed, start$residual, start$omega),
  sxdata, combined_control(1401021L, 500L),
  population = combined_population(state$vine,
    c(state$margins[1:2], list(combined_flexible_covariate_margin(simulation$crp))),
    simulation$crp, "transformed-additive"))
profile <- tempfile(fileext = ".out")
oldEnvironment <- new.env(parent = globalenv())
oldSource <- system2("git", c("show", "3d0faeb:R/naturalMarginSelection.R"),
  stdout = TRUE)
eval(parse(text = oldSource), envir = oldEnvironment)
arguments <- list(object = incumbent, supports = c("positive", "positive"),
  candidates = list("lognormal", c("lognormal", "gamma", "weibull")),
  posteriorDraws = 300L, max.iter = 300L, seed = 1401031L,
  optimizerMaxit = 200L, minimumEssFraction = 0,
  minimumPosteriorEssFraction = 0, maximumMcse = Inf)
oldElapsed <- system.time(oldSelection <- suppressWarnings(do.call(
  oldEnvironment$copulaSelectParameterMargins, arguments)))["elapsed"]
Rprof(profile, interval = .01)
elapsed <- system.time(selection <- suppressWarnings(do.call(
  copulaSelectParameterMargins, arguments)))["elapsed"]
Rprof(NULL)
profileSummary <- summaryRprof(profile)$by.total
print(selection$table[, c("families", "validation_bic_advantage")])
stopifnot(identical(oldSelection$families, selection$families),
  max(abs(oldSelection$table$validation_bic_advantage -
    selection$table$validation_bic_advantage)) < 1e-8)
print(head(profileSummary, 20L))
cat("screen elapsed seconds old/new/speedup:", oldElapsed, elapsed,
  oldElapsed / elapsed, "\n")
write.csv(data.frame(old_seconds = unname(oldElapsed),
  cached_seconds = unname(elapsed), speedup = unname(oldElapsed / elapsed),
  selected_family = paste(selection$families, collapse = "/"),
  maximum_bic_difference = max(abs(
    oldSelection$table$validation_bic_advantage -
      selection$table$validation_bic_advantage))),
  file.path(root, "screening_benchmark.csv"), row.names = FALSE)
