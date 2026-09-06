## Small integration gate: ordinary Gaussian saemix -> automatic joint margin
## learning -> fixed Gaussian-copula FREM population declaration.
suppressPackageStartupMessages({
  library(devtools)
  load_all("C:/package/saemix-copula", quiet = TRUE)
})
sourceFile <- paste0("C:/package/saemix-copula/copula/examples/",
  "gaussian-copula-frem/run_example.R")
code <- readLines(sourceFile, warn = FALSE)
cut <- grep("^simulation_path <-", code)[1L]
eval(parse(text = code[seq_len(cut - 1L)]), envir = environment())
simulation <- readRDS(file.path(out_dir, "simulation.rds"))
control <- list(seed = 82751L, save = FALSE, save.graphs = FALSE,
  print = FALSE, displayProgress = FALSE, warnings = FALSE,
  nbiter.saemix = c(20, 10), nbiter.mcmc = c(1, 1, 1, 0),
  ll.is = FALSE, fim = FALSE, map = FALSE)
gaussian <- saemix(make_saemix_model(), make_saemix_data(simulation$data), control)
learning <- learnGaussianCopulaFremMargins(gaussian,
  simulation$conditioning, nsamp = 3L, etaFamilies = "normal",
  covariateFamilies = c("normal", "lognormal"), cycles = 1L,
  jointMaxit = 10L, seed = 82752L)
stopifnot(inherits(learning, "saemixGaussianCopulaFremLearning"),
  inherits(learning$population, "saemixPopulation"),
  length(learning$selectedFamilies) == 4L,
  all(learning$ranking$backend == "gaussian-copula-frem-mvnormal-em"),
  isTRUE(learning$selection$commonMeasure))
cat("automatic Gaussian-copula FREM workflow check passed\n")
