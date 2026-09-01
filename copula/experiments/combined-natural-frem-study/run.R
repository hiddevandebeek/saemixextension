args <- commandArgs(trailingOnly = TRUE)
replicates <- if (length(args)) as.integer(args[1L]) else 10L
workers <- if (length(args) >= 2L) as.integer(args[2L]) else 4L
mode <- if (length(args) >= 3L) args[3L] else "final"
force <- length(args) >= 4L && identical(tolower(args[4L]), "force")
if (!mode %in% c("pilot", "final")) stop("mode must be pilot or final")
root <- "C:/package/saemix-copula/copula/experiments/combined-natural-frem-study"
resultRoot <- file.path(root, if (mode == "pilot") "pilot" else "replicates")
dir.create(resultRoot, recursive = TRUE, showWarnings = FALSE)
functionFile <- file.path(root, "functions.R")
paths <- file.path(resultRoot, sprintf("replicate_%03d.rds", seq_len(replicates)))
pending <- seq_len(replicates)[force | !file.exists(paths)]
cat("mode:", mode, "replicates:", replicates, "pending:", length(pending),
  "workers:", workers, "\n")
if (!length(pending)) quit(save = "no", status = 0L)
cluster <- parallel::makePSOCKcluster(min(workers, length(pending)),
  outfile = file.path(root, paste0(mode, "-workers.log")))
on.exit(parallel::stopCluster(cluster), add = TRUE)
parallel::clusterExport(cluster, c("functionFile", "resultRoot", "mode", "force"),
  envir = environment())
parallel::clusterEvalQ(cluster, {
  suppressPackageStartupMessages({
    library(devtools); devtools::load_all("C:/package/saemix-copula", quiet = TRUE)
    library(rvinecopulib)
  })
  `%||%` <- function(x, y) if (is.null(x)) y else x
  source(functionFile); NULL
})
started <- Sys.time()
result <- parallel::parLapplyLB(cluster, pending, function(i) tryCatch({
  settings <- if (mode == "pilot")
    list(nSubjects = 180L, nVpc = 20000L, iterations = 1000L,
      screenDraws = 400L, screenRetry = 1000L, likelihoodDraws = 1200L) else
    list(nSubjects = 250L, nVpc = 50000L, iterations = 1500L,
      screenDraws = 700L, screenRetry = 2000L, likelihoodDraws = 2500L)
  do.call(combined_run_replicate, c(list(replicate = i,
    resultRoot = resultRoot, force = force), settings))
  list(ok = TRUE, replicate = i)
}, error = function(e) list(ok = FALSE, replicate = i,
  error = conditionMessage(e), call = paste(deparse(conditionCall(e)),
    collapse = " "))))
failed <- Filter(function(x) !isTRUE(x$ok), result)
saveRDS(result, file.path(root, paste0(mode, "-last-run-status.rds")))
cat("elapsed minutes:", round(difftime(Sys.time(), started, units = "mins"), 2),
  "failures:", length(failed), "\n")
if (length(failed)) { print(failed); quit(save = "no", status = 1L) }
