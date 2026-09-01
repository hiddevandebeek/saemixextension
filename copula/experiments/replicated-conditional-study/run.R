args <- commandArgs(trailingOnly = TRUE)
replicates <- if (length(args)) as.integer(args[1L]) else 50L
workers <- if (length(args) >= 2L) as.integer(args[2L]) else 8L
if (!is.finite(replicates) || replicates < 1L) stop("invalid replicate count")
if (!is.finite(workers) || workers < 1L) stop("invalid worker count")

repo <- normalizePath("C:/package/saemix-copula", winslash = "/")
root <- file.path(repo, "copula", "experiments", "replicated-conditional-study")
outputRoot <- file.path(root, "results")
dir.create(outputRoot, recursive = TRUE, showWarnings = FALSE)
functionFile <- file.path(root, "functions.R")
pending <- seq_len(replicates)[!file.exists(file.path(outputRoot,
  sprintf("replicate_%03d.rds", seq_len(replicates))))]
cat("Replicates:", replicates, "pending:", length(pending),
  "workers:", workers, "\n")
if (!length(pending)) quit(save = "no", status = 0L)

cluster <- parallel::makePSOCKcluster(min(workers, length(pending)),
  outfile = file.path(root, "worker.log"))
on.exit(parallel::stopCluster(cluster), add = TRUE)
parallel::clusterExport(cluster, c("repo", "functionFile", "outputRoot"),
  envir = environment())
parallel::clusterEvalQ(cluster, {
  suppressPackageStartupMessages({
    library(devtools)
    devtools::load_all(repo, quiet = TRUE)
    library(rvinecopulib)
  })
  `%||%` <- function(x, y) if (is.null(x)) y else x
  source(functionFile)
  NULL
})
started <- Sys.time()
result <- parallel::parLapplyLB(cluster, pending, function(i) {
  tryCatch(list(ok = TRUE, replicate = i,
      file = study_run_one(i, outputRoot)),
    error = function(e) list(ok = FALSE, replicate = i,
      error = conditionMessage(e)))
})
failed <- Filter(function(x) !isTRUE(x$ok), result)
saveRDS(result, file.path(root, "last_run_status.rds"))
cat("Completed in", round(difftime(Sys.time(), started, units = "mins"), 2),
  "minutes; failures:", length(failed), "\n")
if (length(failed)) {
  print(failed)
  quit(save = "no", status = 1L)
}
