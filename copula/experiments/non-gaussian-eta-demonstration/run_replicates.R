suppressPackageStartupMessages(library(parallel))

args <- commandArgs(trailingOnly = TRUE)
replicates <- if (length(args)) as.integer(args[1L]) else 10L
workers <- if (length(args) >= 2L) as.integer(args[2L]) else 5L
force <- "--force" %in% args
if (!is.finite(replicates) || replicates < 1L ||
    !is.finite(workers) || workers < 1L) stop("invalid replicate controls")

repo <- "C:/package/saemix-copula"
root <- file.path(repo, "copula", "experiments",
  "non-gaussian-eta-demonstration")
resultRoot <- file.path(root, "replicates")
dir.create(resultRoot, recursive = TRUE, showWarnings = FALSE)

cluster <- makeCluster(min(workers, replicates), outfile = "")
on.exit(try(stopCluster(cluster), silent = TRUE), add = TRUE)
clusterExport(cluster, c("repo", "root", "resultRoot", "force"))
clusterEvalQ(cluster, {
  suppressPackageStartupMessages({
    library(devtools)
    load_all(repo, quiet = TRUE)
  })
  source(file.path(root, "functions.R"))
  NULL
})
result <- parLapplyLB(cluster, seq_len(replicates), function(replicate)
  ng_run_replicate(replicate, resultRoot, force = force))
stopifnot(length(result) == replicates,
  all(vapply(result, `[[`, integer(1), "schema") == 1L))
cat("Completed", replicates, "independent paired datasets in", resultRoot, "\n")
