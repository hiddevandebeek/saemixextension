## Convenience entry point for the validated 10-dataset true-Gamma study.
## Use --force to refit; otherwise existing replicate files are summarized.

root <- "C:/package/saemix-copula/copula/experiments/non-gaussian-eta-demonstration"
force <- "--force" %in% commandArgs(trailingOnly = TRUE)
resultFiles <- list.files(file.path(root, "replicates"),
  pattern = "^replicate_[0-9]{3}[.]rds$", full.names = TRUE)
if (force || length(resultFiles) != 10L) {
  command <- c(file.path(root, "run_replicates.R"), "10", "5")
  if (force) command <- c(command, "--force")
  status <- system2(file.path(R.home("bin"), "Rscript"),
    c("--vanilla", command))
  if (!identical(status, 0L)) stop("replicated fitting failed")
}
status <- system2(file.path(R.home("bin"), "Rscript"),
  c("--vanilla", file.path(root, "summarize_replicates.R")))
if (!identical(status, 0L)) stop("replicate summarization failed")
