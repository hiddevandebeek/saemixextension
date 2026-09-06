## A dangling reference to a variable that no longer exists -- a Cholesky
## factor removed when a kernel stopped needing it, but still named in the
## value the kernel returned -- broke two fits, and only running those fits
## found it. codetools finds that class of mistake statically, so this sweeps
## every function in the namespace for names it uses but never defines.
##
## The allowed names are the ones that are genuinely resolved elsewhere:
## constants R supplies, the variables that local() closures capture, and the
## data-masking pronouns used inside plotting code.

suppressPackageStartupMessages({
  library(devtools); load_all("C:/package/saemix-copula", quiet = TRUE)
})

allowed <- c(".Machine", ".Random.seed", ".data", ".x",
  "cache", "last", "value", "both", "free", "onlyLower", "onlyUpper")
namespace <- asNamespace("saemix")
known <- c(ls(namespace, all.names = TRUE), ls(baseenv()), allowed,
  unlist(lapply(c("stats", "utils", "methods", "graphics", "grDevices"),
    function(package) ls(asNamespace(package)))))

suspects <- character(0)
for (name in ls(namespace, all.names = TRUE)) {
  object <- get(name, envir = namespace)
  if (!is.function(object)) next
  globals <- try(codetools::findGlobals(object, merge = FALSE)$variables,
    silent = TRUE)
  if (inherits(globals, "try-error")) next
  undefined <- setdiff(globals, known)
  if (length(undefined))
    suspects <- c(suspects, sprintf("%s uses %s", name,
      paste(undefined, collapse = ", ")))
}

cat(sprintf("%-54s %s %d\n", "functions referencing an undefined name",
  if (!length(suspects)) "PASS" else "**FAIL**", length(suspects)))
if (length(suspects)) {
  writeLines(suspects)
  quit(status = 1L)
}
cat("\nno undefined names\n")
