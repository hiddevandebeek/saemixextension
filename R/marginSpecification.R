## Saying which margins to estimate, and how.
##
## Three things a user wants to be able to do, and which until now needed
## different entry points or were not possible at all:
##
##   fix a coordinate's family, because the analyst has decided it, while still
##     estimating its parameters;
##   fix it outright, parameters and all, because they are known from
##     elsewhere;
##   let a coordinate be chosen automatically;
##   let it be chosen automatically but only from a stated set -- "lognormal or
##     generalized gamma, nothing else" -- which is the common regulatory or
##     prior-knowledge case;
##
## and to mix all three in one model, coordinate by coordinate. That is what a
## specification is for. It is a named list whose entries are one of
##
##   a margin object            the *family* is fixed; its parameters are still
##                              estimated, and the values supplied are starting
##                              values. `copulaNaturalMarginLognormal(0.3)`
##                              means "lognormal, starting the scale at 0.3",
##                              not "lognormal with scale exactly 0.3"
##   copulaFixMargin(m)         the family *and* the parameters are fixed, so
##                              nothing about that coordinate is estimated
##   "auto"                     the default search for that coordinate
##   copulaAuto(...)            a restricted search, see below
##
## and anything not named keeps the default.

## A restricted automatic choice.
##
## `families` names the candidates to consider, from the support-based registry
## (`copulaNaturalMarginFamilies()` lists them). Giving none considers all of
## them for that coordinate's support, which is the default: a fitted margin
## then has a name a reader recognises -- lognormal, gamma, generalized gamma
## -- rather than a vector of tilt coefficients.
##
## `maxDegree` instead asks for the tilted ladder, searched from degree zero
## up to the value given. That family is strictly more flexible and strictly
## less interpretable: it contains the named shapes (a degree-2 tilt fitted to
## gamma draws recovers the gamma to a KL of 0.0003, against 0.023 for the
## nearest wrong family) but reports coefficients rather than a name, and at
## degree 2 those coefficients are individually near-unidentified even where
## the distribution itself is pinned down. Ask for it when a named family will
## not do.
##
## Giving both is an error: they are different searches and their selection
## criteria are not comparable.
copulaAuto <- function(families = NULL, maxDegree = NULL) {
  if (!is.null(families) && !is.null(maxDegree))
    stop("give either families or maxDegree, not both: they are different ",
      "searches and their selection criteria are not comparable")
  if (!is.null(families)) {
    families <- tolower(as.character(families))
    if (!length(families)) stop("families must name at least one candidate")
    known <- unique(unlist(lapply(c("real", "positive", "unit"),
      copulaNaturalMarginFamilies)))
    unknown <- setdiff(families, known)
    if (length(unknown))
      stop("unknown margin families: ", paste(unknown, collapse = ", "),
        ". Available: ", paste(sort(known), collapse = ", "))
  }
  if (!is.null(maxDegree)) {
    maxDegree <- as.integer(maxDegree)
    if (length(maxDegree) != 1L || is.na(maxDegree) || maxDegree < 0L ||
        maxDegree > 6L)
      stop("maxDegree must be a single integer between 0 and 6")
  }
  structure(list(families = families, maxDegree = maxDegree),
    class = "copulaAutoMargin")
}

copulaIsAutoMargin <- function(x) inherits(x, "copulaAutoMargin") ||
  (is.character(x) && length(x) == 1L && identical(tolower(x), "auto"))

## The specification itself. Entries may be named by coordinate name or given
## positionally; unnamed coordinates default to "auto".
copulaMarginSpec <- function(...) {
  entries <- list(...)
  if (length(entries) == 1L && is.list(entries[[1L]]) &&
      !copulaIsAutoMargin(entries[[1L]]) &&
      is.null(entries[[1L]]$log_density))
    entries <- entries[[1L]]
  for (k in seq_along(entries)) {
    entry <- entries[[k]]
    ok <- copulaIsAutoMargin(entry) ||
      (is.list(entry) && !is.null(entry$log_density))
    if (!ok)
      stop("entry ", k, " of the margin specification must be a margin ",
        "object, \"auto\", or copulaAuto(...)")
  }
  structure(entries, class = "copulaMarginSpec")
}

## Resolve a specification against the coordinates of a fit, giving one plan
## per coordinate. Coordinates absent from the specification are "auto".
copulaResolveMarginSpec <- function(spec, names, supports) {
  d <- length(names)
  supports <- rep_len(as.character(supports), d)
  if (is.null(spec)) spec <- copulaMarginSpec()
  if (!inherits(spec, "copulaMarginSpec"))
    stop("spec must be built by copulaMarginSpec()")
  given <- names(spec)
  if (is.null(given)) given <- rep("", length(spec))
  unknown <- setdiff(given[nzchar(given)], names)
  if (length(unknown))
    stop("the specification names coordinates that are not in the model: ",
      paste(unknown, collapse = ", "), ". Available: ",
      paste(names, collapse = ", "))
  plan <- vector("list", d)
  for (j in seq_len(d)) {
    entry <- if (names[j] %in% given) spec[[names[j]]] else
      if (j <= length(spec) && !nzchar(given[min(j, length(given))]))
        spec[[j]] else "auto"
    if (is.list(entry) && !is.null(entry$log_density)) {
      if (!identical(entry$support, supports[j]))
        stop("the margin given for ", names[j], " has support \"",
          entry$support, "\" but that coordinate is \"", supports[j], "\"")
      plan[[j]] <- list(kind = "fixed", margin = entry,
        estimated = sum(entry$free))
    } else if (inherits(entry, "copulaAutoMargin")) {
      if (!is.null(entry$families)) {
        allowed <- copulaNaturalMarginFamilies(supports[j])
        wrong <- setdiff(entry$families, allowed)
        if (length(wrong))
          stop("for ", names[j], ", which is \"", supports[j], "\", these ",
            "families are not available: ", paste(wrong, collapse = ", "),
            ". Available: ", paste(allowed, collapse = ", "))
        plan[[j]] <- list(kind = "families", families = entry$families)
      } else if (!is.null(entry$maxDegree)) {
        plan[[j]] <- list(kind = "ladder", maxDegree = entry$maxDegree)
      } else {
        plan[[j]] <- list(kind = "families",
          families = copulaNaturalMarginFamilies(supports[j]))
      }
    } else {
      ## The default is a named family, because a reader can act on "gamma" and
      ## cannot act on a tilt coefficient. The ladder is available by asking.
      plan[[j]] <- list(kind = "families",
        families = copulaNaturalMarginFamilies(supports[j]))
    }
  }
  stats::setNames(plan, names)
}

## A one-line human summary, so a user can see what a specification resolved to
## before spending an hour fitting it.
copulaDescribeMarginSpec <- function(plan) {
  vapply(seq_along(plan), function(j) {
    p <- plan[[j]]
    switch(p$kind,
      fixed = paste0(names(plan)[j], ": ", p$margin$name,
        if (p$estimated) paste0(", estimating ", p$estimated, " of ",
          length(p$margin$parameters), " parameters") else
          ", fully fixed"),
      families = paste0(names(plan)[j], ": named family, chosen among ",
        paste(p$families, collapse = ", ")),
      ladder = paste0(names(plan)[j], ": tilted ladder, degree 0 to ",
        p$maxDegree),
      paste0(names(plan)[j], ": ?"))
  }, character(1))
}
