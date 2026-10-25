# Parallel scenarios share a mortality table and the existing scalar engines.
# Only size-one arguments recycle; NULL keeps its function-specific meaning.
.single_life_scenarios <- function(fun, lt, args, tidy) {
  sizes <- vapply(args, function(x) if (is.null(x)) 1L else length(x), integer(1L))
  size <- max(sizes)
  if (any(sizes == 0L) || any(!sizes %in% c(1L, size))) {
    stop("Scenario arguments must have length one or a common positive length.",
         call. = FALSE)
  }
  if (size == 1L) return(NULL)
  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }
  rows <- lapply(seq_len(size), function(j) {
    values <- lapply(args, function(x) {
      if (is.null(x) || length(x) == 1L) x else x[j]
    })
    tryCatch(do.call(fun, c(list(lt = lt), values, list(tidy = tidy))),
             error = function(e) {
               stop("Scenario ", j, ": ", conditionMessage(e), call. = FALSE)
             })
  })
  if (tidy) dplyr::bind_rows(rows) else unlist(rows, use.names = FALSE)
}

# Ratio first preserves ordinary accuracy near one. Log counts retain range
# when the ratio underflows, e.g. l(x+t)=1e-100 and l(x)=1e300.
.single_life_log_ratio <- function(numerator, denominator) {
  ratio <- numerator / denominator
  out <- log(ratio)
  recover <- !is.na(ratio) & ratio == 0 & numerator > 0
  out[recover] <- log(numerator[recover]) - log(denominator)
  out
}
