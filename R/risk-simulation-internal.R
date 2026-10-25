# Internal utilities shared by risk simulation and risk measures.
# These are intentionally not exported.

.risk_check_n_sim <- function(n_sim) {
  if (!is.numeric(n_sim) || length(n_sim) != 1L || is.na(n_sim) ||
      !is.finite(n_sim) || n_sim < 1 || n_sim != floor(n_sim) ||
      n_sim > .Machine$integer.max) {
    stop("`n_sim` must be a positive integer within R's integer range.",
         call. = FALSE)
  }
  as.integer(n_sim)
}

.risk_check_seed <- function(seed) {
  if (is.null(seed)) return(invisible(NULL))
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed != floor(seed) || seed < 0 ||
      seed > .Machine$integer.max) {
    stop("`seed` must be NULL or a nonnegative integer seed.", call. = FALSE)
  }
  invisible(NULL)
}

.risk_with_seed <- function(seed, expr) {
  .risk_check_seed(seed)
  if (is.null(seed)) return(force(expr))
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) previous <- get(".Random.seed", envir = .GlobalEnv,
                                inherits = FALSE)
  on.exit({
    if (existed) {
      assign(".Random.seed", previous, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv,
                      inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))
  force(expr)
}

.risk_simulation_callback <- function(x, kind) {
  registry <- if (identical(kind, "frequency")) .frequency_registry() else
    .severity_registry()
  if (identical(x$dist, "custom")) {
    fn <- x$custom$rng
  } else {
    entry <- registry[[x$dist]]
    if (is.null(entry)) {
      stop("Unsupported distribution in risk simulation: ", x$dist, ".",
           call. = FALSE)
    }
    fn <- entry$rng
  }
  if (!is.function(fn)) {
    stop("Distribution '", x$dist,
         "' has no RNG provider; simulation is unavailable.", call. = FALSE)
  }
  fn
}

.risk_draw <- function(x, n, kind) {
  if (n == 0L) return(numeric())
  fn <- .risk_simulation_callback(x, kind)
  out <- fn(n, x$params)
  if (!is.numeric(out) || length(out) != n || anyNA(out) ||
      any(!is.finite(out)) || any(out < 0)) {
    stop("The ", kind, " RNG returned invalid values.", call. = FALSE)
  }
  if (identical(kind, "frequency") && any(out != floor(out))) {
    stop("The frequency RNG must return nonnegative integer counts.",
         call. = FALSE)
  }
  as.numeric(out)
}

.risk_check_max_claims <- function(max_claims) {
  if (!is.numeric(max_claims) || length(max_claims) != 1L ||
      is.na(max_claims) || !is.finite(max_claims) ||
      max_claims < 1 || max_claims != floor(max_claims) ||
      max_claims > .Machine$integer.max) {
    stop("`max_claims` must be a positive integer within R's integer range.",
         call. = FALSE)
  }
  as.integer(max_claims)
}
