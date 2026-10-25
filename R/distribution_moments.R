#' Calculate moments of a frequency or severity distribution
#'
#' @md
#' @description
#' Calculates raw, central, factorial, or cumulant moments from an existing
#' frequency or severity specification. The calculation is deterministic:
#' it does not fit a model, generate observations, or approximate moments
#' using simulation or numerical integration.
#'
#' @param x A `tidyact_frequency_spec` or `tidyact_severity_spec` object.
#' @param order Nonempty numeric vector of nonnegative integer moment orders.
#'   For cumulants, all orders must be positive.
#' @param type One of `"raw"` (default), `"central"`, `"factorial"`, or
#'   `"cumulant"`. Factorial moments are available only for frequency
#'   specifications.
#'
#' @details
#' For a random variable \eqn{Z}, the raw and central moments are
#' \deqn{m_r = E[Z^r], \qquad \mu_r = E[(Z-E[Z])^r].}
#'
#' For a nonnegative integer-valued frequency \eqn{N}, the factorial moment is
#' \deqn{E[(N)_r], \qquad (N)_r=N(N-1)\cdots(N-r+1).}
#'
#' Cumulants are calculated algebraically from raw moments. In particular,
#' \eqn{\kappa_1=E[Z]}, \eqn{\kappa_2=\mathrm{Var}(Z)}, and
#' \eqn{\kappa_3=E[(Z-E[Z])^3]} when the relevant moments are finite.
#'
#' Moment order zero is one for raw, central, and factorial moments.
#' Cumulants start at order one.
#'
#' All built-in distributions use their registered analytical raw or
#' factorial moment functions. A custom specification must explicitly supply
#' `raw_moment` or (for a frequency) `factorial_moment` to provide a moment
#' calculation; supplying only a PMF or PDF does not trigger numerical
#' integration or summation. Raw and factorial frequency moments can be
#' converted into one another using Stirling numbers when the necessary
#' source moments are available.
#'
#' The `status` column distinguishes:
#'
#' - `"finite"`: a finite numerical result;
#' - `"infinite"`: a mathematically divergent moment;
#' - `"undefined"`: the moment is not defined or the specification does not
#'   provide sufficient information to evaluate it;
#' - `"overflow"`: a theoretically finite moment or an intermediate
#'   algebraic expression exceeds the available floating-point range.
#'
#' For built-in Pareto II, Burr XII, and log-logistic severities, divergence
#' is determined from the distribution's analytical moment-existence
#' condition. For custom distributions, an explicitly returned positive
#' `Inf` is interpreted as mathematical divergence; custom moment functions
#' are responsible for this distinction.
#'
#' A central moment of positive order requires a finite mean. When the mean
#' exists but a higher-order moment diverges, the corresponding central
#' moment is infinite. Higher-order cumulants are reported as undefined if
#' the moments needed for their definition are infinite. The first cumulant
#' follows the first raw moment.
#'
#' Central moments and cumulants derived from raw moments may suffer from
#' cancellation for extremely large parameter values. Such arithmetic is
#' not silently replaced by a simulated or approximated answer.
#'
#' @return A [tibble::tibble()] with one row per requested order and columns:
#'   `type` (character), `order` (integer), `value` (double), and `status`
#'   (character). Input order and repeated orders are preserved. Undefined
#'   values are returned as `NA_real_`.
#'
#' @references
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A., and
#' Nesbitt, C. J. (1997). *Actuarial Mathematics*, 2nd ed.
#' Society of Actuaries.
#'
#' @seealso [frequency_spec()], [severity_spec()]
#'
#' @examples
#' # Raw and central moments of a claim count distribution.
#' counts <- frequency_spec("poisson", lambda = 4)
#' distribution_moments(counts)
#' distribution_moments(counts, order = 1:4, type = "central")
#' distribution_moments(counts, order = 0:4, type = "factorial")
#' distribution_moments(counts, order = 1:4, type = "cumulant")
#'
#' # Raw moments of a nonnegative claim amount.
#' loss <- severity_spec("gamma", shape = 2.5, rate = 0.001)
#' distribution_moments(loss, order = 1:4)
#'
#' # The mean exists but the variance does not.
#' heavy_tail <- severity_spec("pareto", shape = 1.5, scale = 2000)
#' distribution_moments(heavy_tail, order = 1:4)
#' distribution_moments(heavy_tail, order = 1:4, type = "central")
#'
#' # A nonrandom payment has zero central moments of positive order.
#' fixed <- severity_spec("degenerate", value = 1000)
#' distribution_moments(fixed, order = 1:4, type = "central")
#'
#' @export
 distribution_moments <- function(x, order = 1:4, type = "raw") {
  is_frequency <- inherits(x, "tidyact_frequency_spec")
  is_severity <- inherits(x, "tidyact_severity_spec")

  if (!is_frequency && !is_severity) {
    stop(
      "`x` must be a frequency or severity distribution specification.",
      call. = FALSE
    )
  }

  if (!is.numeric(order) || length(order) == 0L || anyNA(order) ||
      any(!is.finite(order)) || any(order < 0) ||
      any(order != floor(order)) ||
      any(order > .Machine$integer.max)) {
    stop("`order` must be a nonempty vector of nonnegative integers.",
         call. = FALSE)
  }

  allowed <- c("raw", "central", "factorial", "cumulant")
  if (!is.character(type) || length(type) != 1L ||
      is.na(type) || !type %in% allowed) {
    stop("`type` must be one of: raw, central, factorial, cumulant.",
         call. = FALSE)
  }

  if (type == "factorial" && !is_frequency) {
    stop("Factorial moments require a frequency specification.",
         call. = FALSE)
  }

  if (type == "cumulant" && any(order == 0)) {
    stop("Cumulant orders must be positive integers.", call. = FALSE)
  }

  order <- as.integer(order)
  registry <- if (is_frequency) .frequency_registry() else .severity_registry()
  entry <- registry[[x$dist]]

  if (is.null(entry)) {
    stop("Unknown distribution identifier in `x`: ", x$dist, ".",
         call. = FALSE)
  }

  calculate <- switch(
    type,
    raw = .distribution_raw_moment,
    central = .distribution_central_moment,
    factorial = .distribution_factorial_moment,
    cumulant = .distribution_cumulant
  )

  results <- lapply(order, function(k) calculate(x, k, entry))

  tibble::tibble(
    type = rep(type, length(order)),
    order = order,
    value = vapply(results, function(z) z$value, numeric(1)),
    status = vapply(results, function(z) z$status, character(1))
  )
}


# Internal moment results --------------------------------------------------

.distribution_moment_result <- function(value, status = "finite") {
  list(value = as.numeric(value), status = status)
}

.distribution_moment_unavailable <- function() {
  .distribution_moment_result(NA_real_, "undefined")
}

.distribution_moment_overflow <- function() {
  .distribution_moment_result(Inf, "overflow")
}

.distribution_moment_callback <- function(x, entry, operation) {
  if (identical(x$dist, "custom")) {
    return(x$custom[[operation]])
  }
  entry[[operation]]
}

# The three heavy-tailed built-in families have an exact moment-existence
# condition. This is deliberately checked before evaluating an expression
# that could also produce Inf from numerical overflow.
.distribution_raw_diverges <- function(x, order) {
  if (!inherits(x, "tidyact_severity_spec") || order == 0L) {
    return(FALSE)
  }

  switch(
    x$dist,
    pareto = order >= x$params$shape,
    burr = order / x$params$shape1 >= x$params$shape2,
    loglogistic = order >= x$params$shape,
    FALSE
  )
}

.distribution_direct_moment <- function(x, order, entry, operation) {
  callback <- .distribution_moment_callback(x, entry, operation)
  if (!is.function(callback)) {
    return(.distribution_moment_unavailable())
  }

  if (operation == "raw_moment" &&
      .distribution_raw_diverges(x, order)) {
    return(.distribution_moment_result(Inf, "infinite"))
  }

  value <- callback(order, x$params)
  if (!is.numeric(value) || length(value) != 1L) {
    stop("`", operation, "` must return one numeric moment value.",
         call. = FALSE)
  }
  if (is.na(value)) {
    return(.distribution_moment_unavailable())
  }
  if (value < 0) {
    stop("`", operation, "` returned a negative nonnegative moment.",
         call. = FALSE)
  }
  if (is.infinite(value)) {
    status <- if (identical(x$dist, "custom")) "infinite" else "overflow"
    return(.distribution_moment_result(Inf, status))
  }
  .distribution_moment_result(value)
}


# Raw and factorial frequency moments -------------------------------------

# Combine an explicit vector of polynomial coefficients with exact source
# moments. No summation over an inferred infinite support is attempted.
.distribution_polynomial_moment <- function(coefficients, source) {
  indices <- which(coefficients != 0)
  if (!length(indices)) return(.distribution_moment_result(0))

  parts <- lapply(indices - 1L, source)
  statuses <- vapply(parts, function(z) z$status, character(1))

  if (any(statuses == "undefined")) return(.distribution_moment_unavailable())
  if (any(statuses == "overflow")) return(.distribution_moment_overflow())
  if (any(statuses == "infinite")) {
    # This helper is used only for raw <-> factorial transformations of a
    # count distribution; both moments of the same positive order have the
    # same existence condition.
    return(.distribution_moment_result(Inf, "infinite"))
  }

  values <- vapply(parts, function(z) z$value, numeric(1))
  answer <- sum(coefficients[indices] * values)
  if (!is.finite(answer)) return(.distribution_moment_overflow())
  if (answer < 0 &&
      abs(answer) <= 1e-12 * max(1, sum(abs(coefficients[indices] * values)))) {
    answer <- 0
  }
  if (answer < 0) {
    return(.distribution_moment_unavailable())
  }
  .distribution_moment_result(answer)
}

.distribution_raw_moment <- function(x, order, entry) {
  if (order == 0L) return(.distribution_moment_result(1))

  raw_fn <- .distribution_moment_callback(x, entry, "raw_moment")
  if (is.function(raw_fn)) {
    return(.distribution_direct_moment(x, order, entry, "raw_moment"))
  }

  if (!inherits(x, "tidyact_frequency_spec") ||
      !is.function(.distribution_moment_callback(x, entry, "factorial_moment"))) {
    return(.distribution_moment_unavailable())
  }

  coefficients <- numeric(order + 1L)
  for (k in seq_len(order)) {
    coefficients[k + 1L] <- .frequency_stirling2(order, k)
  }

  .distribution_polynomial_moment(
    coefficients,
    function(k) .distribution_direct_moment(x, k, entry, "factorial_moment")
  )
}

.distribution_factorial_moment <- function(x, order, entry) {
  if (order == 0L) return(.distribution_moment_result(1))

  factorial_fn <- .distribution_moment_callback(x, entry, "factorial_moment")
  if (is.function(factorial_fn)) {
    return(.distribution_direct_moment(x, order, entry, "factorial_moment"))
  }

  if (!is.function(.distribution_moment_callback(x, entry, "raw_moment"))) {
    return(.distribution_moment_unavailable())
  }

  # Coefficients of (z)_order = z(z-1)...(z-order+1), ordered by z^k.
  coefficients <- 1
  for (k in seq_len(order)) {
    coefficients <- c(0, coefficients) - (k - 1) * c(coefficients, 0)
  }

  .distribution_polynomial_moment(
    coefficients,
    function(k) .distribution_direct_moment(x, k, entry, "raw_moment")
  )
}


# Central moments ----------------------------------------------------------

.distribution_central_moment <- function(x, order, entry) {
  if (order == 0L) return(.distribution_moment_result(1))

  mean <- .distribution_raw_moment(x, 1L, entry)
  if (mean$status != "finite") {
    if (mean$status == "overflow") return(.distribution_moment_overflow())
    return(.distribution_moment_unavailable())
  }
  if (order == 1L) return(.distribution_moment_result(0))

  # Exact identities avoid subtractive cancellation for point masses.
  if (inherits(x, "tidyact_severity_spec") &&
      identical(x$dist, "degenerate")) {
    return(.distribution_moment_result(0))
  }

  # For finite discrete severity distributions, direct centering is more
  # stable than reconstructing central moments from raw powers.
  if (inherits(x, "tidyact_severity_spec") &&
      identical(x$dist, "discrete")) {
    centered <- sum(x$params$probs * (x$params$values - mean$value)^order)
    if (!is.finite(centered)) return(.distribution_moment_overflow())
    return(.distribution_moment_result(centered))
  }

  higher <- .distribution_raw_moment(x, order, entry)
  if (higher$status == "infinite") {
    return(.distribution_moment_result(Inf, "infinite"))
  }
  if (higher$status == "overflow") return(.distribution_moment_overflow())
  if (higher$status == "undefined") return(.distribution_moment_unavailable())

  raw <- lapply(0:order, function(k) .distribution_raw_moment(x, k, entry))
  statuses <- vapply(raw, function(z) z$status, character(1))
  if (any(statuses == "undefined")) return(.distribution_moment_unavailable())
  if (any(statuses == "overflow")) return(.distribution_moment_overflow())
  if (any(statuses == "infinite")) return(.distribution_moment_unavailable())

  moments <- vapply(raw, function(z) z$value, numeric(1))
  k <- 0:order
  terms <- choose(order, k) * (-mean$value)^(order - k) * moments
  result <- sum(terms)
  if (!is.finite(result)) return(.distribution_moment_overflow())

  if (order %% 2L == 0L && result < 0) {
    if (abs(result) <= 1e-12 * max(1, sum(abs(terms)))) {
      result <- 0
    } else {
      return(.distribution_moment_unavailable())
    }
  }
  .distribution_moment_result(result)
}


# Cumulants ----------------------------------------------------------------

.distribution_cumulant <- function(x, order, entry) {
  if (order == 1L) return(.distribution_raw_moment(x, 1L, entry))

  if (inherits(x, "tidyact_severity_spec") &&
      identical(x$dist, "degenerate")) {
    return(.distribution_moment_result(0))
  }

  raw <- lapply(0:order, function(k) .distribution_raw_moment(x, k, entry))
  statuses <- vapply(raw, function(z) z$status, character(1))
  if (any(statuses == "undefined") || any(statuses == "infinite")) {
    return(.distribution_moment_unavailable())
  }
  if (any(statuses == "overflow")) return(.distribution_moment_overflow())

  moments <- vapply(raw, function(z) z$value, numeric(1))
  cumulants <- numeric(order)

  for (n in seq_len(order)) {
    kappa <- moments[n + 1L]
    if (n > 1L) {
      for (j in seq_len(n - 1L)) {
        kappa <- kappa - choose(n - 1L, j - 1L) *
          cumulants[j] * moments[n - j + 1L]
      }
    }
    if (!is.finite(kappa)) return(.distribution_moment_overflow())
    cumulants[n] <- kappa
  }

  .distribution_moment_result(cumulants[order])
}
