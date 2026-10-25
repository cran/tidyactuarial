#' Specify an empirical claim frequency distribution
#'
#' @md
#' @description
#' Constructs the exact discrete distribution implied by observed claim counts
#' or a weighted count table. No parametric distribution is fitted, no smoothing
#' is applied, and unobserved counts receive probability zero.
#'
#' @param data A nonempty numeric vector of observed nonnegative integer counts,
#'   or a data frame containing those counts.
#' @param count For a data frame, an unquoted column containing the observed
#'   counts. Required for data-frame input and omitted for vector input.
#' @param weight Optional nonnegative observation weights. For a data frame,
#'   provide an unquoted column; for a vector input, provide a numeric vector
#'   of the same length. When omitted, every row has weight one. Weights in a
#'   summarized frequency table may represent the number of policies having
#'   each claim count; they are **not** exposure offsets for fitting a rate.
#' @param label Optional nonempty character scalar identifying the distribution.
#'
#' @details
#' If the observed counts are \eqn{n_1,\ldots,n_m} and their nonnegative weights
#' are \eqn{w_1,\ldots,w_m}, the empirical probability mass function is
#' \deqn{\widehat p_k = \frac{\sum_{i=1}^m w_i\,I(n_i=k)}{\sum_{i=1}^m w_i}.}
#'
#' At least one weight must be positive. Rows with zero weight do not contribute
#' to the empirical support. Missing, negative, fractional, or nonfinite claim
#' counts are rejected rather than silently discarded. Counts must not exceed
#' \eqn{2^{53}}, the largest consecutive integer representable exactly by an
#' IEEE-754 double. Weights must be finite,
#' nonnegative, and have a strictly positive finite total. Counts are kept in
#' increasing order; no gaps in the observed set are filled with artificial
#' positive probabilities.
#'
#' This is the empirical distribution of the **claim count per observational
#' unit** (such as a policy-year), not a table of individual claims. A vector
#' of all ones describes a degenerate count distribution at one, not the total
#' number of claims in a portfolio.
#'
#' The result has the same lightweight S3 contract as [frequency_spec()]. Its
#' `dist` field is `"empirical"`; `params` holds the finite count values and
#' their exact empirical probabilities, and `data` contains the aggregated
#' observed weights and probabilities. Only this compact frequency table is
#' retained, not the original individual records.
#'
#' The internal frequency registry provides the PMF, CDF, quantile function,
#' RNG, raw and factorial moments, and probability generating function for
#' this model. A general finite empirical frequency distribution is not
#' assumed to belong to the Panjer \eqn{(a,b,0)} class; its `panjer` operation
#' is therefore unavailable. Other methods can use the same specification
#' without silently fitting or approximating its probabilities.
#'
#' @return An object of classes `"tidyact_frequency_spec"` and
#'   `"tidyact_distribution_spec"`, with the ordered fields `schema_version`,
#'   `dist`, `params`, `support`, `data`, `custom`, `label`, and `metadata`.
#'   The `data` field is a [tibble::tibble()] with columns `count`, `weight`,
#'   and `probability`, one row per count with positive empirical mass.
#'
#' @seealso [frequency_spec()], [distribution_moments()]
#'
#' @examples
#' # One observation per policy-year: 0, 1, or 2 claims.
#' counts <- c(0, 0, 0, 1, 1, 2)
#' n <- empirical_frequency(counts)
#' n
#' n$data
#' distribution_moments(n, order = 1:3)
#' distribution_moments(n, order = 1:3, type = "factorial")
#'
#' # The same empirical distribution as a summarized portfolio table.
#' portfolio <- tibble::tibble(
#'   claims = c(0, 1, 2),
#'   policies = c(3, 2, 1)
#' )
#' portfolio |>
#'   empirical_frequency(count = claims, weight = policies)
#'
#' # Unweighted data-frame input and a descriptive label.
#' observations <- tibble::tibble(claim_count = c(0, 1, 0, 2, 1))
#' observations |>
#'   empirical_frequency(count = claim_count, label = "Observed policy counts")
#'
#' @export
empirical_frequency <- function(data, count = NULL, weight = NULL, label = NULL) {
  .frequency_validate_label(label)

  count_expr <- rlang::enquo(count)
  weight_expr <- rlang::enquo(weight)

  if (is.data.frame(data)) {
    if (rlang::quo_is_null(count_expr)) {
      stop("`count` must identify a column when `data` is a data frame.",
           call. = FALSE)
    }
    counts <- rlang::eval_tidy(count_expr, data = data)
    if (length(counts) != nrow(data)) {
      stop("`count` must supply exactly one observation per data-frame row.",
           call. = FALSE)
    }
    weights <- if (rlang::quo_is_null(weight_expr)) {
      rep(1, length(counts))
    } else {
      rlang::eval_tidy(weight_expr, data = data)
    }
  } else {
    if (!rlang::quo_is_null(count_expr)) {
      stop("`count` must be omitted when `data` is a numeric vector.",
           call. = FALSE)
    }
    counts <- data
    weights <- if (rlang::quo_is_null(weight_expr)) {
      rep(1, length(counts))
    } else {
      rlang::eval_tidy(weight_expr)
    }
  }

  if (!is.numeric(counts) || !is.null(dim(counts)) || length(counts) == 0L ||
      anyNA(counts) || any(!is.finite(counts)) || any(counts < 0) ||
      any(counts != floor(counts)) || any(counts > 2^53)) {
    stop("`count` observations must be a nonempty numeric vector of finite nonnegative integers.",
         call. = FALSE)
  }

  if (!is.numeric(weights) || !is.null(dim(weights)) ||
      length(weights) != length(counts) || anyNA(weights) ||
      any(!is.finite(weights)) || any(weights < 0)) {
    stop("`weight` must contain one finite nonnegative number per count observation.",
         call. = FALSE)
  }

  positive <- weights > 0
  if (!any(positive)) {
    stop("At least one observation must have positive weight.",
         call. = FALSE)
  }

  counts <- as.numeric(counts[positive])
  weights <- as.numeric(weights[positive])
  values <- sort(unique(counts))

  # Numeric group indices keep the sorting of large integer-valued counts;
  # no string conversion or expansion of summarized observations is needed.
  group <- match(counts, values)
  grouped <- rowsum(weights, group = group, reorder = TRUE)
  mass <- as.numeric(grouped[as.character(seq_along(values)), 1L])

  total_weight <- sum(mass)
  if (!is.finite(total_weight) || total_weight <= 0) {
    stop("The total observation weight must be positive and finite.",
         call. = FALSE)
  }

  probs <- mass / total_weight
  if (any(probs <= 0)) {
    stop("The weight scale causes numerical underflow in empirical probabilities.",
         call. = FALSE)
  }
  if (any(!is.finite(probs)) || abs(sum(probs) - 1) > 1e-12) {
    stop("Empirical probabilities could not be normalized accurately.",
         call. = FALSE)
  }

  frequency_table <- tibble::tibble(
    count = values,
    weight = mass,
    probability = probs
  )

  structure(
    list(
      schema_version = 1L,
      dist = "empirical",
      params = list(values = values, probs = probs),
      support = list(
        type = "count",
        lower = min(values),
        upper = max(values),
        step = 1
      ),
      data = frequency_table,
      custom = NULL,
      label = label,
      metadata = list()
    ),
    class = c("tidyact_frequency_spec", "tidyact_distribution_spec")
  )
}


# Internal registry entry --------------------------------------------------
# Integration: add `empirical = .empirical_frequency_registry_entry(),`
# to the list returned by .frequency_registry() in frequency-spec.R.
# This entry is not a separate public constructor.

.empirical_frequency_registry_entry <- function() {
  list(
    validate = function(params, custom = NULL) {
      stop("Use `empirical_frequency()` to specify an empirical count distribution.",
           call. = FALSE)
    },
    support = function(params, custom = NULL) {
      list(type = "count", lower = min(params$values),
           upper = max(params$values), step = 1)
    },
    pmf = function(x, params, log = FALSE) {
      index <- match(x, params$values)
      out <- numeric(length(x))
      found <- !is.na(index)
      out[found] <- params$probs[index[found]]
      out[is.na(x)] <- NA_real_
      if (log) log(out) else out
    },
    cdf = function(q, params, lower.tail = TRUE, log.p = FALSE) {
      out <- rep(NA_real_, length(q))
      ok <- !is.na(q)
      if (any(ok)) {
        cumulative <- cumsum(params$probs)
        cumulative[length(cumulative)] <- 1
        indices <- findInterval(q[ok], params$values)
        result <- numeric(length(indices))
        found <- indices > 0L
        result[found] <- cumulative[indices[found]]
        out[ok] <- result
      }
      if (!lower.tail) out <- 1 - out
      if (log.p) log(out) else out
    },
    quantile = function(p, params, lower.tail = TRUE, log.p = FALSE) {
      if (!is.numeric(p)) {
        stop("`p` must be numeric.", call. = FALSE)
      }

      invalid <- !is.na(p) & if (log.p) p > 0 else (p < 0 | p > 1)
      pp <- if (log.p) exp(p) else p
      pp[invalid] <- NaN
      if (!lower.tail) pp <- 1 - pp

      cumulative <- cumsum(params$probs)
      cumulative[length(cumulative)] <- 1

      vapply(pp, function(target) {
        if (is.na(target)) return(NaN)
        if (target <= 0) return(params$values[1L])
        if (target >= 1) return(params$values[length(params$values)])

        tol <- 8 * .Machine$double.eps * max(1, abs(target))
        index <- which(cumulative >= target - tol)[1L]

        if (is.na(index)) {
          params$values[length(params$values)]
        } else {
          params$values[index]
        }
      }, numeric(1))
    },
    rng = function(n, params) {
      if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
          !is.finite(n) || n < 0 || n != floor(n) ||
          n > .Machine$integer.max) {
        stop("`n` must be a nonnegative integer scalar.", call. = FALSE)
      }
      if (n == 0) return(numeric())
      index <- sample.int(length(params$values), size = as.integer(n),
                          replace = TRUE, prob = params$probs)
      params$values[index]
    },
    raw_moment = function(order, params) {
      .empirical_frequency_check_order(order)
      sum(params$probs * params$values^order)
    },
    factorial_moment = function(order, params) {
      .empirical_frequency_check_order(order)
      if (order == 0) return(1)
      keep <- params$values >= order
      if (!any(keep)) return(0)
      values <- params$values[keep]
      probs <- params$probs[keep]
      sum(probs * exp(lgamma(values + 1) - lgamma(values - order + 1)))
    },
    pgf = function(z, params) {
      if (!is.numeric(z) && !is.complex(z)) {
        stop("`z` must be numeric or complex.", call. = FALSE)
      }
      vapply(
        z,
        function(t) sum(params$probs * t^params$values),
        if (is.complex(z)) 0 + 0i else 0
      )
    },
    panjer = NULL
  )
}


.empirical_frequency_check_order <- function(order) {
  if (!is.numeric(order) || length(order) != 1L || is.na(order) ||
      !is.finite(order) || order < 0 || order != floor(order) ||
      order > .Machine$integer.max) {
    stop("`order` must be a nonnegative integer scalar.", call. = FALSE)
  }
  invisible(TRUE)
}
