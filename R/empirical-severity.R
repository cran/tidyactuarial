#' Specify an empirical claim severity distribution
#'
#' @md
#' @description
#' Constructs the discrete severity distribution implied by observed individual
#' claim amounts or a weighted table of claim amounts. No parametric family is
#' fitted, observations are not rounded, and no smoothing or discretization
#' is performed.
#'
#' @param data A nonempty numeric vector of nonnegative observed claim amounts,
#'   or a data frame containing those amounts.
#' @param amount For data-frame input, an unquoted column containing individual
#'   claim amounts. Required for data frames; omit for numeric vectors.
#' @param weight Optional nonnegative observation weights. For data frames,
#'   supply an unquoted column; for vectors, supply a numeric vector of equal
#'   length. If omitted, all observations have weight one. In a summarized
#'   table, weights can represent the multiplicities of each observed amount.
#'   Weights are not exposure offsets or inputs to a fitted model.
#' @param label Optional nonempty character scalar identifying the distribution.
#'
#' @details
#' Given observed amounts \eqn{x_1,\ldots,x_n} and nonnegative weights
#' \eqn{w_1,\ldots,w_n}, the exact empirical probability assigned to amount
#' \eqn{x} is
#' \deqn{\widehat{p}(x) =
#'   \frac{\sum_{i=1}^n w_i I(x_i=x)}{\sum_{i=1}^n w_i}.}
#'
#' Only strictly positive-weight observations contribute to the distribution;
#' at least one positive weight is required. Repeated amounts are combined and
#' stored in increasing order. Every unobserved amount has probability zero.
#' Amounts must be finite, nonnegative numeric values. An observed amount of
#' zero is valid and produces a genuine probability mass at zero; it is not
#' removed or converted into a positive-payment severity.
#'
#' The result is a *finite discrete* severity specification, even when the
#' observed amounts are decimal or come from an underlying continuous loss
#' process. It is **not** a continuous density estimate and has no PDF.
#' Values are compared exactly as represented numerically: close but unequal
#' amounts are not merged. The function also does not infer a loss/payment
#' basis, deductibles, limits, or monetary rounding rules.
#'
#' The S3 contract matches [severity_spec()]: `dist = "empirical"`, with
#' `params = list(values = ..., probs = ...)` and a compact `data` table.
#' The original individual observations are not retained. The internal
#' severity registry makes the PMF, CDF, quantile, RNG, raw moments,
#' limited moments, MGF, and Laplace transform available. No simulation
#' occurs during construction.
#'
#' For these finite observations all positive integer-order raw moments exist
#' mathematically, although large powers can exceed floating-point range.
#'
#' @return An object with classes `"tidyact_severity_spec"` and
#'   `"tidyact_distribution_spec"`, containing `schema_version`, `dist`,
#'   `params`, `support`, `data`, `custom`, `label`, and `metadata`. The `data`
#'   field is a [tibble::tibble()] with `amount`, `weight`, and `probability`
#'   columns, sorted by amount. Only observed amounts with positive empirical
#'   probability are retained.
#'
#' @seealso [severity_spec()], [empirical_frequency()],
#'   [distribution_moments()]
#'
#' @examples
#' # One claim amount per observation, including a zero payment.
#' losses <- c(0, 500, 500, 1250.50, 3000)
#' x <- empirical_severity(losses)
#' x
#' x$data
#' distribution_moments(x, order = 1:4)
#'
#' # Exactly the same empirical distribution from a summarized claim table.
#' claims <- tibble::tibble(
#'   payment = c(0, 500, 1250.50, 3000),
#'   claims = c(1, 2, 1, 1)
#' )
#' claims |>
#'   empirical_severity(amount = payment, weight = claims)
#'
#' # A tidyverse pipeline with individual claim observations.
#' tibble::tibble(claim_amount = c(150, 150, 375, 900)) |>
#'   empirical_severity(amount = claim_amount, label = "Reported severities")
#'
#' @export
empirical_severity <- function(data, amount = NULL, weight = NULL, label = NULL) {
  .severity_validate_label(label)

  amount_expr <- rlang::enquo(amount)
  weight_expr <- rlang::enquo(weight)

  if (is.data.frame(data)) {
    if (rlang::quo_is_null(amount_expr)) {
      stop("`amount` must identify a column when `data` is a data frame.",
           call. = FALSE)
    }
    amounts <- rlang::eval_tidy(amount_expr, data = data)
    if (length(amounts) != nrow(data)) {
      stop("`amount` must supply exactly one observation per data-frame row.",
           call. = FALSE)
    }
    weights <- if (rlang::quo_is_null(weight_expr)) {
      rep(1, length(amounts))
    } else {
      rlang::eval_tidy(weight_expr, data = data)
    }
  } else {
    if (!rlang::quo_is_null(amount_expr)) {
      stop("`amount` must be omitted when `data` is a numeric vector.",
           call. = FALSE)
    }
    amounts <- data
    weights <- if (rlang::quo_is_null(weight_expr)) {
      rep(1, length(amounts))
    } else {
      rlang::eval_tidy(weight_expr)
    }
  }

  if (!is.numeric(amounts) || !is.null(dim(amounts)) || !length(amounts) ||
      anyNA(amounts) || any(!is.finite(amounts)) || any(amounts < 0)) {
    stop("`amount` observations must be a nonempty numeric vector of finite nonnegative claim amounts.",
         call. = FALSE)
  }

  if (!is.numeric(weights) || !is.null(dim(weights)) ||
      length(weights) != length(amounts) || anyNA(weights) ||
      any(!is.finite(weights)) || any(weights < 0)) {
    stop("`weight` must contain one finite nonnegative value per claim amount.",
         call. = FALSE)
  }

  positive <- weights > 0
  if (!any(positive)) {
    stop("At least one observation must have positive weight.",
         call. = FALSE)
  }

  amounts <- as.numeric(amounts[positive])
  weights <- as.numeric(weights[positive])
  values <- sort(unique(amounts))

  # Group by numeric support index, without expanding claim multiplicities.
  group <- match(amounts, values)
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

  severity_table <- tibble::tibble(
    amount = values,
    weight = mass,
    probability = probs
  )

  structure(
    list(
      schema_version = 1L,
      dist = "empirical",
      params = list(values = values, probs = probs),
      support = list(
        type = "discrete",
        lower = min(values),
        upper = max(values)
      ),
      data = severity_table,
      custom = NULL,
      label = label,
      metadata = list()
    ),
    class = c("tidyact_severity_spec", "tidyact_distribution_spec")
  )
}


# Internal severity registry entry -----------------------------------------
# This entry is not another public constructor. Add
# `empirical = .empirical_severity_registry_entry(),` to .severity_registry().

.empirical_severity_registry_entry <- function() {
  list(
    validate = function(params) {
      stop("Use `empirical_severity()` to specify empirical claim amounts.",
           call. = FALSE)
    },
    support = function(params) {
      list(type = "discrete", lower = min(params$values),
           upper = max(params$values))
    },
    pdf = NULL,
    pmf = function(x, params, log = FALSE) {
      if (!is.numeric(x)) {
        stop("`x` must be numeric.", call. = FALSE)
      }
      index <- match(x, params$values)
      out <- numeric(length(x))
      found <- !is.na(index)
      out[found] <- params$probs[index[found]]
      out[is.na(x)] <- NA_real_
      if (log) log(out) else out
    },
    cdf = function(q, params, lower.tail = TRUE, log.p = FALSE) {
      if (!is.numeric(q)) {
        stop("`q` must be numeric.", call. = FALSE)
      }
      out <- rep(NA_real_, length(q))
      ok <- !is.na(q)
      if (any(ok)) {
        cumulative <- cumsum(params$probs)
        cumulative[length(cumulative)] <- 1
        indices <- findInterval(q[ok], params$values)
        result <- numeric(length(indices))
        inside <- indices > 0L
        result[inside] <- cumulative[indices[inside]]
        out[ok] <- result
      }
      if (!lower.tail) out <- 1 - out
      if (log.p) log(out) else out
    },
    quantile = function(prob, params, lower.tail = TRUE, log.p = FALSE) {
      pp <- .severity_valid_prob(prob, log.p)
      if (log.p) pp <- exp(pp)
      if (!lower.tail) pp <- 1 - pp
      cumulative <- cumsum(params$probs)
      cumulative[length(cumulative)] <- 1
      vapply(pp, function(target) {
        if (is.na(target)) return(NaN)
        if (target <= 0) return(params$values[1L])
        if (target >= 1) return(params$values[length(params$values)])
        index <- which(cumulative >= target)[1L]
        if (is.na(index)) params$values[length(params$values)] else
          params$values[index]
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
      .severity_check_order(order)
      if (order == 0) return(1)
      sum(params$probs * params$values^order)
    },
    limited_moment = function(order, limit, params) {
      .severity_check_order(order)
      .severity_check_limit(limit)
      if (order == 0) return(1)
      sum(params$probs * pmin(params$values, limit)^order)
    },
    mgf = function(t, params) {
      if (!is.numeric(t)) stop("`t` must be numeric.", call. = FALSE)
      vapply(t, function(z) sum(params$probs * exp(z * params$values)),
             numeric(1))
    },
    laplace = function(t, params) {
      .severity_check_laplace_argument(t)
      vapply(t, function(z) sum(params$probs * exp(-z * params$values)),
             numeric(1))
    }
  )
}
