#' Arithmetic annuity factor
#'
#' Computes the actuarial present value or accumulated value factor for an
#' arithmetic annuity.
#'
#' This function covers increasing, decreasing, and custom arithmetic payment
#' patterns. It replaces the more specific increasing and decreasing annuity
#' functions.
#'
#' Supported payment patterns:
#' \itemize{
#'   \item \code{"increasing"}: payments \eqn{1, 2, \ldots, n}.
#'   \item \code{"decreasing"}: payments \eqn{n, n-1, \ldots, 1}.
#'   \item \code{"custom"}: payments following
#'     \eqn{P_r = P_1 + (r - 1)g}.
#' }
#'
#' Supported timing conventions:
#' \itemize{
#'   \item \code{"immediate"}: payments at the end of each period.
#'   \item \code{"due"}: payments at the beginning of each period.
#' }
#'
#' @param n_years Numeric vector of payment durations in years.
#'   Each value must be positive and finite.
#' @param payments_per_year Positive integer vector giving the number of
#'   payments per year.
#' @param rate Numeric vector of rate values.
#' @param rate_type Character vector indicating the rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the compounding frequency
#'   for nominal rates.
#' @param deferral_years Numeric vector of deferral times in years.
#'   Must be greater than or equal to 0.
#' @param timing Character vector. One of \code{"immediate"} or \code{"due"}.
#' @param pattern Character vector. One of \code{"increasing"},
#'   \code{"decreasing"}, or \code{"custom"}.
#' @param first_payment Numeric vector. First payment for
#'   \code{pattern = "custom"}. Ignored for \code{"increasing"} and
#'   \code{"decreasing"}.
#' @param increment Numeric vector. Arithmetic increment for
#'   \code{pattern = "custom"}. Ignored for \code{"increasing"} and
#'   \code{"decreasing"}.
#' @param valuation Character string. Use \code{"present"} for present value
#'   or \code{"accumulated"} for accumulated value at the end of the term.
#' @param output Character string. Use \code{"value"} to return a numeric
#'   factor, or \code{"table"} to return a tibble with intermediate quantities.
#'
#' @return
#' If \code{output = "value"}, a numeric vector of arithmetic annuity factors.
#'
#' If \code{output = "table"}, a tibble with input values, equivalent rates,
#' period quantities, and both present and accumulated value factors.
#'
#' @details
#' The function first converts the supplied rate to the equivalent annual
#' effective interest rate using \code{\link{standardize_interest}}.
#'
#' For each scenario, the total number of payments is
#' \deqn{N = n \times k}
#' where \eqn{n} is \code{n_years} and \eqn{k} is \code{payments_per_year}.
#'
#' The present value is computed by summing the discounted payment stream.
#' The accumulated value is computed at the end of the annuity term. Under
#' this convention, \code{deferral_years} affects the present value factor but
#' not the accumulated value factor.
#'
#' @seealso \code{\link{a_angle}}, \code{\link{s_angle}},
#'   \code{\link{standardize_interest}}
#'
#' @family annuities
#'
#' @examples
#' # Increasing arithmetic annuity
#' annuity_arith(n_years = 10, rate = 0.05, pattern = "increasing")
#'
#' # Decreasing arithmetic annuity
#' annuity_arith(n_years = 10, rate = 0.05, pattern = "decreasing")
#'
#' # Custom arithmetic annuity
#' annuity_arith(
#'   n_years = 10,
#'   rate = 0.05,
#'   pattern = "custom",
#'   first_payment = 100,
#'   increment = 25
#' )
#'
#' # Accumulated value factor
#' annuity_arith(
#'   n_years = 10,
#'   rate = 0.05,
#'   pattern = "increasing",
#'   valuation = "accumulated"
#' )
#'
#' # Tibble output
#' annuity_arith(
#'   n_years = 10,
#'   rate = 0.05,
#'   pattern = "decreasing",
#'   output = "table"
#' )
#'
#' @export
annuity_arith <- function(
    n_years,
    payments_per_year = 1L,
    rate,
    rate_type = "effective",
    m = 1L,
    deferral_years = 0,
    timing = "immediate",
    pattern = "increasing",
    first_payment = 1,
    increment = 1,
    valuation = c("present", "accumulated"),
    output = c("value", "table")
) {
  valuation <- match.arg(valuation)
  output <- match.arg(output)

  if (missing(n_years)) {
    stop("`n_years` must be provided.", call. = FALSE)
  }
  if (missing(rate)) {
    stop("`rate` must be provided.", call. = FALSE)
  }

  # --- Early type validation ---
  if (!is.numeric(n_years)) {
    stop("`n_years` must be a numeric vector.", call. = FALSE)
  }
  if (!is.numeric(payments_per_year)) {
    stop("`payments_per_year` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(rate)) {
    stop("`rate` must be a numeric vector.", call. = FALSE)
  }
  if (!is.character(rate_type)) {
    stop("`rate_type` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(m)) {
    stop("`m` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(deferral_years)) {
    stop("`deferral_years` must be numeric.", call. = FALSE)
  }
  if (!is.character(timing)) {
    stop("`timing` must be a character vector.", call. = FALSE)
  }
  if (!is.character(pattern)) {
    stop("`pattern` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(first_payment)) {
    stop("`first_payment` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(increment)) {
    stop("`increment` must be numeric.", call. = FALSE)
  }

  # --- Determine common size ---
  size <- max(
    length(n_years),
    length(payments_per_year),
    length(rate),
    length(rate_type),
    length(m),
    length(deferral_years),
    length(timing),
    length(pattern),
    length(first_payment),
    length(increment),
    1L
  )

  valid_size <- function(x) length(x) %in% c(1L, size)

  if (!valid_size(n_years) ||
      !valid_size(payments_per_year) ||
      !valid_size(rate) ||
      !valid_size(rate_type) ||
      !valid_size(m) ||
      !valid_size(deferral_years) ||
      !valid_size(timing) ||
      !valid_size(pattern) ||
      !valid_size(first_payment) ||
      !valid_size(increment)) {
    stop(
      "`n_years`, `payments_per_year`, `rate`, `rate_type`, `m`, ",
      "`deferral_years`, `timing`, `pattern`, `first_payment`, and ",
      "`increment` must have length 1 or a common length.",
      call. = FALSE
    )
  }

  # --- Recycle ---
  n_years <- rep_len(n_years, size)
  payments_per_year <- rep_len(payments_per_year, size)
  rate <- rep_len(rate, size)
  rate_type <- rep_len(rate_type, size)
  m <- rep_len(m, size)
  deferral_years <- rep_len(deferral_years, size)
  timing <- rep_len(timing, size)
  pattern <- rep_len(pattern, size)
  first_payment <- rep_len(first_payment, size)
  increment <- rep_len(increment, size)

  # --- Value-level validation ---
  bad_n <- !is.na(n_years) & (!is.finite(n_years) | n_years <= 0)
  if (any(bad_n)) {
    stop("`n_years` must contain only finite values greater than 0 or NA.", call. = FALSE)
  }

  bad_deferral <- !is.na(deferral_years) &
    (!is.finite(deferral_years) | deferral_years < 0)

  if (any(bad_deferral)) {
    stop("`deferral_years` must contain only finite values >= 0 or NA.", call. = FALSE)
  }

  bad_first <- !is.na(first_payment) & !is.finite(first_payment)
  if (any(bad_first)) {
    stop("`first_payment` must contain only finite numeric values or NA.", call. = FALSE)
  }

  bad_increment <- !is.na(increment) & !is.finite(increment)
  if (any(bad_increment)) {
    stop("`increment` must contain only finite numeric values or NA.", call. = FALSE)
  }

  timing <- tolower(timing)
  valid_timing <- c("immediate", "due")

  bad_timing <- !is.na(timing) & !(timing %in% valid_timing)
  if (any(bad_timing)) {
    stop("`timing` must be 'immediate' or 'due'.", call. = FALSE)
  }

  pattern <- tolower(pattern)
  valid_pattern <- c("increasing", "decreasing", "custom")

  bad_pattern <- !is.na(pattern) & !(pattern %in% valid_pattern)
  if (any(bad_pattern)) {
    stop("`pattern` must be 'increasing', 'decreasing', or 'custom'.", call. = FALSE)
  }

  # --- Convert to effective annual rate ---
  i_effective <- standardize_interest(type = rate_type, rate = rate, m = m)

  present_value_factor <- rep(NA_real_, size)
  accumulated_value_factor <- rep(NA_real_, size)
  selected_factor <- rep(NA_real_, size)

  i_period <- rep(NA_real_, size)
  v_period <- rep(NA_real_, size)
  n_periods <- rep(NA_real_, size)
  deferral_periods <- rep(NA_real_, size)
  first_payment_used <- rep(NA_real_, size)
  increment_used <- rep(NA_real_, size)

  is_na_elem <- is.na(n_years) |
    is.na(payments_per_year) |
    is.na(rate) |
    is.na(i_effective) |
    is.na(deferral_years) |
    is.na(timing) |
    is.na(pattern)

  computable_idx <- which(!is_na_elem)

  for (j in computable_idx) {
    k <- payments_per_year[j]

    if (is.na(k) || !is.finite(k) || k < 1 || k != floor(k)) {
      stop("`payments_per_year` must be a positive integer.", call. = FALSE)
    }

    k <- as.integer(k)

    n_raw <- n_years[j] * k
    n <- round(n_raw)

    if (abs(n_raw - n) > 1e-10) {
      stop(
        "For discrete annuities, `n_years * payments_per_year` must be an integer.",
        call. = FALSE
      )
    }

    h_raw <- deferral_years[j] * k
    h <- round(h_raw)

    if (abs(h_raw - h) > 1e-10) {
      stop(
        "For discrete annuities, `deferral_years * payments_per_year` must be an integer.",
        call. = FALSE
      )
    }

    i_period[j] <- (1 + i_effective[j])^(1 / k) - 1
    v_period[j] <- 1 / (1 + i_period[j])
    n_periods[j] <- n
    deferral_periods[j] <- h

    if (pattern[j] == "increasing") {
      p1 <- 1
      g <- 1
    } else if (pattern[j] == "decreasing") {
      p1 <- n
      g <- -1
    } else {
      if (is.na(first_payment[j]) || is.na(increment[j])) {
        stop(
          "`first_payment` and `increment` must not be NA when `pattern = 'custom'`.",
          call. = FALSE
        )
      }

      p1 <- first_payment[j]
      g <- increment[j]
    }

    first_payment_used[j] <- p1
    increment_used[j] <- g

    payment_index <- seq_len(n)
    payments <- p1 + (payment_index - 1) * g

    local_times <- if (timing[j] == "immediate") {
      payment_index
    } else {
      payment_index - 1
    }

    present_value_factor[j] <- sum(payments * v_period[j]^(h + local_times))

    accumulated_value_factor[j] <- sum(
      payments * (1 + i_period[j])^(n - local_times)
    )

    selected_factor[j] <- if (valuation == "present") {
      present_value_factor[j]
    } else {
      accumulated_value_factor[j]
    }
  }

  if (output == "value") {
    return(selected_factor)
  }

  tibble::tibble(
    n_years = n_years,
    payments_per_year = as.integer(payments_per_year),
    deferral_years = deferral_years,
    timing = timing,
    pattern = pattern,
    valuation = valuation,
    rate_input = rate,
    rate_type = rate_type,
    m = m,
    i_effective = i_effective,
    i_period = i_period,
    v_period = v_period,
    n_periods = n_periods,
    deferral_periods = deferral_periods,
    first_payment = first_payment_used,
    increment = increment_used,
    present_value_factor = present_value_factor,
    accumulated_value_factor = accumulated_value_factor,
    annuity_factor = selected_factor
  )
}
