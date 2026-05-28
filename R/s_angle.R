#' Level annuity accumulation factor s-angle-n
#'
#' Computes the actuarial accumulation factor for a level annuity.
#'
#' Supported timing conventions:
#' \itemize{
#'   \item \code{"immediate"}: annuity-immediate with discrete payments.
#'   \item \code{"due"}: annuity-due with discrete payments.
#'   \item \code{"continuous"}: continuous annuity.
#' }
#'
#' For discrete annuities, \code{payments_per_year = k} means payments are made
#' every \eqn{1/k} year. The function returns the accumulation factor, assuming
#' a unit payment at each payment time.
#'
#' Horizon convention:
#' the future value is measured at the time of the last payment. Under this
#' convention, a pure deferral that shifts the entire payment block forward in
#' time does not change the accumulation factor when the payment pattern is
#' otherwise unchanged. Therefore, \code{deferral_years} is recorded and
#' validated, but it does not modify the factor.
#'
#' The future value of a perpetuity diverges, so perpetuities are not supported
#' in \code{s_angle()}.
#'
#' @param n_years Numeric vector of payment durations in years.
#'   Each value must be positive and finite.
#' @param payments_per_year Positive integer vector giving the number of
#'   discrete payments per year. Ignored for continuous annuities.
#' @param rate Numeric vector of rate values.
#' @param rate_type Character vector indicating the rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the compounding frequency
#'   for nominal rates. Ignored for \code{"effective"} and \code{"force"}.
#' @param deferral_years Numeric vector of deferral times in years.
#'   Must be greater than or equal to 0. Under the adopted horizon convention,
#'   this is metadata only for accumulation factors.
#' @param timing Character vector. One of \code{"immediate"},
#'   \code{"due"}, or \code{"continuous"}.
#' @param payment Numeric vector of level payment amounts. Used only when
#'   \code{output = "table"} to report the corresponding future value.
#'   The accumulation factor itself is always computed for unit payments.
#' @param output Character string. Use \code{"value"} to return a numeric
#'   accumulation factor, or \code{"table"} to return a tibble with intermediate
#'   calculations.
#'
#' @return
#' If \code{output = "value"}, a numeric vector of accumulation factors.
#'
#' If \code{output = "table"}, a tibble with input values, equivalent rates,
#' accumulation factors, payment amounts, and future values.
#'
#' @details
#' The function first converts the supplied rate to the equivalent annual
#' effective interest rate using \code{\link{standardize_interest}}.
#'
#' For finite discrete annuities:
#' \deqn{s_{\overline{n|}} = \frac{(1+i)^n - 1}{i}}
#'
#' For due annuities:
#' \deqn{\ddot{s}_{\overline{n|}} = (1+i)s_{\overline{n|}}}
#'
#' For continuous annuities:
#' \deqn{\bar{s}_{\overline{n|}} = \frac{e^{\delta n} - 1}{\delta}}
#'
#' Input vectors must have length 1 or a common length. Missing values are
#' propagated.
#'
#' @seealso \code{\link{a_angle}}, \code{\link{standardize_interest}},
#'   \code{\link{future_value}}
#'
#' @family annuities
#'
#' @examples
#' # Numeric accumulation factor
#' s_angle(n_years = 10, rate = 0.05)
#'
#' # Nominal interest converted monthly, with monthly payments
#' s_angle(
#'   n_years = 10,
#'   rate = 0.06,
#'   rate_type = "nominal_interest",
#'   m = 12,
#'   payments_per_year = 12
#' )
#'
#' # Continuous annuity
#' s_angle(
#'   n_years = 15,
#'   rate = 0.04,
#'   rate_type = "force",
#'   timing = "continuous"
#' )
#'
#' # Tibble output for teaching or auditing
#' s_angle(
#'   n_years = 10,
#'   rate = 0.05,
#'   payment = 1000,
#'   output = "table"
#' )
#'
#' # Vectorized example
#' s_angle(
#'   n_years = c(5, 10, 20),
#'   payments_per_year = c(1, 12, 1),
#'   rate = c(0.05, 0.06, 0.04),
#'   rate_type = c("effective", "nominal_interest", "force"),
#'   m = c(1, 12, 1),
#'   deferral_years = c(0, 2, 3),
#'   timing = c("immediate", "due", "continuous")
#' )
#'
#' @export
s_angle <- function(
    n_years,
    payments_per_year = 1L,
    rate,
    rate_type = "effective",
    m = 1L,
    deferral_years = 0,
    timing = "immediate",
    payment = 1,
    output = c("value", "table")
) {
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
  if (!is.numeric(m)) {
    stop("`m` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(deferral_years)) {
    stop("`deferral_years` must be numeric.", call. = FALSE)
  }
  if (!is.character(rate_type)) {
    stop("`rate_type` must be a character vector.", call. = FALSE)
  }
  if (!is.character(timing)) {
    stop("`timing` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(payment)) {
    stop("`payment` must be a numeric vector.", call. = FALSE)
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
    length(payment),
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
      !valid_size(payment)) {
    stop(
      "`n_years`, `payments_per_year`, `rate`, `rate_type`, `m`, ",
      "`deferral_years`, `timing`, and `payment` must have length 1 ",
      "or a common length.",
      call. = FALSE
    )
  }

  # --- Recycle ---
  n_years           <- rep_len(n_years, size)
  payments_per_year <- rep_len(payments_per_year, size)
  rate              <- rep_len(rate, size)
  rate_type         <- rep_len(rate_type, size)
  m                 <- rep_len(m, size)
  deferral_years    <- rep_len(deferral_years, size)
  timing            <- rep_len(timing, size)
  payment           <- rep_len(payment, size)

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

  bad_payment <- !is.na(payment) & !is.finite(payment)
  if (any(bad_payment)) {
    stop("`payment` must contain only finite numeric values or NA.", call. = FALSE)
  }

  timing <- tolower(timing)
  valid_timing <- c("immediate", "due", "continuous")

  bad_timing <- !is.na(timing) & !(timing %in% valid_timing)
  if (any(bad_timing)) {
    stop("`timing` must be 'immediate', 'due', or 'continuous'.", call. = FALSE)
  }

  # --- Convert rates ---
  i_effective <- standardize_interest(type = rate_type, rate = rate, m = m)
  delta <- log1p(i_effective)

  accumulation_factor <- rep(NA_real_, size)
  eps <- 1e-12

  # --- Intermediate output columns ---
  i_period <- rep(NA_real_, size)
  v_period <- rep(NA_real_, size)
  n_periods <- rep(NA_real_, size)
  deferral_periods <- rep(NA_real_, size)

  # --- Identify computable elements ---
  is_na_elem <- is.na(n_years) |
    is.na(rate) |
    is.na(i_effective) |
    is.na(deferral_years) |
    is.na(timing)

  computable_idx <- which(!is_na_elem)

  for (j in computable_idx) {
    if (timing[j] == "continuous") {
      accumulation_factor[j] <- if (abs(delta[j]) < eps) {
        n_years[j]
      } else {
        expm1(delta[j] * n_years[j]) / delta[j]
      }

      next
    }

    k <- payments_per_year[j]

    if (is.na(k) || !is.finite(k) || k < 1 || k != floor(k)) {
      stop(
        "`payments_per_year` must be a positive integer for discrete annuities.",
        call. = FALSE
      )
    }

    k <- as.integer(k)

    n_periods_raw <- n_years[j] * k
    n_periods[j] <- round(n_periods_raw)

    if (abs(n_periods_raw - n_periods[j]) > 1e-10) {
      stop(
        "For discrete annuities, `n_years * payments_per_year` must be an integer.",
        call. = FALSE
      )
    }

    deferral_periods_raw <- deferral_years[j] * k
    deferral_periods[j] <- round(deferral_periods_raw)

    if (abs(deferral_periods_raw - deferral_periods[j]) > 1e-10) {
      stop(
        "For discrete annuities, `deferral_years * payments_per_year` must be an integer.",
        call. = FALSE
      )
    }

    i_period[j] <- (1 + i_effective[j])^(1 / k) - 1
    v_period[j] <- 1 / (1 + i_period[j])

    base <- if (abs(i_period[j]) < eps) {
      n_periods[j]
    } else {
      ((1 + i_period[j])^n_periods[j] - 1) / i_period[j]
    }

    if (timing[j] == "due") {
      base <- (1 + i_period[j]) * base
    }

    accumulation_factor[j] <- base
  }

  if (output == "value") {
    return(accumulation_factor)
  }

  payments_per_year_out <- as.integer(payments_per_year)
  payments_per_year_out[timing == "continuous"] <- NA_integer_

  tibble::tibble(
    n_years = n_years,
    payments_per_year = payments_per_year_out,
    deferral_years = deferral_years,
    timing = timing,
    rate_input = rate,
    rate_type = rate_type,
    m = m,
    i_effective = i_effective,
    delta = delta,
    i_period = i_period,
    v_period = v_period,
    n_periods = n_periods,
    deferral_periods = deferral_periods,
    accumulation_factor = accumulation_factor,
    payment = payment,
    future_value = payment * accumulation_factor
  )
}
