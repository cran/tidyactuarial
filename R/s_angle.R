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
#' every \eqn{1/k} year. The function returns the accumulation factor only,
#' assuming a unit payment at each payment time.
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
#' @param type Character vector indicating the rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the compounding frequency
#'   for nominal rates. Ignored for \code{"effective"} and \code{"force"}.
#' @param deferral_years Numeric vector of deferral times in years.
#'   Must be greater than or equal to 0. Under the adopted horizon convention,
#'   this is metadata only for accumulation factors.
#' @param timing Character vector. One of \code{"immediate"},
#'   \code{"due"}, or \code{"continuous"}.
#'
#' @return Numeric vector of accumulation factors.
#'
#' @details
#' The function first converts the supplied rate to the equivalent annual
#' effective interest rate using \code{\link{standardize_interest}}.
#'
#' For finite discrete annuities:
#' \deqn{s_n = \frac{(1+i)^n - 1}{i}}{s_n = ((1+i)^n - 1) / i}
#'
#' For due annuities:
#' \deqn{\ddot{s}_n = (1+i) s_n}{s_n_due = (1+i) * s_n}
#'
#' For continuous annuities:
#' \deqn{\bar{s}_n = \frac{e^{\delta n} - 1}{\delta}}{s_n_bar = (exp(delta*n) - 1) / delta}
#'
#' Input vectors must have length 1 or a common length.
#' Missing values are propagated.
#'
#' @seealso \code{\link{s_angle_tbl}}, \code{\link{a_angle}},
#'   \code{\link{standardize_interest}}
#'
#' @family annuities
#'
#' @examples
#' # Simple scalar examples
#' s_angle(n_years = 10, rate = 0.05, type = "effective")
#' s_angle(n_years = 10, rate = 0.06, type = "nominal_interest", m = 12,
#'         payments_per_year = 12)
#' s_angle(n_years = 15, rate = 0.04, type = "force", timing = "continuous")
#'
#' # Medium vectorized example
#' s_angle(
#'   n_years = c(5, 10, 20),
#'   payments_per_year = c(1, 12, 1),
#'   rate = c(0.05, 0.06, 0.04),
#'   type = c("effective", "nominal_interest", "force"),
#'   m = c(1, 12, 1),
#'   deferral_years = c(0, 2, 3),
#'   timing = c("immediate", "due", "continuous")
#' )
#'
#' # Use inside a data pipeline
#' if (requireNamespace("dplyr", quietly = TRUE) &&
#'     requireNamespace("tibble", quietly = TRUE)) {
#'   contracts <- tibble::tibble(
#'     n_years = c(10, 15, 20),
#'     rate = c(0.05, 0.06, 0.04),
#'     type = c("effective", "nominal_interest", "force"),
#'     m = c(1, 12, 1),
#'     payments_per_year = c(1, 12, NA),
#'     deferral_years = c(0, 2, 3),
#'     timing = c("immediate", "due", "continuous")
#'   )
#'
#'   dplyr::mutate(
#'     contracts,
#'     factor = s_angle(
#'       n_years = n_years,
#'       payments_per_year = dplyr::coalesce(payments_per_year, 1L),
#'       rate = rate,
#'       type = type,
#'       m = m,
#'       deferral_years = deferral_years,
#'       timing = timing
#'     )
#'   )
#' }
#'
#' @export
s_angle <- function(
    n_years,
    payments_per_year = 1L,
    rate,
    type = "effective",
    m = 1L,
    deferral_years = 0,
    timing = "immediate"
) {
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
  if (!is.character(type)) {
    stop("`type` must be a character vector.", call. = FALSE)
  }
  if (!is.character(timing)) {
    stop("`timing` must be a character vector.", call. = FALSE)
  }

  # --- Determine common size ---
  size <- max(
    length(n_years),
    length(payments_per_year),
    length(rate),
    length(type),
    length(m),
    length(deferral_years),
    length(timing)
  )

  valid_size <- function(x) length(x) %in% c(1L, size)

  if (!valid_size(n_years) || !valid_size(payments_per_year) || !valid_size(rate) ||
      !valid_size(type) || !valid_size(m) || !valid_size(deferral_years) ||
      !valid_size(timing)) {
    stop(
      "`n_years`, `payments_per_year`, `rate`, `type`, `m`, ",
      "`deferral_years`, and `timing` must have length 1 or a common length.",
      call. = FALSE
    )
  }

  # --- Recycle ---
  n_years           <- rep_len(n_years, size)
  payments_per_year <- rep_len(payments_per_year, size)
  rate              <- rep_len(rate, size)
  type              <- rep_len(type, size)
  m                 <- rep_len(m, size)
  deferral_years    <- rep_len(deferral_years, size)
  timing            <- rep_len(timing, size)

  # --- Value-level validation ---
  bad_n <- !is.na(n_years) & (!is.finite(n_years) | n_years <= 0)
  if (any(bad_n)) {
    stop("`n_years` must contain only finite values greater than 0 or NA.", call. = FALSE)
  }

  bad_deferral <- !is.na(deferral_years) & (!is.finite(deferral_years) | deferral_years < 0)
  if (any(bad_deferral)) {
    stop("`deferral_years` must contain only finite values >= 0 or NA.", call. = FALSE)
  }

  timing <- tolower(timing)
  valid_timing <- c("immediate", "due", "continuous")
  bad_timing <- !is.na(timing) & !(timing %in% valid_timing)
  if (any(bad_timing)) {
    stop("`timing` must be 'immediate', 'due', or 'continuous'.", call. = FALSE)
  }

  # --- Convert rates ---
  i_effective <- standardize_interest(type = type, rate = rate, m = m)
  delta <- log1p(i_effective)

  out <- rep(NA_real_, size)
  eps <- 1e-12

  # --- Identify computable elements ---
  is_na_elem <- is.na(n_years) | is.na(rate) | is.na(i_effective) |
    is.na(deferral_years) | is.na(timing)
  disc_idx <- which(!is_na_elem)

  for (j in disc_idx) {
    if (timing[j] == "continuous") {
      out[j] <- if (abs(delta[j]) < eps) {
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
    n_periods <- round(n_periods_raw)
    if (abs(n_periods_raw - n_periods) > 1e-10) {
      stop(
        "For discrete annuities, `n_years * payments_per_year` must be an integer.",
        call. = FALSE
      )
    }

    deferral_periods_raw <- deferral_years[j] * k
    deferral_periods <- round(deferral_periods_raw)
    if (abs(deferral_periods_raw - deferral_periods) > 1e-10) {
      stop(
        "For discrete annuities, `deferral_years * payments_per_year` must be an integer.",
        call. = FALSE
      )
    }

    i_period <- (1 + i_effective[j])^(1 / k) - 1

    base <- if (abs(i_period) < eps) {
      n_periods
    } else {
      ((1 + i_period)^n_periods - 1) / i_period
    }

    if (timing[j] == "due") {
      base <- (1 + i_period) * base
    }

    out[j] <- base
  }

  out
}


#' Level annuity accumulation details in tibble form
#'
#' Computes the actuarial accumulation factor for a level annuity and returns
#' a tibble with the main input values, implied rates, accumulation factor,
#' payment amount, and future value.
#'
#' This is a reporting wrapper around \code{\link{s_angle}}. The accumulation
#' factor assumes unit payments. The future value is then computed as
#' \deqn{FV = R \times s}{FV = R * s}
#' where \eqn{R} is the payment amount and \eqn{s} is the accumulation factor.
#'
#' Under the adopted horizon convention, \code{deferral_years} is metadata only
#' and does not affect the accumulation factor.
#'
#' @param n_years Numeric vector of payment durations in years.
#' @param payments_per_year Positive integer vector giving the number of
#'   discrete payments per year. Ignored for continuous annuities.
#' @param rate Numeric vector of rate values.
#' @param type Character vector indicating the rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the compounding frequency
#'   for nominal rates.
#' @param deferral_years Numeric vector of deferral times in years.
#' @param timing Character vector. One of \code{"immediate"},
#'   \code{"due"}, or \code{"continuous"}.
#' @param payment Numeric vector of level payment amounts.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{n_years}{Payment duration in years.}
#'   \item{payments_per_year}{Number of payments per year.}
#'   \item{deferral_years}{Deferral period in years.}
#'   \item{timing}{Payment timing convention.}
#'   \item{rate_input}{Original supplied rate.}
#'   \item{rate_type}{Type of supplied rate.}
#'   \item{m}{Compounding frequency for nominal rates.}
#'   \item{i_effective}{Equivalent annual effective rate.}
#'   \item{delta}{Equivalent force of interest.}
#'   \item{i_period}{Equivalent per-payment effective rate for discrete annuities.}
#'   \item{v_period}{Equivalent per-payment discount factor for discrete annuities.}
#'   \item{n_periods}{Number of payment periods for discrete annuities.}
#'   \item{deferral_periods}{Number of deferred periods for discrete annuities.}
#'   \item{accumulation_factor}{Computed accumulation factor.}
#'   \item{payment}{Level payment amount.}
#'   \item{future_value}{Future value of the annuity.}
#' }
#'
#' @seealso \code{\link{s_angle}}, \code{\link{a_angle_tbl}},
#'   \code{\link{standardize_interest}}, \code{\link{future_value}}
#'
#' @family annuities
#'
#' @examples
#' # Simple scalar example
#' s_angle_tbl(n_years = 10, rate = 0.05, payment = 1000)
#'
#' # Medium vectorized example
#' s_angle_tbl(
#'   n_years = c(10, 15, 20),
#'   payments_per_year = c(1, 12, 1),
#'   rate = c(0.05, 0.06, 0.04),
#'   type = c("effective", "nominal_interest", "force"),
#'   m = c(1, 12, 1),
#'   deferral_years = c(0, 2, 3),
#'   timing = c("immediate", "due", "continuous"),
#'   payment = c(1000, 200, 5000)
#' )
#'
#' @export
s_angle_tbl <- function(
    n_years,
    payments_per_year = 1L,
    rate,
    type = "effective",
    m = 1L,
    deferral_years = 0,
    timing = "immediate",
    payment = 1
) {
  if (missing(n_years)) {
    stop("`n_years` must be provided.", call. = FALSE)
  }
  if (missing(rate)) {
    stop("`rate` must be provided.", call. = FALSE)
  }
  if (!is.numeric(payment)) {
    stop("`payment` must be a numeric vector.", call. = FALSE)
  }

  size <- max(
    length(n_years),
    length(payments_per_year),
    length(rate),
    length(type),
    length(m),
    length(deferral_years),
    length(timing),
    length(payment)
  )

  valid_size <- function(x) length(x) %in% c(1L, size)

  if (!valid_size(n_years) || !valid_size(payments_per_year) || !valid_size(rate) ||
      !valid_size(type) || !valid_size(m) || !valid_size(deferral_years) ||
      !valid_size(timing) || !valid_size(payment)) {
    stop(
      "`n_years`, `payments_per_year`, `rate`, `type`, `m`, ",
      "`deferral_years`, `timing`, and `payment` must have length 1 ",
      "or a common length.",
      call. = FALSE
    )
  }

  n_years           <- rep_len(n_years, size)
  payments_per_year <- rep_len(payments_per_year, size)
  rate              <- rep_len(rate, size)
  type              <- rep_len(type, size)
  m                 <- rep_len(m, size)
  deferral_years    <- rep_len(deferral_years, size)
  timing            <- rep_len(timing, size)
  payment           <- rep_len(payment, size)

  bad_payment <- !is.na(payment) & !is.finite(payment)
  if (any(bad_payment)) {
    stop("`payment` must contain only finite numeric values or NA.", call. = FALSE)
  }

  accumulation_factor <- s_angle(
    n_years = n_years,
    payments_per_year = payments_per_year,
    rate = rate,
    type = type,
    m = m,
    deferral_years = deferral_years,
    timing = timing
  )

  i_effective <- standardize_interest(type = type, rate = rate, m = m)
  delta <- log1p(i_effective)

  timing_lower <- tolower(timing)
  discrete <- timing_lower %in% c("immediate", "due")

  i_period <- rep(NA_real_, size)
  v_period <- rep(NA_real_, size)
  n_periods <- rep(NA_real_, size)
  deferral_periods <- rep(NA_real_, size)

  if (any(discrete, na.rm = TRUE)) {
    k <- payments_per_year[discrete]
    i_period[discrete] <- (1 + i_effective[discrete])^(1 / k) - 1
    v_period[discrete] <- 1 / (1 + i_period[discrete])
    n_periods[discrete] <- n_years[discrete] * k
    deferral_periods[discrete] <- deferral_years[discrete] * k
  }

  payments_per_year_out <- as.integer(payments_per_year)
  payments_per_year_out[timing_lower == "continuous"] <- NA_integer_

  tibble::tibble(
    n_years = n_years,
    payments_per_year = payments_per_year_out,
    deferral_years = deferral_years,
    timing = timing_lower,
    rate_input = rate,
    rate_type = type,
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
