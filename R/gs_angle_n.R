#' Geometric annuity accumulation factor gs-angle-n
#'
#' Computes the actuarial accumulation factor for a geometric annuity.
#'
#' The payment pattern is geometric. For a unit first payment, the payments are:
#' \itemize{
#'   \item immediate: \eqn{1, (1+g_p), (1+g_p)^2, \ldots}{1, (1+gp), (1+gp)^2, ...}
#'   \item due: same geometric pattern, but shifted one period earlier
#' }
#' where \eqn{g_p}{gp} is the effective growth rate per payment period.
#'
#' Supported timing conventions:
#' \itemize{
#'   \item \code{"immediate"}: payments at the end of each period.
#'   \item \code{"due"}: payments at the beginning of each period.
#' }
#'
#' Horizon convention:
#' the future value is measured at the standard terminal horizon for the annuity
#' accumulation factor. Under this convention, a pure deferral that shifts the
#' entire payment block forward in time does not change the accumulation factor
#' when the payment pattern is otherwise unchanged. Therefore, \code{deferral_years}
#' is recorded and validated, but it does not modify the factor.
#'
#' The interest rate and the growth rate may each be supplied in FM-style notation:
#' \itemize{
#'   \item annual effective rate,
#'   \item nominal annual interest rate,
#'   \item nominal annual discount rate,
#'   \item force of interest.
#' }
#'
#' Internally, both are first converted to annual effective rates and then to
#' effective rates per payment period.
#'
#' @param n_years Numeric vector of payment durations in years.
#'   Each value must be positive and finite.
#' @param payments_per_year Positive integer vector giving the number of
#'   payments per year.
#' @param rate Numeric vector of annual rate values for accumulation.
#' @param type Character vector indicating the annual accumulation-rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the compounding frequency for
#'   nominal accumulation-rate inputs.
#' @param growth Numeric vector of annual growth-rate values.
#' @param growth_type Character vector indicating the annual growth-rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param growth_m Positive integer vector giving the compounding frequency for
#'   nominal growth-rate inputs.
#' @param deferral_years Numeric vector of deferral times in years.
#'   Must be greater than or equal to 0. Under the adopted horizon convention,
#'   this is metadata only for accumulation factors.
#' @param timing Character vector. One of \code{"immediate"} or \code{"due"}.
#'
#' @return Numeric vector of geometric accumulation factors.
#'
#' @details
#' Let \eqn{i_p}{ip} be the effective accumulation rate per payment period,
#' and \eqn{g_p}{gp} the effective growth rate per payment period.
#'
#' For a finite geometric annuity-immediate with \eqn{n} payment periods:
#' \deqn{gs_n = \sum_{t=1}^{n} (1+g_p)^{t-1}(1+i_p)^{n-t}}{gs_n = sum((1+gp)^(t-1) * (1+ip)^(n-t), t=1..n)}
#'
#' If \eqn{i_p \neq g_p}{ip != gp}, then
#' \deqn{gs_n = \frac{(1+i_p)^n - (1+g_p)^n}{i_p-g_p}}{gs_n = ((1+ip)^n - (1+gp)^n) / (ip - gp)}
#'
#' If \eqn{i_p = g_p}{ip = gp}, the limiting formula is
#' \deqn{gs_n = n(1+i_p)^{n-1}}{gs_n = n * (1+ip)^(n-1)}
#'
#' The annuity-due version is obtained by multiplying the immediate factor by
#' \eqn{1+i_p}{1+ip}.
#'
#' Input vectors must have length 1 or a common length.
#' Missing values are propagated.
#'
#' Geometric perpetuities are not supported by this accumulation-factor
#' function.
#'
#' @seealso \code{\link{gs_angle_tbl}}, \code{\link{s_angle}},
#'   \code{\link{ga_angle}}, \code{\link{standardize_interest}}
#'
#' @family annuities
#'
#' @examples
#' # Simple scalar example
#' gs_angle(
#'   n_years = 10,
#'   rate = 0.05,
#'   type = "effective",
#'   growth = 0.02,
#'   growth_type = "effective",
#'   timing = "immediate"
#' )
#'
#' # Growth = 0 collapses to a level accumulation factor
#' gs_angle(
#'   n_years = 10,
#'   rate = 0.05,
#'   type = "effective",
#'   growth = 0,
#'   growth_type = "effective",
#'   timing = "immediate"
#' )
#'
#' # Medium vectorized example
#' gs_angle(
#'   n_years = c(10, 20),
#'   payments_per_year = c(1, 12),
#'   rate = c(0.05, 0.08),
#'   type = c("effective", "nominal_interest"),
#'   m = c(1, 12),
#'   growth = c(0.02, 0.03),
#'   growth_type = c("effective", "effective"),
#'   growth_m = c(1, 1),
#'   deferral_years = c(0, 2),
#'   timing = c("immediate", "due")
#' )
#'
#' @export
gs_angle <- function(
    n_years,
    payments_per_year = 1L,
    rate,
    type = "effective",
    m = 1L,
    growth = 0,
    growth_type = "effective",
    growth_m = 1L,
    deferral_years = 0,
    timing = c("immediate", "due")
) {
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
  if (!is.character(type)) {
    stop("`type` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(m)) {
    stop("`m` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(growth)) {
    stop("`growth` must be a numeric vector.", call. = FALSE)
  }
  if (!is.character(growth_type)) {
    stop("`growth_type` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(growth_m)) {
    stop("`growth_m` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(deferral_years)) {
    stop("`deferral_years` must be numeric.", call. = FALSE)
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
    length(growth),
    length(growth_type),
    length(growth_m),
    length(deferral_years),
    length(timing)
  )

  valid_size <- function(x) length(x) %in% c(1L, size)

  if (!valid_size(n_years) || !valid_size(payments_per_year) || !valid_size(rate) ||
      !valid_size(type) || !valid_size(m) || !valid_size(growth) ||
      !valid_size(growth_type) || !valid_size(growth_m) ||
      !valid_size(deferral_years) || !valid_size(timing)) {
    stop(
      "`n_years`, `payments_per_year`, `rate`, `type`, `m`, `growth`, ",
      "`growth_type`, `growth_m`, `deferral_years`, and `timing` must have ",
      "length 1 or a common length.",
      call. = FALSE
    )
  }

  # --- Recycle ---
  n_years           <- rep_len(n_years, size)
  payments_per_year <- rep_len(payments_per_year, size)
  rate              <- rep_len(rate, size)
  type              <- rep_len(type, size)
  m                 <- rep_len(m, size)
  growth            <- rep_len(growth, size)
  growth_type       <- rep_len(growth_type, size)
  growth_m          <- rep_len(growth_m, size)
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
  valid_timing <- c("immediate", "due")
  bad_timing <- !is.na(timing) & !(timing %in% valid_timing)
  if (any(bad_timing)) {
    stop("`timing` must be 'immediate' or 'due'.", call. = FALSE)
  }

  # --- Convert rates ---
  i_effective <- standardize_interest(type = type, rate = rate, m = m)
  g_effective <- standardize_interest(type = growth_type, rate = growth, m = growth_m)

  out <- rep(NA_real_, size)
  eps <- 1e-12

  # --- Identify computable elements ---
  is_na_elem <- is.na(n_years) | is.na(rate) | is.na(growth) |
    is.na(i_effective) | is.na(g_effective) |
    is.na(deferral_years) | is.na(timing)
  disc_idx <- which(!is_na_elem)

  for (j in disc_idx) {
    k <- payments_per_year[j]
    if (is.na(k) || !is.finite(k) || k < 1 || k != floor(k)) {
      stop("`payments_per_year` must be a positive integer.", call. = FALSE)
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
    g_period <- (1 + g_effective[j])^(1 / k) - 1

    base <- if (abs(i_period - g_period) < eps) {
      n_periods * (1 + i_period)^(n_periods - 1)
    } else {
      ((1 + i_period)^n_periods - (1 + g_period)^n_periods) / (i_period - g_period)
    }

    if (timing[j] == "due") {
      base <- (1 + i_period) * base
    }

    out[j] <- base
  }

  out
}


#' Geometric annuity accumulation details in tibble form
#'
#' Computes the actuarial accumulation factor for a geometric annuity and
#' returns a tibble with the main input values, implied rates, factor,
#' payment amount, and future value.
#'
#' This is a reporting wrapper around \code{\link{gs_angle}}. The factor assumes
#' that the first payment equals 1. The actual future value is then computed as
#' \deqn{FV = P_1 \times gs}{FV = P1 * gs}
#' where \eqn{P_1}{P1} is the first payment of the geometric sequence.
#'
#' Under the adopted horizon convention, \code{deferral_years} is metadata only
#' and does not affect the accumulation factor.
#'
#' @param n_years Numeric vector of payment durations in years.
#' @param payments_per_year Positive integer vector giving the number of
#'   payments per year.
#' @param rate Numeric vector of annual rate values for accumulation.
#' @param type Character vector indicating the annual accumulation-rate type.
#' @param m Positive integer vector giving the compounding frequency for
#'   nominal accumulation-rate inputs.
#' @param growth Numeric vector of annual growth-rate values.
#' @param growth_type Character vector indicating the annual growth-rate type.
#' @param growth_m Positive integer vector giving the compounding frequency for
#'   nominal growth-rate inputs.
#' @param deferral_years Numeric vector of deferral times in years.
#' @param timing Character vector. One of \code{"immediate"} or \code{"due"}.
#' @param payment Numeric vector giving the first payment of the geometric sequence.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{n_years}{Payment duration in years.}
#'   \item{payments_per_year}{Number of payments per year.}
#'   \item{deferral_years}{Deferral period in years.}
#'   \item{timing}{Payment timing convention.}
#'   \item{rate_input}{Original supplied accumulation rate.}
#'   \item{rate_type}{Type of supplied accumulation rate.}
#'   \item{m}{Compounding frequency for nominal accumulation-rate inputs.}
#'   \item{growth_input}{Original supplied growth rate.}
#'   \item{growth_type}{Type of supplied growth rate.}
#'   \item{growth_m}{Compounding frequency for nominal growth-rate inputs.}
#'   \item{i_effective}{Equivalent annual effective accumulation rate.}
#'   \item{g_effective}{Equivalent annual effective growth rate.}
#'   \item{i_period}{Equivalent per-payment accumulation rate.}
#'   \item{g_period}{Equivalent per-payment growth rate.}
#'   \item{n_periods}{Number of payment periods.}
#'   \item{deferral_periods}{Number of deferred periods.}
#'   \item{gs_factor}{Computed geometric accumulation factor.}
#'   \item{payment}{First payment of the geometric sequence.}
#'   \item{future_value}{Future value of the geometric annuity.}
#' }
#'
#' @seealso \code{\link{gs_angle}}, \code{\link{s_angle}},
#'   \code{\link{ga_angle_tbl}}, \code{\link{standardize_interest}}
#'
#' @family annuities
#'
#' @examples
#' # Simple scalar example
#' gs_angle_tbl(
#'   n_years = 10,
#'   rate = 0.05,
#'   type = "effective",
#'   growth = 0.02,
#'   growth_type = "effective",
#'   payment = 100
#' )
#'
#' # Medium vectorized example
#' gs_angle_tbl(
#'   n_years = c(10, 20),
#'   payments_per_year = c(1, 12),
#'   rate = c(0.05, 0.08),
#'   type = c("effective", "nominal_interest"),
#'   m = c(1, 12),
#'   growth = c(0.02, 0.03),
#'   growth_type = c("effective", "effective"),
#'   growth_m = c(1, 1),
#'   deferral_years = c(0, 2),
#'   timing = c("immediate", "due"),
#'   payment = c(100, 50)
#' )
#'
#' @export
gs_angle_tbl <- function(
    n_years,
    payments_per_year = 1L,
    rate,
    type = "effective",
    m = 1L,
    growth = 0,
    growth_type = "effective",
    growth_m = 1L,
    deferral_years = 0,
    timing = c("immediate", "due"),
    payment = 1
) {
  if (!is.numeric(payment)) {
    stop("`payment` must be a numeric vector.", call. = FALSE)
  }

  size <- max(
    length(n_years),
    length(payments_per_year),
    length(rate),
    length(type),
    length(m),
    length(growth),
    length(growth_type),
    length(growth_m),
    length(deferral_years),
    length(timing),
    length(payment)
  )

  valid_size <- function(x) length(x) %in% c(1L, size)

  if (!valid_size(n_years) || !valid_size(payments_per_year) || !valid_size(rate) ||
      !valid_size(type) || !valid_size(m) || !valid_size(growth) ||
      !valid_size(growth_type) || !valid_size(growth_m) ||
      !valid_size(deferral_years) || !valid_size(timing) ||
      !valid_size(payment)) {
    stop(
      "`n_years`, `payments_per_year`, `rate`, `type`, `m`, `growth`, ",
      "`growth_type`, `growth_m`, `deferral_years`, `timing`, and `payment` ",
      "must have length 1 or a common length.",
      call. = FALSE
    )
  }

  n_years           <- rep_len(n_years, size)
  payments_per_year <- rep_len(payments_per_year, size)
  rate              <- rep_len(rate, size)
  type              <- rep_len(type, size)
  m                 <- rep_len(m, size)
  growth            <- rep_len(growth, size)
  growth_type       <- rep_len(growth_type, size)
  growth_m          <- rep_len(growth_m, size)
  deferral_years    <- rep_len(deferral_years, size)
  timing            <- rep_len(timing, size)
  payment           <- rep_len(payment, size)

  bad_payment <- !is.na(payment) & !is.finite(payment)
  if (any(bad_payment)) {
    stop("`payment` must contain only finite numeric values or NA.", call. = FALSE)
  }

  gs_factor <- gs_angle(
    n_years = n_years,
    payments_per_year = payments_per_year,
    rate = rate,
    type = type,
    m = m,
    growth = growth,
    growth_type = growth_type,
    growth_m = growth_m,
    deferral_years = deferral_years,
    timing = timing
  )

  i_effective <- standardize_interest(type = type, rate = rate, m = m)
  g_effective <- standardize_interest(type = growth_type, rate = growth, m = growth_m)

  k <- as.integer(payments_per_year)
  i_period <- (1 + i_effective)^(1 / k) - 1
  g_period <- (1 + g_effective)^(1 / k) - 1

  n_periods <- n_years * payments_per_year
  deferral_periods <- deferral_years * payments_per_year

  tibble::tibble(
    n_years = n_years,
    payments_per_year = as.integer(payments_per_year),
    deferral_years = deferral_years,
    timing = tolower(timing),
    rate_input = rate,
    rate_type = type,
    m = m,
    growth_input = growth,
    growth_type = growth_type,
    growth_m = growth_m,
    i_effective = i_effective,
    g_effective = g_effective,
    i_period = i_period,
    g_period = g_period,
    n_periods = n_periods,
    deferral_periods = deferral_periods,
    gs_factor = gs_factor,
    payment = payment,
    future_value = payment * gs_factor
  )
}
