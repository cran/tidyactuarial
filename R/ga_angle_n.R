#' Geometric annuity factor ga-angle-n
#'
#' Computes the actuarial present value factor for a geometric annuity.
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
#' Deferral is supported through \code{deferral_years = h}. For discrete
#' annuities, the deferral must align with the payment grid, that is,
#' \eqn{h k} must be an integer, where \eqn{k} is the number of payments per year.
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
#'   Ignored when \code{perpetuity = TRUE}. If \code{perpetuity = FALSE},
#'   each value must be positive and finite.
#' @param payments_per_year Positive integer vector giving the number of
#'   payments per year.
#' @param rate Numeric vector of annual rate values for discounting.
#' @param type Character vector indicating the annual discount-rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the compounding frequency for
#'   nominal discount-rate inputs.
#' @param growth Numeric vector of annual growth-rate values.
#' @param growth_type Character vector indicating the annual growth-rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param growth_m Positive integer vector giving the compounding frequency for
#'   nominal growth-rate inputs.
#' @param deferral_years Numeric vector of deferral times in years.
#'   Must be greater than or equal to 0.
#' @param timing Character vector. One of \code{"immediate"} or \code{"due"}.
#' @param perpetuity Logical vector. If \code{TRUE}, computes the perpetuity factor.
#'
#' @return Numeric vector of geometric annuity factors.
#'
#' @details
#' Let \eqn{i_p}{ip} be the effective interest rate per payment period,
#' \eqn{g_p}{gp} the effective growth rate per payment period, and
#' \eqn{v_p = (1+i_p)^{-1}}{vp = 1/(1+ip)}.
#'
#' For a finite geometric annuity-immediate with \eqn{n} payment periods:
#' \deqn{ga_n = \sum_{t=1}^{n}\frac{(1+g_p)^{t-1}}{(1+i_p)^t}}{ga_n = sum((1+gp)^(t-1) / (1+ip)^t, t=1..n)}
#'
#' If \eqn{i_p \neq g_p}{ip != gp}, then
#' \deqn{ga_n = \frac{1-\left(\frac{1+g_p}{1+i_p}\right)^n}{i_p-g_p}}{ga_n = (1 - ((1+gp)/(1+ip))^n) / (ip - gp)}
#'
#' If \eqn{i_p = g_p}{ip = gp}, the limiting formula is
#' \deqn{ga_n = \frac{n}{1+i_p}}{ga_n = n / (1+ip)}
#'
#' The annuity-due version is obtained by multiplying the immediate factor by
#' \eqn{1+i_p}{1+ip}. Deferred versions are obtained by multiplying by
#' \eqn{v_p^h}{vp^h}, where \eqn{h} is the number of deferred periods.
#'
#' For perpetuities, convergence requires \eqn{i_p > g_p}{ip > gp}. In that case:
#' \deqn{ga_\infty = \frac{1}{i_p-g_p}}{ga_inf = 1 / (ip - gp)}
#' and for due:
#' \deqn{\ddot{ga}_\infty = \frac{1+i_p}{i_p-g_p}}{ga_inf_due = (1+ip) / (ip - gp)}
#'
#' Input vectors must have length 1 or a common length.
#' Missing values are propagated.
#'
#' @seealso \code{\link{ga_angle_tbl}}, \code{\link{a_angle}},
#'   \code{\link{gs_angle}}, \code{\link{standardize_interest}}
#'
#' @family annuities
#'
#' @examples
#' # Simple scalar example
#' ga_angle(
#'   n_years = 10,
#'   rate = 0.05,
#'   type = "effective",
#'   growth = 0.02,
#'   growth_type = "effective",
#'   timing = "immediate"
#' )
#'
#' # Growth = 0 collapses to a level annuity
#' ga_angle(
#'   n_years = 10,
#'   rate = 0.05,
#'   type = "effective",
#'   growth = 0,
#'   growth_type = "effective",
#'   timing = "immediate"
#' )
#'
#' # Medium vectorized example
#' ga_angle(
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
ga_angle <- function(
    n_years = NULL,
    payments_per_year = 1L,
    rate,
    type = "effective",
    m = 1L,
    growth = 0,
    growth_type = "effective",
    growth_m = 1L,
    deferral_years = 0,
    timing = c("immediate", "due"),
    perpetuity = FALSE
) {
  # --- Early type validation ---
  if (!is.null(n_years) && !is.numeric(n_years)) {
    stop("`n_years` must be numeric or NULL.", call. = FALSE)
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
  if (!is.logical(perpetuity)) {
    stop("`perpetuity` must be a logical vector.", call. = FALSE)
  }

  # --- Determine common size ---
  size <- max(
    if (is.null(n_years)) 0L else length(n_years),
    length(payments_per_year),
    length(rate),
    length(type),
    length(m),
    length(growth),
    length(growth_type),
    length(growth_m),
    length(deferral_years),
    length(timing),
    length(perpetuity),
    1L
  )

  valid_size <- function(x) is.null(x) || length(x) %in% c(1L, size)

  if (!valid_size(n_years) || !valid_size(payments_per_year) || !valid_size(rate) ||
      !valid_size(type) || !valid_size(m) || !valid_size(growth) ||
      !valid_size(growth_type) || !valid_size(growth_m) ||
      !valid_size(deferral_years) || !valid_size(timing) ||
      !valid_size(perpetuity)) {
    stop(
      "`n_years`, `payments_per_year`, `rate`, `type`, `m`, `growth`, ",
      "`growth_type`, `growth_m`, `deferral_years`, `timing`, and `perpetuity` ",
      "must have length 1 or a common length.",
      call. = FALSE
    )
  }

  # --- Recycle ---
  if (is.null(n_years)) {
    n_years <- rep(NA_real_, size)
  } else {
    n_years <- rep_len(n_years, size)
  }

  payments_per_year <- rep_len(payments_per_year, size)
  rate              <- rep_len(rate, size)
  type              <- rep_len(type, size)
  m                 <- rep_len(m, size)
  growth            <- rep_len(growth, size)
  growth_type       <- rep_len(growth_type, size)
  growth_m          <- rep_len(growth_m, size)
  deferral_years    <- rep_len(deferral_years, size)
  timing            <- rep_len(timing, size)
  perpetuity        <- rep_len(perpetuity, size)

  # --- Value-level validation ---
  bad_n <- !is.na(n_years) & (!is.finite(n_years) | n_years <= 0)
  if (any(bad_n & !perpetuity)) {
    stop(
      "When `perpetuity = FALSE`, `n_years` must contain only finite values greater than 0.",
      call. = FALSE
    )
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

  if (any(is.na(perpetuity))) {
    stop("`perpetuity` must not contain NA.", call. = FALSE)
  }

  # --- Convert rates ---
  i_effective <- standardize_interest(type = type, rate = rate, m = m)
  g_effective <- standardize_interest(type = growth_type, rate = growth, m = growth_m)

  out <- rep(NA_real_, size)
  eps <- 1e-12

  # --- Identify computable elements ---
  is_na_elem <- is.na(rate) | is.na(growth) | is.na(i_effective) |
    is.na(g_effective) | is.na(deferral_years) | is.na(timing)
  disc_idx <- which(!is_na_elem)

  for (j in disc_idx) {
    k <- payments_per_year[j]
    if (is.na(k) || !is.finite(k) || k < 1 || k != floor(k)) {
      stop("`payments_per_year` must be a positive integer.", call. = FALSE)
    }
    k <- as.integer(k)

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
    v_period <- 1 / (1 + i_period)

    if (perpetuity[j]) {
      if (i_period <= g_period + eps) {
        stop(
          "For a geometric perpetuity, the per-period interest rate must be strictly greater than the per-period growth rate.",
          call. = FALSE
        )
      }

      base <- 1 / (i_period - g_period)

      if (timing[j] == "due") {
        base <- (1 + i_period) * base
      }

      out[j] <- base * v_period^deferral_periods
      next
    }

    n_periods_raw <- n_years[j] * k
    n_periods <- round(n_periods_raw)
    if (abs(n_periods_raw - n_periods) > 1e-10) {
      stop(
        "For discrete annuities, `n_years * payments_per_year` must be an integer.",
        call. = FALSE
      )
    }

    if (abs(i_period - g_period) < eps) {
      base <- n_periods / (1 + i_period)
    } else {
      q <- (1 + g_period) / (1 + i_period)
      base <- (1 - q^n_periods) / (i_period - g_period)
    }

    if (timing[j] == "due") {
      base <- (1 + i_period) * base
    }

    out[j] <- base * v_period^deferral_periods
  }

  out
}


#' Geometric annuity details in tibble form
#'
#' Computes the actuarial present value factor for a geometric annuity and
#' returns a tibble with the main input values, implied rates, factor,
#' payment amount, and present value.
#'
#' This is a reporting wrapper around \code{\link{ga_angle}}. The factor assumes
#' that the first payment equals 1. The actual present value is then computed as
#' \deqn{PV = P_1 \times ga}{PV = P1 * ga}
#' where \eqn{P_1}{P1} is the first payment of the geometric sequence.
#'
#' @param n_years Numeric vector of payment durations in years.
#' @param payments_per_year Positive integer vector giving the number of
#'   payments per year.
#' @param rate Numeric vector of annual rate values for discounting.
#' @param type Character vector indicating the annual discount-rate type.
#' @param m Positive integer vector giving the compounding frequency for
#'   nominal discount-rate inputs.
#' @param growth Numeric vector of annual growth-rate values.
#' @param growth_type Character vector indicating the annual growth-rate type.
#' @param growth_m Positive integer vector giving the compounding frequency for
#'   nominal growth-rate inputs.
#' @param deferral_years Numeric vector of deferral times in years.
#' @param timing Character vector. One of \code{"immediate"} or \code{"due"}.
#' @param payment Numeric vector giving the first payment of the geometric sequence.
#' @param perpetuity Logical vector. If \code{TRUE}, computes the perpetuity factor.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{n_years}{Payment duration in years.}
#'   \item{payments_per_year}{Number of payments per year.}
#'   \item{deferral_years}{Deferral period in years.}
#'   \item{timing}{Payment timing convention.}
#'   \item{perpetuity}{Whether the annuity is perpetual.}
#'   \item{rate_input}{Original supplied discount rate.}
#'   \item{rate_type}{Type of supplied discount rate.}
#'   \item{m}{Compounding frequency for nominal discount-rate inputs.}
#'   \item{growth_input}{Original supplied growth rate.}
#'   \item{growth_type}{Type of supplied growth rate.}
#'   \item{growth_m}{Compounding frequency for nominal growth-rate inputs.}
#'   \item{i_effective}{Equivalent annual effective discount rate.}
#'   \item{g_effective}{Equivalent annual effective growth rate.}
#'   \item{i_period}{Equivalent per-payment discount rate.}
#'   \item{g_period}{Equivalent per-payment growth rate.}
#'   \item{v_period}{Equivalent per-payment discount factor.}
#'   \item{n_periods}{Number of payment periods for finite annuities.}
#'   \item{deferral_periods}{Number of deferred periods.}
#'   \item{ga_factor}{Computed geometric annuity factor.}
#'   \item{payment}{First payment of the geometric sequence.}
#'   \item{present_value}{Present value of the geometric annuity.}
#' }
#'
#' @seealso \code{\link{ga_angle}}, \code{\link{a_angle_tbl}},
#'   \code{\link{standardize_interest}}
#'
#' @family annuities
#'
#' @examples
#' # Simple scalar example
#' ga_angle_tbl(
#'   n_years = 10,
#'   rate = 0.05,
#'   type = "effective",
#'   growth = 0.02,
#'   growth_type = "effective",
#'   payment = 100
#' )
#'
#' # Medium vectorized example
#' ga_angle_tbl(
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
ga_angle_tbl <- function(
    n_years = NULL,
    payments_per_year = 1L,
    rate,
    type = "effective",
    m = 1L,
    growth = 0,
    growth_type = "effective",
    growth_m = 1L,
    deferral_years = 0,
    timing = c("immediate", "due"),
    payment = 1,
    perpetuity = FALSE
) {
  if (!is.numeric(payment)) {
    stop("`payment` must be a numeric vector.", call. = FALSE)
  }

  # --- Determine common size (fixed: use length, not value, for n_years) ---
  size <- max(
    if (is.null(n_years)) 0L else length(n_years),
    length(payments_per_year),
    length(rate),
    length(type),
    length(m),
    length(growth),
    length(growth_type),
    length(growth_m),
    length(deferral_years),
    length(timing),
    length(payment),
    length(perpetuity),
    1L
  )

  valid_size <- function(x) is.null(x) || length(x) %in% c(1L, size)

  if (!valid_size(n_years) || !valid_size(payments_per_year) || !valid_size(rate) ||
      !valid_size(type) || !valid_size(m) || !valid_size(growth) ||
      !valid_size(growth_type) || !valid_size(growth_m) ||
      !valid_size(deferral_years) || !valid_size(timing) ||
      !valid_size(payment) || !valid_size(perpetuity)) {
    stop(
      "`n_years`, `payments_per_year`, `rate`, `type`, `m`, `growth`, ",
      "`growth_type`, `growth_m`, `deferral_years`, `timing`, `payment`, ",
      "and `perpetuity` must have length 1 or a common length.",
      call. = FALSE
    )
  }

  if (is.null(n_years)) {
    n_years <- rep(NA_real_, size)
  } else {
    n_years <- rep_len(n_years, size)
  }

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
  perpetuity        <- rep_len(perpetuity, size)

  bad_payment <- !is.na(payment) & !is.finite(payment)
  if (any(bad_payment)) {
    stop("`payment` must contain only finite numeric values or NA.", call. = FALSE)
  }

  ga_factor <- ga_angle(
    n_years = n_years,
    payments_per_year = payments_per_year,
    rate = rate,
    type = type,
    m = m,
    growth = growth,
    growth_type = growth_type,
    growth_m = growth_m,
    deferral_years = deferral_years,
    timing = timing,
    perpetuity = perpetuity
  )

  i_effective <- standardize_interest(type = type, rate = rate, m = m)
  g_effective <- standardize_interest(type = growth_type, rate = growth, m = growth_m)

  k <- as.integer(payments_per_year)
  i_period <- (1 + i_effective)^(1 / k) - 1
  g_period <- (1 + g_effective)^(1 / k) - 1
  v_period <- 1 / (1 + i_period)

  n_periods <- rep(NA_real_, size)
  finite_idx <- !perpetuity
  n_periods[finite_idx] <- n_years[finite_idx] * payments_per_year[finite_idx]

  deferral_periods <- deferral_years * payments_per_year

  tibble::tibble(
    n_years = n_years,
    payments_per_year = as.integer(payments_per_year),
    deferral_years = deferral_years,
    timing = tolower(timing),
    perpetuity = perpetuity,
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
    v_period = v_period,
    n_periods = n_periods,
    deferral_periods = deferral_periods,
    ga_factor = ga_factor,
    payment = payment,
    present_value = payment * ga_factor
  )
}
