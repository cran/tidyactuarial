#' Decreasing annuity factor (Da-angle-n)
#'
#' Computes the actuarial present value factor for a decreasing annuity.
#'
#' The payment pattern is arithmetic decreasing. For a unit pattern,
#' the payments over \eqn{n} periods are:
#' \deqn{n, n-1, \ldots, 2, 1.}{n, n-1, ..., 2, 1.}
#'
#' Supported timing conventions:
#' \itemize{
#'   \item \code{"immediate"}: payments at the end of each period.
#'   \item \code{"due"}: payments at the beginning of each period.
#' }
#'
#' For discrete payments with \code{payments_per_year = k}, the total number
#' of payment periods is \eqn{n k}. Deferral is supported through
#' \code{deferral_years}; for discrete annuities, the deferral must align
#' with the payment grid.
#'
#' The input interest rate may be supplied as:
#' \itemize{
#'   \item annual effective interest rate \eqn{i},
#'   \item nominal annual interest rate \eqn{j^{(m)}},
#'   \item nominal annual discount rate \eqn{d^{(m)}},
#'   \item force of interest \eqn{\delta}.
#' }
#'
#' Internally, all rate specifications are first converted to the
#' equivalent annual effective interest rate using \code{\link{standardize_interest}}.
#'
#' @param n_years Numeric vector of payment durations in years.
#'   Each value must be positive and finite.
#' @param payments_per_year Positive integer vector giving the number of
#'   payments per year.
#' @param rate Numeric vector of rate values.
#' @param type Character vector indicating the rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the compounding frequency
#'   for nominal rates.
#' @param deferral_years Numeric vector of deferral times in years.
#'   Must be greater than or equal to 0.
#' @param timing Character vector. One of \code{"immediate"} or \code{"due"}.
#'
#' @return Numeric vector of decreasing annuity factors.
#'
#' @details
#' Let \eqn{i_p} be the effective rate per payment period,
#' \eqn{v_p = (1+i_p)^{-1}}, and let \eqn{n} be the total number of
#' payment periods. Then the decreasing annuity-immediate factor is
#' \deqn{(Da)_{\overline{n|}} = \sum_{k=1}^n (n+1-k)v_p^k.}{(Da)_n| = sum_(k=1)^n (n+1-k)v_p^k.}
#'
#' For \eqn{i_p \neq 0}, a closed-form expression is
#' \deqn{(Da)_{\overline{n|}} = \frac{n}{i_p} - \frac{1}{i_p^2} + \frac{v_p^n}{i_p^2}.}{(Da)_n| = (n)/(i_p) - (1)/(i_p^2) + (v_p^n)/(i_p^2).}
#'
#' For \eqn{i_p = 0}, the limit is
#' \deqn{(Da)_{\overline{n|}} = \frac{n(n+1)}{2}.}{(Da)_n| = (n(n+1))/(2).}
#'
#' The annuity-due version is obtained by multiplying the immediate factor
#' by \eqn{1+i_p}. Deferred versions are obtained by multiplying by
#' \eqn{v_p^h}, where \eqn{h} is the number of deferred periods.
#'
#' A useful identity relates the decreasing and increasing annuity factors:
#' \deqn{(Da)_{\overline{n|}} + (Ia)_{\overline{n|}} = (n+1) \, a_{\overline{n|}}.}{(Da)_n| + (Ia)_n| = (n+1) a_n|.}
#'
#' Input vectors must have length 1 or a common length.
#' Missing values are propagated.
#'
#' @seealso \code{\link{Da_angle_tbl}}, \code{\link{Ia_angle}}, \code{\link{a_angle}},
#'   \code{\link{s_angle}}, \code{\link{standardize_interest}}
#'
#' @family annuities
#'
#' @examples
#' # Simple scalar example
#' Da_angle(
#'   n_years = 10,
#'   rate = 0.05,
#'   type = "effective",
#'   timing = "immediate"
#' )
#'
#' # Medium vectorized example
#' Da_angle(
#'   n_years = c(10, 20),
#'   payments_per_year = c(1, 12),
#'   rate = c(0.05, 0.08),
#'   type = c("effective", "nominal_interest"),
#'   m = c(1, 12),
#'   deferral_years = c(0, 2),
#'   timing = c("immediate", "due")
#' )
#'
#' @export
Da_angle <- function(
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

  # --- Early type validation (before recycling) ---
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

  # --- Determine common size and validate lengths ---
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

  # --- Recycle to common size ---
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
  valid_timing <- c("immediate", "due")
  bad_timing <- !is.na(timing) & !(timing %in% valid_timing)
  if (any(bad_timing)) {
    stop("`timing` must be 'immediate' or 'due'.", call. = FALSE)
  }

  # --- Convert to effective annual rate ---
  i_effective <- standardize_interest(type = type, rate = rate, m = m)

  out <- rep(NA_real_, size)
  eps <- 1e-12

  # --- Identify NA elements ---
  is_na_elem <- is.na(n_years) | is.na(rate) | is.na(i_effective) |
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
    v_period <- 1 / (1 + i_period)

    base <- if (abs(i_period) < eps) {
      n_periods * (n_periods + 1) / 2
    } else {
      n_periods / i_period - 1 / (i_period^2) + v_period^n_periods / (i_period^2)
    }

    if (timing[j] == "due") {
      base <- (1 + i_period) * base
    }

    out[j] <- base * v_period^deferral_periods
  }

  out
}


#' Decreasing annuity details in tibble form
#'
#' Computes the actuarial present value factor for a decreasing annuity and
#' returns a tibble with the main input values, implied rates, decreasing
#' annuity factor, payment scale, and present value.
#'
#' This is a reporting wrapper around \code{\link{Da_angle}}. The decreasing annuity
#' factor assumes a unit pattern. The present value is then computed as
#' \deqn{PV = R \times (Da)}{PV = R * (Da)}
#' where \eqn{R} is the payment scale.
#'
#' @param n_years Numeric vector of payment durations in years.
#' @param payments_per_year Positive integer vector giving the number of
#'   payments per year.
#' @param rate Numeric vector of rate values.
#' @param type Character vector indicating the rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the compounding frequency
#'   for nominal rates.
#' @param deferral_years Numeric vector of deferral times in years.
#' @param timing Character vector. One of \code{"immediate"} or \code{"due"}.
#' @param payment Numeric vector of payment scale factors.
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
#'   \item{i_period}{Equivalent per-payment effective rate.}
#'   \item{v_period}{Equivalent per-payment discount factor.}
#'   \item{n_periods}{Number of payment periods.}
#'   \item{deferral_periods}{Number of deferred periods.}
#'   \item{Da_factor}{Computed decreasing annuity factor.}
#'   \item{payment}{Payment scale factor.}
#'   \item{present_value}{Present value of the decreasing annuity.}
#' }
#'
#' @seealso \code{\link{Da_angle}}, \code{\link{Ia_angle_tbl}}, \code{\link{a_angle_tbl}},
#'   \code{\link{standardize_interest}}
#'
#' @family annuities
#'
#' @examples
#' # Simple scalar example
#' Da_angle_tbl(
#'   n_years = 10,
#'   rate = 0.05,
#'   type = "effective",
#'   payment = 100
#' )
#'
#' # Medium vectorized example
#' Da_angle_tbl(
#'   n_years = c(10, 20),
#'   payments_per_year = c(1, 12),
#'   rate = c(0.05, 0.08),
#'   type = c("effective", "nominal_interest"),
#'   m = c(1, 12),
#'   deferral_years = c(0, 2),
#'   timing = c("immediate", "due"),
#'   payment = c(100, 50)
#' )
#'
#' @export
Da_angle_tbl <- function(
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

  Da_factor <- Da_angle(
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

  k <- as.integer(payments_per_year)
  i_period <- (1 + i_effective)^(1 / k) - 1
  v_period <- 1 / (1 + i_period)
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
    i_effective = i_effective,
    delta = delta,
    i_period = i_period,
    v_period = v_period,
    n_periods = n_periods,
    deferral_periods = deferral_periods,
    Da_factor = Da_factor,
    payment = payment,
    present_value = payment * Da_factor
  )
}
