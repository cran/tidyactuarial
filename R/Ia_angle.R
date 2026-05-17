#' Increasing annuity factor (Ia-angle-n)
#'
#' Computes the actuarial present value factor for an increasing annuity.
#'
#' The payment pattern is arithmetic increasing. For a unit gradient,
#' the payments over \eqn{n} periods are:
#' \deqn{1, 2, \ldots, n-1, n.}{1, 2, ..., n-1, n.}
#'
#' Supported timing conventions:
#' \itemize{
#'   \item \code{"immediate"}: payments at the end of each period.
#'   \item \code{"due"}: payments at the beginning of each period.
#' }
#'
#' For an annuity-immediate, the payments \eqn{1,2,\ldots,n}
#' occur at times \eqn{1,2,\ldots,n}. For an annuity-due, the payments
#' \eqn{1,2,\ldots,n} occur at times \eqn{0,1,\ldots,n-1}.
#'
#' For discrete payments with \code{payments_per_year = k}, the total number
#' of payment periods is \eqn{n k}. Deferral is supported through
#' \code{deferral_years}; for discrete annuities, the deferral must align
#' with the payment grid.
#'
#' The input interest rate may be supplied as:
#' \itemize{
#'   \item annual effective interest rate \eqn{i},
#'   \item nominal annual interest rate \eqn{j^{(m)}}{j(m)},
#'   \item nominal annual discount rate \eqn{d^{(m)}}{d(m)},
#'   \item force of interest \eqn{\delta}{delta}.
#' }
#'
#' Internally, all rate specifications are first converted to the
#' equivalent annual effective interest rate using
#' \code{\link{standardize_interest}}.
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
#' @return Numeric vector of increasing annuity factors.
#'
#' @details
#' Let \eqn{i_p}{ip} be the effective rate per payment period,
#' \eqn{v_p = (1+i_p)^{-1}}{vp = 1/(1+ip)}, and let \eqn{n} be the total
#' number of payment periods. Also let
#' \deqn{a_n = \frac{1 - v_p^n}{i_p}}{a_n = (1 - vp^n) / ip}
#' and
#' \deqn{\ddot{a}_n = (1+i_p)a_n.}{a_due_n = (1 + ip) * a_n.}
#'
#' Then the increasing annuity-immediate factor is
#' \deqn{
#' (Ia)_n = \frac{\ddot{a}_n - n v_p^n}{i_p}
#' }{
#' (Ia)_n = (a_due_n - n * vp^n) / ip
#' }
#'
#' Equivalently,
#' \deqn{
#' (Ia)_n = \sum_{k=1}^{n} k v_p^k.
#' }{
#' (Ia)_n = sum_{k=1}^n k * vp^k.
#' }
#'
#' For \eqn{i_p = 0}{ip = 0}, the limit is
#' \deqn{(Ia)_n = \frac{n(n+1)}{2}}{(Ia)_n = n(n+1)/2}
#'
#' The annuity-due version is obtained by multiplying the immediate factor
#' by \eqn{1+i_p}{1+ip}. Deferred versions are obtained by multiplying by
#' \eqn{v_p^h}{vp^h}, where \eqn{h} is the number of deferred periods.
#'
#' A useful identity relates the increasing and decreasing annuity factors:
#' \deqn{(Ia)_n + (Da)_n = (n+1) \, a_n}{(Ia)_n + (Da)_n = (n+1) * a_n}
#'
#' Input vectors must have length 1 or a common length.
#' Missing values are propagated.
#'
#' @seealso \code{Ia_angle_tbl}, \code{\link{Da_angle}},
#'   \code{\link{a_angle}}, \code{\link{s_angle}},
#'   \code{\link{standardize_interest}}
#'
#' @family annuities
#'
#' @examples
#' Ia_angle(
#'   n_years = 10,
#'   rate = 0.05,
#'   type = "effective",
#'   timing = "immediate"
#' )
#'
#' Ia_angle(
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
Ia_angle <- function(
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

  if (!valid_size(n_years) ||
      !valid_size(payments_per_year) ||
      !valid_size(rate) ||
      !valid_size(type) ||
      !valid_size(m) ||
      !valid_size(deferral_years) ||
      !valid_size(timing)) {
    stop(
      "`n_years`, `payments_per_year`, `rate`, `type`, `m`, ",
      "`deferral_years`, and `timing` must have length 1 or a common length.",
      call. = FALSE
    )
  }

  # --- Recycle ---
  n_years <- rep_len(n_years, size)
  payments_per_year <- rep_len(payments_per_year, size)
  rate <- rep_len(rate, size)
  type <- rep_len(type, size)
  m <- rep_len(m, size)
  deferral_years <- rep_len(deferral_years, size)
  timing <- rep_len(timing, size)

  # --- Value-level validation ---
  bad_n <- !is.na(n_years) & (!is.finite(n_years) | n_years <= 0)
  if (any(bad_n)) {
    stop(
      "`n_years` must contain only finite values greater than 0 or NA.",
      call. = FALSE
    )
  }

  bad_deferral <- !is.na(deferral_years) &
    (!is.finite(deferral_years) | deferral_years < 0)
  if (any(bad_deferral)) {
    stop(
      "`deferral_years` must contain only finite values >= 0 or NA.",
      call. = FALSE
    )
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

  # --- Identify computable elements ---
  ok <- which(
    !is.na(n_years) &
      !is.na(payments_per_year) &
      !is.na(rate) &
      !is.na(i_effective) &
      !is.na(deferral_years) &
      !is.na(timing)
  )

  for (j in ok) {
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

    i_period <- (1 + i_effective[j])^(1 / k) - 1
    v_period <- 1 / (1 + i_period)

    base <- if (abs(i_period) < eps) {
      n * (n + 1) / 2
    } else {
      a_n <- (1 - v_period^n) / i_period
      a_due_n <- (1 + i_period) * a_n

      (a_due_n - n * v_period^n) / i_period
    }

    if (timing[j] == "due") {
      base <- (1 + i_period) * base
    }

    out[j] <- base * v_period^h
  }

  out
}
