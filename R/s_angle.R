#' Level annuity accumulation factor s-angle-n
#'
#' Computes the actuarial accumulation factor for a level annuity using compact
#' actuarial notation.
#'
#' Supported timing conventions:
#' \itemize{
#'   \item \code{"immediate"}: annuity-immediate with discrete payments.
#'   \item \code{"due"}: annuity-due with discrete payments.
#'   \item \code{"continuous"}: continuous annuity.
#' }
#'
#' For discrete annuities, \code{k} is the number of payments per year, so
#' payments are made every \eqn{1/k} year. The function returns the accumulation
#' factor, assuming a unit payment at each payment time.
#'
#' Horizon convention:
#' the future value is measured at the time of the last payment. Under this
#' convention, a pure deferment that shifts the entire payment block forward in
#' time does not change the accumulation factor when the payment pattern is
#' otherwise unchanged. Therefore, \code{h} is recorded and validated, but it
#' does not modify the factor.
#'
#' The future value of a perpetuity diverges, so perpetuities are not supported
#' in \code{s_angle()}.
#'
#' @param n Numeric vector of payment durations in years. Each value must be
#'   positive and finite.
#' @param k Positive integer vector giving the number of discrete payments per
#'   year. Ignored for continuous annuities.
#' @param i Numeric vector of interest-rate values.
#' @param i_type Character vector indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the conversion frequency for nominal
#'   rates. Ignored for \code{"effective"} and \code{"force"}.
#' @param h Numeric vector of deferment times in years. Must be greater than or
#'   equal to 0. Under the adopted horizon convention, this is metadata only for
#'   accumulation factors.
#' @param timing Character vector. One of \code{"immediate"}, \code{"due"}, or
#'   \code{"continuous"}.
#' @param payment Numeric vector of level payment amounts. Used only when
#'   \code{tidy = TRUE} to report the corresponding future value. The
#'   accumulation factor itself is always computed for unit payments.
#' @param tidy Logical scalar. If \code{FALSE}, returns a numeric accumulation
#'   factor. If \code{TRUE}, returns a tibble with intermediate calculations.
#'
#' @return
#' If \code{tidy = FALSE}, a numeric vector of accumulation factors.
#'
#' If \code{tidy = TRUE}, a tibble with input values, equivalent rates,
#' accumulation factors, payment amounts, and future values.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}:
#'
#' \itemize{
#'   \item \code{n}: annuity term;
#'   \item \code{k}: payment frequency;
#'   \item \code{i}: interest rate;
#'   \item \code{i_type}: interest-rate type;
#'   \item \code{m}: conversion frequency for nominal rates;
#'   \item \code{h}: deferment period.
#' }
#'
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
#' s_angle(n = 10, i = 0.05)
#'
#' # Nominal interest converted monthly, with monthly payments
#' s_angle(
#'   n = 10,
#'   i = 0.06,
#'   i_type = "nominal_interest",
#'   m = 12,
#'   k = 12
#' )
#'
#' # Continuous annuity
#' s_angle(
#'   n = 15,
#'   i = 0.04,
#'   i_type = "force",
#'   timing = "continuous"
#' )
#'
#' # Tibble output for teaching or auditing
#' s_angle(
#'   n = 10,
#'   i = 0.05,
#'   payment = 1000,
#'   tidy = TRUE
#' )
#'
#' # Vectorized example
#' s_angle(
#'   n = c(5, 10, 20),
#'   k = c(1, 12, 1),
#'   i = c(0.05, 0.06, 0.04),
#'   i_type = c("effective", "nominal_interest", "force"),
#'   m = c(1, 12, 1),
#'   h = c(0, 2, 3),
#'   timing = c("immediate", "due", "continuous")
#' )
#'
#' @export
s_angle <- function(
    n,
    k = 1L,
    i,
    i_type = "effective",
    m = 1L,
    h = 0,
    timing = "immediate",
    payment = 1,
    tidy = FALSE
) {
  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  if (missing(n)) {
    stop("`n` must be provided.", call. = FALSE)
  }
  if (missing(i)) {
    stop("`i` must be provided.", call. = FALSE)
  }

  # --- Early type validation ---
  if (!is.numeric(n)) {
    stop("`n` must be a numeric vector.", call. = FALSE)
  }
  if (!is.numeric(k)) {
    stop("`k` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(i)) {
    stop("`i` must be a numeric vector.", call. = FALSE)
  }
  if (!is.numeric(m)) {
    stop("`m` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(h)) {
    stop("`h` must be numeric.", call. = FALSE)
  }
  if (!is.character(i_type)) {
    stop("`i_type` must be a character vector.", call. = FALSE)
  }
  if (!is.character(timing)) {
    stop("`timing` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(payment)) {
    stop("`payment` must be a numeric vector.", call. = FALSE)
  }

  # --- Determine common size ---
  size <- max(
    length(n),
    length(k),
    length(i),
    length(i_type),
    length(m),
    length(h),
    length(timing),
    length(payment),
    1L
  )

  valid_size <- function(x) length(x) %in% c(1L, size)

  if (!valid_size(n) ||
      !valid_size(k) ||
      !valid_size(i) ||
      !valid_size(i_type) ||
      !valid_size(m) ||
      !valid_size(h) ||
      !valid_size(timing) ||
      !valid_size(payment)) {
    stop(
      "`n`, `k`, `i`, `i_type`, `m`, `h`, `timing`, and `payment` ",
      "must have length 1 or a common length.",
      call. = FALSE
    )
  }

  # --- Recycle ---
  n       <- rep_len(n, size)
  k       <- rep_len(k, size)
  i       <- rep_len(i, size)
  i_type  <- rep_len(i_type, size)
  m       <- rep_len(m, size)
  h       <- rep_len(h, size)
  timing  <- rep_len(timing, size)
  payment <- rep_len(payment, size)

  # --- Value-level validation ---
  bad_n <- !is.na(n) & (!is.finite(n) | n <= 0)
  if (any(bad_n)) {
    stop("`n` must contain only finite values greater than 0 or NA.", call. = FALSE)
  }

  bad_h <- !is.na(h) & (!is.finite(h) | h < 0)

  if (any(bad_h)) {
    stop("`h` must contain only finite values >= 0 or NA.", call. = FALSE)
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
  i_effective <- standardize_interest(i_type = i_type, i = i, m = m)
  delta <- log1p(i_effective)

  accumulation_factor <- rep(NA_real_, size)
  eps <- 1e-12

  # --- Intermediate output columns ---
  i_period <- rep(NA_real_, size)
  v_period <- rep(NA_real_, size)
  n_periods <- rep(NA_real_, size)
  h_periods <- rep(NA_real_, size)

  # --- Identify computable elements ---
  is_na_elem <- is.na(n) |
    is.na(i) |
    is.na(i_effective) |
    is.na(h) |
    is.na(timing)

  computable_idx <- which(!is_na_elem)

  for (idx in computable_idx) {
    if (timing[idx] == "continuous") {
      accumulation_factor[idx] <- if (abs(delta[idx]) < eps) {
        n[idx]
      } else {
        expm1(delta[idx] * n[idx]) / delta[idx]
      }

      next
    }

    kk <- k[idx]

    if (is.na(kk) || !is.finite(kk) || kk < 1 || kk != floor(kk)) {
      stop(
        "`k` must be a positive integer for discrete annuities.",
        call. = FALSE
      )
    }

    kk <- as.integer(kk)

    n_periods_raw <- n[idx] * kk
    n_periods[idx] <- round(n_periods_raw)

    if (abs(n_periods_raw - n_periods[idx]) > 1e-10) {
      stop(
        "For discrete annuities, `n * k` must be an integer.",
        call. = FALSE
      )
    }

    h_periods_raw <- h[idx] * kk
    h_periods[idx] <- round(h_periods_raw)

    if (abs(h_periods_raw - h_periods[idx]) > 1e-10) {
      stop(
        "For discrete annuities, `h * k` must be an integer.",
        call. = FALSE
      )
    }

    i_period[idx] <- (1 + i_effective[idx])^(1 / kk) - 1
    v_period[idx] <- 1 / (1 + i_period[idx])

    base <- if (abs(i_period[idx]) < eps) {
      n_periods[idx]
    } else {
      ((1 + i_period[idx])^n_periods[idx] - 1) / i_period[idx]
    }

    if (timing[idx] == "due") {
      base <- (1 + i_period[idx]) * base
    }

    accumulation_factor[idx] <- base
  }

  if (!tidy) {
    return(accumulation_factor)
  }

  k_out <- as.integer(k)
  k_out[timing == "continuous"] <- NA_integer_

  tibble::tibble(
    n = n,
    k = k_out,
    h = h,
    timing = timing,
    i_input = i,
    i_type = i_type,
    m = m,
    i_effective = i_effective,
    delta = delta,
    i_period = i_period,
    v_period = v_period,
    n_periods = n_periods,
    h_periods = h_periods,
    accumulation_factor = accumulation_factor,
    payment = payment,
    future_value = payment * accumulation_factor
  )
}
