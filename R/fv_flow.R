#' Future value of a general cash flow
#'
#' Computes the future value of a cash-flow vector under either:
#' \itemize{
#'   \item a constant interest-rate specification, or
#'   \item a term structure of spot rates, one rate per cash flow.
#' }
#'
#' The cash flow is supplied explicitly through \code{cf}. Its timing is
#' supplied either through \code{t} (in years) or \code{date} (calendar dates).
#' If \code{date} is supplied, the earliest date is taken as time 0.
#'
#' The future value is accumulated to the latest payment time or latest date.
#'
#' Interest-rate input:
#' \itemize{
#'   \item If \code{i} has length 1, the same rate is used for all cash flows.
#'   \item If \code{i} has the same length as \code{cf}, each rate is
#'         interpreted as the annual effective spot rate associated with the
#'         corresponding cash-flow time.
#' }
#'
#' Rate types may be supplied in FM-style notation:
#' \itemize{
#'   \item annual effective rate \eqn{i},
#'   \item nominal annual interest rate \eqn{j^{(m)}}{j(m)},
#'   \item nominal annual discount rate \eqn{d^{(m)}}{d(m)},
#'   \item force of interest \eqn{\delta}{delta}.
#' }
#'
#' Internally, all supplied rates are converted to annual effective rates using
#' \code{\link{standardize_interest}}. When \code{i} is a vector, the
#' accumulation uses the implied discount-factor ratio under the spot-rate
#' interpretation.
#'
#' @param cf Numeric vector of cash flows.
#' @param i Numeric scalar or numeric vector of interest-rate values.
#' @param i_type Character vector indicating the interest-rate type:
#'   \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, or \code{"force"}.
#'   May have length 1 or the same length as \code{cf}.
#' @param m Positive integer vector giving the conversion frequency for nominal
#'   rates. May have length 1 or the same length as \code{cf}.
#' @param t Optional numeric vector of cash-flow times in years.
#' @param date Optional vector of cash-flow dates. If supplied, the earliest
#'   date is treated as time 0.
#' @param day_count Day-count convention used to convert dates to year
#'   fractions. One of \code{"act/365"} or \code{"act/360"}.
#' @param tol Numeric tolerance used in internal checks.
#'
#' @return Numeric scalar: the future value of the cash flow accumulated to the
#' latest cash-flow time.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{cf} denotes cash flows, \code{t} denotes time,
#' \code{i} denotes the interest rate, \code{i_type} denotes the interest-rate
#' type, and \code{m} denotes the conversion frequency for nominal rates.
#'
#' Let \eqn{T} be the latest cash-flow time, the accumulation horizon. For each
#' cash flow \eqn{C_k} at time \eqn{t_k} with annual effective spot rate
#' \eqn{i_k}, the future value contribution is
#' \deqn{FV_k = C_k \cdot \frac{(1+i_T)^T}{(1+i_k)^{t_k}}.}
#'
#' The total future value is \eqn{FV = \sum_k FV_k}. When a single constant rate
#' is supplied, \eqn{i_k = i_T = i} for all \eqn{k}, and the formula simplifies
#' to
#' \deqn{FV = \sum_k C_k (1+i)^{T-t_k}.}
#'
#' @seealso \code{\link{pv_flow}}, \code{\link{future_value}},
#'   \code{\link{standardize_interest}}
#'
#' @family time-value
#'
#' @examples
#' # Constant annual effective rate
#' fv_flow(
#'   cf = c(100, 150, 200),
#'   i = 0.08,
#'   i_type = "effective",
#'   t = c(0, 1, 2)
#' )
#'
#' # Spot rates, one per cash flow
#' fv_flow(
#'   cf = c(100, 150, 200),
#'   i = c(0.05, 0.055, 0.06),
#'   i_type = "effective",
#'   t = c(1, 2, 3)
#' )
#'
#' # Using dates; earliest date is taken as t = 0
#' fv_flow(
#'   cf = c(100, 150, 200),
#'   i = c(0.05, 0.055, 0.06),
#'   i_type = "effective",
#'   date = as.Date(c("2026-01-10", "2027-01-10", "2028-01-10"))
#' )
#'
#' # Nominal rates by cash flow
#' fv_flow(
#'   cf = c(100, 100, 100),
#'   i = c(0.12, 0.12, 0.12),
#'   i_type = "nominal_interest",
#'   m = c(12, 12, 12),
#'   t = c(1, 2, 3)
#' )
#'
#' @export
fv_flow <- function(
    cf,
    i,
    i_type = "effective",
    m = 1L,
    t = NULL,
    date = NULL,
    day_count = c("act/365", "act/360"),
    tol = 1e-10
) {
  day_count <- match.arg(day_count)

  if (missing(cf)) {
    stop("`cf` must be provided.", call. = FALSE)
  }
  if (missing(i)) {
    stop("`i` must be provided.", call. = FALSE)
  }

  # --- Early type validation ---
  if (!is.numeric(cf)) {
    stop("`cf` must be a numeric vector.", call. = FALSE)
  }
  if (length(cf) == 0L) {
    return(0)
  }
  if (any(is.na(cf)) || any(!is.finite(cf))) {
    stop("`cf` must contain only finite numeric values.", call. = FALSE)
  }

  if (!is.numeric(i)) {
    stop("`i` must be numeric.", call. = FALSE)
  }
  if (any(is.na(i)) || any(!is.finite(i))) {
    stop("`i` must contain only finite numeric values.", call. = FALSE)
  }

  if (!is.character(i_type)) {
    stop("`i_type` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(m)) {
    stop("`m` must be numeric.", call. = FALSE)
  }

  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("`tol` must be a single positive numeric value.", call. = FALSE)
  }

  if (!is.null(t) && !is.null(date)) {
    stop("Provide only one of `t` or `date`, not both.", call. = FALSE)
  }

  n_cf <- length(cf)

  # ---- Build times ----
  if (!is.null(t)) {
    if (!is.numeric(t)) {
      stop("`t` must be a numeric vector.", call. = FALSE)
    }
    if (length(t) != n_cf) {
      stop("`t` and `cf` must have the same length.", call. = FALSE)
    }
    if (any(is.na(t)) || any(!is.finite(t)) || any(t < 0)) {
      stop("`t` must contain only finite values >= 0.", call. = FALSE)
    }
    t_cf <- as.numeric(t)

  } else if (!is.null(date)) {
    date <- as.Date(date)
    if (any(is.na(date))) {
      stop("`date` must contain valid dates.", call. = FALSE)
    }
    if (length(date) != n_cf) {
      stop("`date` and `cf` must have the same length.", call. = FALSE)
    }

    origin <- min(date)
    denom <- switch(
      day_count,
      "act/365" = 365,
      "act/360" = 360
    )
    t_cf <- as.numeric(date - origin) / denom

  } else {
    stop("You must provide either `t` or `date`.", call. = FALSE)
  }

  # ---- Validate compatible lengths for rate inputs ----
  valid_len <- function(x) length(x) %in% c(1L, n_cf)

  if (!valid_len(i)) {
    stop("`i` must have length 1 or the same length as `cf`.", call. = FALSE)
  }
  if (!valid_len(i_type)) {
    stop("`i_type` must have length 1 or the same length as `cf`.", call. = FALSE)
  }
  if (!valid_len(m)) {
    stop("`m` must have length 1 or the same length as `cf`.", call. = FALSE)
  }

  i      <- rep_len(i, n_cf)
  i_type <- rep_len(i_type, n_cf)
  m      <- rep_len(m, n_cf)

  i_eff <- standardize_interest(type = i_type, rate = i, m = m)

  # ---- Accumulate to the latest cash-flow time ----
  # The horizon rate i_T is the spot rate at T. When multiple cash flows share
  # the horizon time, their spot rates must agree.
  T_horizon <- max(t_cf)

  idx_h <- which(abs(t_cf - T_horizon) <= tol)
  i_h <- i_eff[idx_h]

  if (length(i_h) == 0L) {
    stop("Could not determine the horizon rate.", call. = FALSE)
  }
  if (max(i_h) - min(i_h) > tol) {
    stop(
      "When multiple cash flows occur at the horizon, their associated spot rates must agree.",
      call. = FALSE
    )
  }

  i_T <- i_h[1]

  sum(cf * ((1 + i_T)^T_horizon / (1 + i_eff)^t_cf))
}
