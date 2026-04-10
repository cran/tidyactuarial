#' Future value of a general cash flow
#'
#' Computes the future value of a cash-flow vector under either:
#' \itemize{
#'   \item a constant interest-rate specification, or
#'   \item a term structure of spot rates, one rate per cash flow.
#' }
#'
#' The cash flow is supplied explicitly through \code{payment}. Its timing is
#' supplied either through \code{time} (in years) or \code{date} (calendar
#' dates). If \code{date} is supplied, the earliest date is taken as time 0.
#'
#' The future value is accumulated to the latest payment time (or latest date).
#'
#' Interest-rate input:
#' \itemize{
#'   \item If \code{rate} has length 1, the same rate is used for all payments.
#'   \item If \code{rate} has the same length as \code{payment}, each rate is
#'         interpreted as the annual effective spot rate associated with the
#'         corresponding payment time.
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
#' \code{\link{standardize_interest}}. When \code{rate} is a vector, the
#' accumulation uses the implied discount-factor ratio under the spot-rate
#' interpretation.
#'
#' @param payment Numeric vector of cash flows.
#' @param rate Numeric scalar or numeric vector of rate values.
#' @param type Character vector indicating the rate type:
#'   \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, or \code{"force"}.
#'   May have length 1 or the same length as \code{payment}.
#' @param m Positive integer vector giving the compounding frequency for
#'   nominal rates. May have length 1 or the same length as \code{payment}.
#' @param time Optional numeric vector of payment times in years.
#' @param date Optional vector of payment dates. If supplied, the earliest
#'   date is treated as time 0.
#' @param day_count Day-count convention used to convert dates to year fractions.
#'   One of \code{"act/365"} or \code{"act/360"}.
#' @param tol Numeric tolerance used in internal checks.
#'
#' @return Numeric scalar: the future value of the cash flow accumulated to the
#'   latest payment time.
#'
#' @details
#' Let \eqn{T} be the latest payment time (the accumulation horizon). For
#' each payment \eqn{C_k} at time \eqn{t_k} with annual effective spot rate
#' \eqn{i_k}, the future value contribution is
#' \deqn{FV_k = C_k \cdot \frac{(1+i_T)^T}{(1+i_k)^{t_k}}}{FV_k = C_k * (1+i_T)^T / (1+i_k)^t_k}
#' and the total future value is \eqn{FV = \sum_k FV_k}{FV = sum(FV_k)}.
#'
#' When a single constant rate is supplied, \eqn{i_k = i_T = i} for all
#' \eqn{k} and the formula simplifies to
#' \eqn{FV = \sum_k C_k (1+i)^{T - t_k}}{FV = sum(C_k * (1+i)^(T - t_k))}.
#'
#' @seealso \code{\link{pv_flow}}, \code{\link{future_value}},
#'   \code{\link{standardize_interest}}
#'
#' @family time-value
#'
#' @examples
#' # Constant annual effective rate
#' fv_flow(
#'   payment = c(100, 150, 200),
#'   rate = 0.08,
#'   type = "effective",
#'   time = c(0, 1, 2)
#' )
#'
#' # Spot rates, one per payment
#' fv_flow(
#'   payment = c(100, 150, 200),
#'   rate = c(0.05, 0.055, 0.06),
#'   type = "effective",
#'   time = c(1, 2, 3)
#' )
#'
#' # Using dates; earliest date is taken as t = 0
#' fv_flow(
#'   payment = c(100, 150, 200),
#'   rate = c(0.05, 0.055, 0.06),
#'   type = "effective",
#'   date = as.Date(c("2026-01-10", "2027-01-10", "2028-01-10"))
#' )
#'
#' # Nominal rates by payment
#' fv_flow(
#'   payment = c(100, 100, 100),
#'   rate = c(0.12, 0.12, 0.12),
#'   type = "nominal_interest",
#'   m = c(12, 12, 12),
#'   time = c(1, 2, 3)
#' )
#'
#' @export
fv_flow <- function(
    payment,
    rate,
    type = "effective",
    m = 1L,
    time = NULL,
    date = NULL,
    day_count = c("act/365", "act/360"),
    tol = 1e-10
) {
  day_count <- match.arg(day_count)

  if (missing(payment)) {
    stop("`payment` must be provided.", call. = FALSE)
  }
  if (missing(rate)) {
    stop("`rate` must be provided.", call. = FALSE)
  }

  # --- Early type validation ---
  if (!is.numeric(payment)) {
    stop("`payment` must be a numeric vector.", call. = FALSE)
  }
  if (length(payment) == 0L) {
    return(0)
  }
  if (any(is.na(payment)) || any(!is.finite(payment))) {
    stop("`payment` must contain only finite numeric values.", call. = FALSE)
  }

  if (!is.numeric(rate)) {
    stop("`rate` must be numeric.", call. = FALSE)
  }
  if (any(is.na(rate)) || any(!is.finite(rate))) {
    stop("`rate` must contain only finite numeric values.", call. = FALSE)
  }

  if (!is.character(type)) {
    stop("`type` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(m)) {
    stop("`m` must be numeric.", call. = FALSE)
  }

  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("`tol` must be a single positive numeric value.", call. = FALSE)
  }

  if (!is.null(time) && !is.null(date)) {
    stop("Provide only one of `time` or `date`, not both.", call. = FALSE)
  }

  n_cf <- length(payment)

  # ---- Build times ----
  if (!is.null(time)) {
    if (!is.numeric(time)) {
      stop("`time` must be a numeric vector.", call. = FALSE)
    }
    if (length(time) != n_cf) {
      stop("`time` and `payment` must have the same length.", call. = FALSE)
    }
    if (any(is.na(time)) || any(!is.finite(time)) || any(time < 0)) {
      stop("`time` must contain only finite values >= 0.", call. = FALSE)
    }
    t <- as.numeric(time)

  } else if (!is.null(date)) {
    date <- as.Date(date)
    if (any(is.na(date))) {
      stop("`date` must contain valid dates.", call. = FALSE)
    }
    if (length(date) != n_cf) {
      stop("`date` and `payment` must have the same length.", call. = FALSE)
    }

    origin <- min(date)
    denom <- switch(
      day_count,
      "act/365" = 365,
      "act/360" = 360
    )
    t <- as.numeric(date - origin) / denom

  } else {
    stop("You must provide either `time` or `date`.", call. = FALSE)
  }

  # ---- Validate compatible lengths for rate inputs ----
  valid_len <- function(x) length(x) %in% c(1L, n_cf)

  if (!valid_len(rate)) {
    stop("`rate` must have length 1 or the same length as `payment`.", call. = FALSE)
  }
  if (!valid_len(type)) {
    stop("`type` must have length 1 or the same length as `payment`.", call. = FALSE)
  }
  if (!valid_len(m)) {
    stop("`m` must have length 1 or the same length as `payment`.", call. = FALSE)
  }

  rate <- rep_len(rate, n_cf)
  type <- rep_len(type, n_cf)
  m    <- rep_len(m, n_cf)

  i_eff <- standardize_interest(type = type, rate = rate, m = m)

  # ---- Accumulate to the latest payment time ----
  # The horizon rate i_T is the spot rate at T. When multiple payments share
  # the horizon time, their spot rates must agree.
  T_horizon <- max(t)

  idx_h <- which(abs(t - T_horizon) <= tol)
  i_h <- i_eff[idx_h]

  if (length(i_h) == 0L) {
    stop("Could not determine the horizon rate.", call. = FALSE)
  }
  if (max(i_h) - min(i_h) > tol) {
    stop(
      "When multiple payments occur at the horizon, their associated spot rates must agree.",
      call. = FALSE
    )
  }

  i_T <- i_h[1]

  sum(payment * ((1 + i_T)^T_horizon / (1 + i_eff)^t))
}
