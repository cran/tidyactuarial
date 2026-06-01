#' Price of a callable bond at a target minimum yield
#'
#' Computes the maximum price an investor should pay for a callable bond in
#' order to guarantee a specified minimum yield, using compact actuarial
#' notation.
#'
#' The bond is evaluated under each possible redemption scenario:
#' \itemize{
#'   \item each callable date with its associated call price, and
#'   \item final maturity with its final redemption value.
#' }
#'
#' For each scenario, the bond price is computed using the target yield. The
#' callable-bond price returned by this function is the smallest of those
#' scenario prices, that is, the maximum price consistent with the target yield
#' under the least favorable redemption scenario for the investor.
#'
#' This follows the standard actuarial/financial interpretation used in
#' introductory fixed-income mathematics: when a bond is callable at the
#' issuer's option, the investor must protect against the redemption scenario
#' that is least favorable to the investor at the required yield.
#'
#' Assumptions:
#' \itemize{
#'   \item Coupons are paid in arrears at regular intervals.
#'   \item \code{n * k} must be an integer.
#'   \item Each \code{call_t * k} must be an integer.
#'   \item Stub periods are not supported.
#'   \item Pricing is performed at a coupon date; no accrued interest is included.
#' }
#'
#' Yield input conventions:
#' \itemize{
#'   \item If \code{y_effective_per_period} is supplied, it takes precedence over
#'         \code{y}, \code{y_type}, and \code{y_m}, and is interpreted as the
#'         effective yield per coupon period.
#'   \item Otherwise, \code{y}, \code{y_type}, and \code{y_m} define an annual
#'         yield specification, which is converted first to annual effective
#'         yield and then to effective yield per coupon period.
#' }
#'
#' @param face Numeric scalar. Face or par value of the bond.
#' @param c Numeric scalar. Annual coupon rate as a proportion.
#' @param n Numeric scalar. Final maturity in years. Must be strictly positive.
#' @param k Positive integer. Number of coupon payments per year.
#' @param call_t Numeric vector of callable times in years. Each value must be
#'   strictly between \code{0} and \code{n}, and must align with coupon dates.
#' @param call_R Numeric vector of call prices corresponding to \code{call_t}.
#' @param y_effective_per_period Optional numeric scalar. Effective yield per
#'   coupon period. If supplied, it is used directly.
#' @param y Optional numeric scalar. Annual yield rate value.
#' @param y_type Character string indicating the annual yield type:
#'   \code{"effective"}, \code{"nominal_interest"}, \code{"nominal_discount"},
#'   or \code{"force"}.
#' @param y_m Positive integer. Conversion frequency for nominal annual yields.
#' @param R Numeric scalar. Redemption value at final maturity. If \code{NULL},
#'   defaults to \code{face}.
#' @param tol Numeric scalar. Tolerance used in alignment checks.
#' @param check Logical scalar. If \code{TRUE}, performs input validation.
#' @param tidy Logical scalar. If \code{FALSE}, returns the worst-case
#'   callable-bond price. If \code{TRUE}, returns a tibble with all redemption
#'   scenarios.
#'
#' @return
#' If \code{tidy = FALSE}, a numeric scalar: the worst-case callable-bond price
#' consistent with the target yield.
#'
#' If \code{tidy = TRUE}, a tibble with one row per redemption scenario,
#' including scenario prices and the worst-case indicator.
#'
#' @details
#' This function follows the compact bond notation used in
#' \code{tidyactuarial}: \code{face} is the face value, \code{c} is the annual
#' coupon rate, \code{n} is maturity, \code{k} is coupon frequency, \code{y} is
#' the target yield, \code{R} is the final redemption value, \code{call_t} is
#' the vector of call times, and \code{call_R} is the vector of call prices.
#'
#' Let the callable bond have possible redemption scenarios indexed by
#' \eqn{j = 1, \dots, J}, where each scenario corresponds either to a call date
#' or to final maturity. For scenario \eqn{j}, let \eqn{P_j(y)} denote the bond
#' price computed at the target yield \eqn{y} assuming redemption occurs at that
#' scenario time and value.
#'
#' Then this function returns
#' \deqn{\min_j P_j(y)}{min_j P_j(y)}
#' when \code{tidy = FALSE}.
#'
#' This is the maximum price an investor can pay while still guaranteeing at
#' least the target yield under the least favorable redemption scenario.
#'
#' @seealso \code{\link{bond_price}}, \code{\link{bond_cash_flows}},
#'   \code{\link{bond_book_value}}, \code{\link{bond_ytm}}
#'
#' @family bonds
#'
#' @examples
#' # Callable bond with two possible call dates
#' bond_callable_price(
#'   face = 100,
#'   c = 0.08,
#'   n = 10,
#'   k = 2,
#'   call_t = c(5, 7),
#'   call_R = c(105, 102),
#'   y = 0.06,
#'   y_type = "effective"
#' )
#'
#' # Tidy output with all redemption scenarios
#' bond_callable_price(
#'   face = 100,
#'   c = 0.08,
#'   n = 10,
#'   k = 2,
#'   call_t = c(5, 7),
#'   call_R = c(105, 102),
#'   y = 0.06,
#'   y_type = "effective",
#'   tidy = TRUE
#' )
#'
#' # Target yield given directly per coupon period
#' bond_callable_price(
#'   face = 1000,
#'   c = 0.05,
#'   n = 12,
#'   k = 2,
#'   call_t = c(4, 8),
#'   call_R = c(1030, 1015),
#'   y_effective_per_period = 0.028
#' )
#'
#' @export
bond_callable_price <- function(
    face,
    c,
    n,
    k = 1L,
    call_t,
    call_R,
    y_effective_per_period = NULL,
    y = NULL,
    y_type = "effective",
    y_m = 1L,
    R = NULL,
    tol = 1e-10,
    check = TRUE,
    tidy = FALSE
) {
  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  out <- .bond_callable_price_scenarios(
    face = face,
    c = c,
    n = n,
    k = k,
    call_t = call_t,
    call_R = call_R,
    y_effective_per_period = y_effective_per_period,
    y = y,
    y_type = y_type,
    y_m = y_m,
    R = R,
    tol = tol,
    check = check
  )

  if (!tidy) {
    return(min(out$price_at_target_yield))
  }

  out
}


# Internal helper: builds scenario table and scenario prices
.bond_callable_price_scenarios <- function(
    face,
    c,
    n,
    k = 1L,
    call_t,
    call_R,
    y_effective_per_period = NULL,
    y = NULL,
    y_type = "effective",
    y_m = 1L,
    R = NULL,
    tol = 1e-10,
    check = TRUE
) {
  if (is.null(R)) {
    R <- face
  }

  if (isTRUE(check)) {
    .validate_bond_core(
      face = face,
      c = c,
      n = n,
      k = k,
      y_m = y_m,
      R = R,
      tol = tol
    )

    if (n <= 0) {
      stop(
        "`n` must be strictly positive for callable bonds.",
        call. = FALSE
      )
    }

    if (missing(call_t) || missing(call_R)) {
      stop("`call_t` and `call_R` must be provided.", call. = FALSE)
    }

    if (!is.numeric(call_t) ||
        any(is.na(call_t)) ||
        any(!is.finite(call_t))) {
      stop("`call_t` must be a finite numeric vector.", call. = FALSE)
    }

    if (!is.numeric(call_R) ||
        any(is.na(call_R)) ||
        any(!is.finite(call_R))) {
      stop("`call_R` must be a finite numeric vector.", call. = FALSE)
    }

    if (length(call_t) != length(call_R)) {
      stop("`call_t` and `call_R` must have the same length.", call. = FALSE)
    }

    if (length(call_t) == 0L) {
      stop("`call_t` must have length at least 1.", call. = FALSE)
    }

    if (any(call_t <= 0) || any(call_t >= n)) {
      stop(
        "`call_t` must be strictly between 0 and `n`.",
        call. = FALSE
      )
    }

    if (any(call_R < 0)) {
      stop("`call_R` must be nonnegative.", call. = FALSE)
    }

    if (any(duplicated(call_t))) {
      stop("`call_t` must not contain duplicates.", call. = FALSE)
    }
  }

  k <- as.integer(round(k))
  y_m <- as.integer(round(y_m))

  N_raw <- n * k

  if (abs(N_raw - round(N_raw)) > tol) {
    stop(
      "`n * k` must be an integer ",
      "(stub periods are not supported).",
      call. = FALSE
    )
  }

  call_N_raw <- call_t * k

  if (any(abs(call_N_raw - round(call_N_raw)) > tol)) {
    stop(
      "Each `call_t * k` must be an integer ",
      "(call dates must align with coupon dates).",
      call. = FALSE
    )
  }

  # --- Resolve yield ---
  yield <- .resolve_bond_yield(
    y_effective_per_period = y_effective_per_period,
    y = y,
    y_type = y_type,
    y_m = y_m,
    k = k
  )

  ip <- yield$ip
  i_annual <- yield$i_annual

  # --- Build scenario table ---
  scenario_tbl <- tibble::tibble(
    scenario_id = seq_along(call_t),
    scenario_type = "call",
    t = as.numeric(call_t),
    redemption_value = as.numeric(call_R)
  )

  maturity_tbl <- tibble::tibble(
    scenario_id = length(call_t) + 1L,
    scenario_type = "maturity",
    t = n,
    redemption_value = R
  )

  out <- dplyr::bind_rows(scenario_tbl, maturity_tbl)

  out <- out[order(out$t, out$scenario_id), , drop = FALSE]
  rownames(out) <- NULL

  # --- Price each scenario ---
  out$price_at_target_yield <- vapply(
    seq_len(nrow(out)),
    function(j) {
      bond_price(
        face = face,
        c = c,
        n = out$t[[j]],
        k = k,
        y_effective_per_period = ip,
        R = out$redemption_value[[j]],
        check = FALSE,
        tol = tol
      )
    },
    numeric(1L)
  )

  worst_price <- min(out$price_at_target_yield)

  out$is_worst_case <- abs(out$price_at_target_yield - worst_price) <= sqrt(tol)

  out$yield_per_period <- ip
  out$yield_effective_annual <- i_annual
  out$k <- k
  out$face <- face
  out$c <- c
  out$n <- n
  out$R <- R

  tibble::as_tibble(out)
}
