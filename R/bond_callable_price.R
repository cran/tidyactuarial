#' Price of a callable bond at a target minimum yield
#'
#' Computes the maximum price an investor should pay for a callable bond in
#' order to guarantee a specified minimum yield.
#'
#' The bond is evaluated under each possible redemption scenario:
#' \itemize{
#'   \item each callable date with its associated call price, and
#'   \item final maturity with its final redemption value.
#' }
#'
#' For each scenario, the bond price is computed using the target yield.
#' The callable-bond price returned by this function is the smallest of those
#' scenario prices, that is, the maximum price consistent with the target
#' yield under the least favorable redemption scenario for the investor.
#'
#' This follows the standard actuarial/financial interpretation used in
#' introductory fixed-income mathematics: when a bond is callable at the
#' issuer's option, the investor must protect against the redemption scenario
#' that is least favorable to the investor at the required yield.
#'
#' Assumptions:
#' \itemize{
#'   \item Coupons are paid in arrears at regular intervals.
#'   \item \code{years_to_maturity * coupons_per_year} must be an integer.
#'   \item Each \code{call_times * coupons_per_year} must be an integer.
#'   \item Stub periods are not supported.
#'   \item Pricing is performed at a coupon date; no accrued interest is included.
#' }
#'
#' Yield input conventions:
#' \itemize{
#'   \item If \code{yield_effective_per_period} is supplied, it takes precedence over
#'         \code{yield_rate}, \code{yield_rate_type}, and \code{yield_m}, and is interpreted
#'         as the effective yield per coupon period.
#'   \item Otherwise, \code{yield_rate}, \code{yield_rate_type}, and \code{yield_m} define
#'         an annual yield specification, which is converted first to annual effective
#'         yield and then to effective yield per coupon period.
#' }
#'
#' @param face Numeric scalar. Face or par value of the bond.
#' @param coupon_rate Numeric scalar. Annual coupon rate as a proportion.
#' @param years_to_maturity Numeric scalar. Final maturity in years.
#'   Must be strictly positive.
#' @param coupons_per_year Positive integer. Number of coupon payments per year.
#' @param call_times Numeric vector of callable times in years.
#'   Each value must be strictly between \code{0} and \code{years_to_maturity},
#'   and must align with coupon dates.
#' @param call_prices Numeric vector of call prices corresponding to
#'   \code{call_times}.
#' @param yield_effective_per_period Optional numeric scalar. Effective yield per
#'   coupon period. If supplied, it is used directly.
#' @param yield_rate Optional numeric scalar. Annual yield rate value.
#' @param yield_rate_type Character string indicating the annual yield type:
#'   \code{"effective"}, \code{"nominal_interest"}, \code{"nominal_discount"}, or
#'   \code{"force"}.
#' @param yield_m Positive integer. Compounding frequency for nominal annual yields.
#' @param redemption Numeric scalar. Redemption value at final maturity.
#'   If \code{NULL}, defaults to \code{face}.
#' @param tol Numeric scalar. Tolerance used in alignment checks.
#' @param check Logical scalar. If \code{TRUE}, performs input validation.
#' @param output Character string. Use \code{"value"} to return the worst-case
#'   callable-bond price, or \code{"table"} to return a tibble with all redemption
#'   scenarios.
#'
#' @return
#' If \code{output = "value"}, a numeric scalar: the worst-case callable-bond
#' price consistent with the target yield.
#'
#' If \code{output = "table"}, a tibble with one row per redemption scenario,
#' including scenario prices and the worst-case indicator.
#'
#' @details
#' Let the callable bond have possible redemption scenarios indexed by
#' \eqn{j = 1, \dots, J}, where each scenario corresponds either to a call date
#' or to final maturity. For scenario \eqn{j}, let \eqn{P_j(y)} denote the bond
#' price computed at the target yield \eqn{y} assuming redemption occurs at that
#' scenario time and value.
#'
#' Then this function returns
#' \deqn{\min_j P_j(y)}{min_j P_j(y)}
#' when \code{output = "value"}.
#'
#' This is the maximum price an investor can pay while still guaranteeing
#' at least the target yield under the least favorable redemption scenario.
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
#'   coupon_rate = 0.08,
#'   years_to_maturity = 10,
#'   coupons_per_year = 2,
#'   call_times = c(5, 7),
#'   call_prices = c(105, 102),
#'   yield_rate = 0.06,
#'   yield_rate_type = "effective"
#' )
#'
#' # Table output with all redemption scenarios
#' bond_callable_price(
#'   face = 100,
#'   coupon_rate = 0.08,
#'   years_to_maturity = 10,
#'   coupons_per_year = 2,
#'   call_times = c(5, 7),
#'   call_prices = c(105, 102),
#'   yield_rate = 0.06,
#'   yield_rate_type = "effective",
#'   output = "table"
#' )
#'
#' # Target yield given directly per coupon period
#' bond_callable_price(
#'   face = 1000,
#'   coupon_rate = 0.05,
#'   years_to_maturity = 12,
#'   coupons_per_year = 2,
#'   call_times = c(4, 8),
#'   call_prices = c(1030, 1015),
#'   yield_effective_per_period = 0.028
#' )
#'
#' @export
bond_callable_price <- function(
    face,
    coupon_rate,
    years_to_maturity,
    coupons_per_year = 1L,
    call_times,
    call_prices,
    yield_effective_per_period = NULL,
    yield_rate = NULL,
    yield_rate_type = "effective",
    yield_m = 1L,
    redemption = NULL,
    tol = 1e-10,
    check = TRUE,
    output = c("value", "table")
) {
  output <- match.arg(output)

  out <- .bond_callable_price_scenarios(
    face = face,
    coupon_rate = coupon_rate,
    years_to_maturity = years_to_maturity,
    coupons_per_year = coupons_per_year,
    call_times = call_times,
    call_prices = call_prices,
    y_effective_per_period = yield_effective_per_period,
    y_rate = yield_rate,
    y_type = yield_rate_type,
    y_m = yield_m,
    redemption = redemption,
    tol = tol,
    check = check
  )

  if (output == "value") {
    return(min(out$price_at_target_yield))
  }

  out
}


# Internal helper: builds scenario table and scenario prices
.bond_callable_price_scenarios <- function(
    face,
    coupon_rate,
    years_to_maturity,
    coupons_per_year = 1L,
    call_times,
    call_prices,
    y_effective_per_period = NULL,
    y_rate = NULL,
    y_type = "effective",
    y_m = 1L,
    redemption = NULL,
    tol = 1e-10,
    check = TRUE
) {
  if (is.null(redemption)) {
    redemption <- face
  }

  if (isTRUE(check)) {
    .validate_bond_core(
      face = face,
      coupon_rate = coupon_rate,
      years_to_maturity = years_to_maturity,
      coupons_per_year = coupons_per_year,
      y_m = y_m,
      redemption = redemption,
      tol = tol
    )

    if (years_to_maturity <= 0) {
      stop(
        "`years_to_maturity` must be strictly positive for callable bonds.",
        call. = FALSE
      )
    }

    if (missing(call_times) || missing(call_prices)) {
      stop("`call_times` and `call_prices` must be provided.", call. = FALSE)
    }

    if (!is.numeric(call_times) ||
        any(is.na(call_times)) ||
        any(!is.finite(call_times))) {
      stop("`call_times` must be a finite numeric vector.", call. = FALSE)
    }

    if (!is.numeric(call_prices) ||
        any(is.na(call_prices)) ||
        any(!is.finite(call_prices))) {
      stop("`call_prices` must be a finite numeric vector.", call. = FALSE)
    }

    if (length(call_times) != length(call_prices)) {
      stop("`call_times` and `call_prices` must have the same length.", call. = FALSE)
    }

    if (length(call_times) == 0L) {
      stop("`call_times` must have length at least 1.", call. = FALSE)
    }

    if (any(call_times <= 0) || any(call_times >= years_to_maturity)) {
      stop(
        "`call_times` must be strictly between 0 and `years_to_maturity`.",
        call. = FALSE
      )
    }

    if (any(call_prices < 0)) {
      stop("`call_prices` must be nonnegative.", call. = FALSE)
    }

    if (any(duplicated(call_times))) {
      stop("`call_times` must not contain duplicates.", call. = FALSE)
    }
  }

  coupons_per_year <- as.integer(round(coupons_per_year))
  y_m <- as.integer(round(y_m))

  N_raw <- years_to_maturity * coupons_per_year

  if (abs(N_raw - round(N_raw)) > tol) {
    stop(
      "`years_to_maturity * coupons_per_year` must be an integer ",
      "(stub periods are not supported).",
      call. = FALSE
    )
  }

  call_N_raw <- call_times * coupons_per_year

  if (any(abs(call_N_raw - round(call_N_raw)) > tol)) {
    stop(
      "Each `call_times * coupons_per_year` must be an integer ",
      "(call dates must align with coupon dates).",
      call. = FALSE
    )
  }

  # --- Resolve yield ---
  yield <- .resolve_bond_yield(
    y_effective_per_period = y_effective_per_period,
    y_rate = y_rate,
    y_type = y_type,
    y_m = y_m,
    coupons_per_year = coupons_per_year
  )

  ip <- yield$ip
  i_annual <- yield$i_annual

  # --- Build scenario table ---
  scenario_tbl <- tibble::tibble(
    scenario_id = seq_along(call_times),
    scenario_type = "call",
    scenario_time = as.numeric(call_times),
    redemption_value = as.numeric(call_prices)
  )

  maturity_tbl <- tibble::tibble(
    scenario_id = length(call_times) + 1L,
    scenario_type = "maturity",
    scenario_time = years_to_maturity,
    redemption_value = redemption
  )

  out <- dplyr::bind_rows(scenario_tbl, maturity_tbl) |>
    dplyr::arrange(.data$scenario_time, .data$scenario_id)

  # --- Price each scenario ---
  out$price_at_target_yield <- vapply(
    seq_len(nrow(out)),
    function(j) {
      bond_price(
        face = face,
        coupon_rate = coupon_rate,
        years_to_maturity = out$scenario_time[[j]],
        coupons_per_year = coupons_per_year,
        y_effective_per_period = ip,
        redemption = out$redemption_value[[j]],
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
  out$coupons_per_year <- coupons_per_year
  out$face <- face
  out$coupon_rate <- coupon_rate
  out$years_to_maturity <- years_to_maturity

  out
}
