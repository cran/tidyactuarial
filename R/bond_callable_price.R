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
#'   \item Pricing is performed at a coupon date (dirty price, no accrued interest).
#' }
#'
#' Yield input conventions:
#' \itemize{
#'   \item If \code{y_effective_per_period} is supplied, it takes precedence over
#'         \code{y_rate}, \code{y_type}, and \code{y_m}, and is interpreted as the effective
#'         yield per coupon period.
#'   \item Otherwise, \code{y_rate}, \code{y_type}, and \code{y_m} define an annual yield
#'         specification, which is converted first to annual effective yield and
#'         then to effective yield per coupon period.
#' }
#'
#' @param face Numeric scalar. Face (par) value of the bond.
#' @param coupon_rate Numeric scalar. Annual coupon rate as a proportion.
#' @param years_to_maturity Numeric scalar. Final maturity in years.
#'   Must be strictly positive (a callable bond needs at least one call date
#'   before maturity).
#' @param coupons_per_year Positive integer. Number of coupon payments per year.
#' @param call_times Numeric vector of callable times in years.
#'   Each value must be strictly between \code{0} and \code{years_to_maturity},
#'   and must align with coupon dates.
#' @param call_prices Numeric vector of call prices corresponding to \code{call_times}.
#' @param y_effective_per_period Optional numeric scalar. Effective yield per
#'   coupon period. If supplied, it is used directly.
#' @param y_rate Optional numeric scalar. Annual yield rate value.
#' @param y_type Character string indicating the annual yield type:
#'   \code{"effective"}, \code{"nominal_interest"}, \code{"nominal_discount"}, or \code{"force"}.
#' @param y_m Positive integer. Compounding frequency for nominal annual yields.
#' @param redemption Numeric scalar. Redemption value at final maturity.
#'   If \code{NULL}, defaults to \code{face}.
#' @param tol Numeric scalar. Tolerance used in alignment checks.
#' @param check Logical scalar. If \code{TRUE}, performs input validation.
#'
#' @return Numeric scalar: the worst-case callable-bond price consistent with
#'   the target yield.
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
#' that is, the smallest price across all redemption scenarios.
#'
#' This is the maximum price an investor can pay while still guaranteeing
#' at least the target yield under the least favorable redemption scenario.
#'
#' @seealso \code{\link{bond_callable_price_tbl}}, \code{\link{bond_price}}, \code{\link{bond_cash_flows}},
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
#'   y_rate = 0.06,
#'   y_type = "effective"
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
#'   y_effective_per_period = 0.028
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
    y_effective_per_period = NULL,
    y_rate = NULL,
    y_type = "effective",
    y_m = 1L,
    redemption = NULL,
    tol = 1e-10,
    check = TRUE
) {
  out <- .bond_callable_price_scenarios(
    face = face,
    coupon_rate = coupon_rate,
    years_to_maturity = years_to_maturity,
    coupons_per_year = coupons_per_year,
    call_times = call_times,
    call_prices = call_prices,
    y_effective_per_period = y_effective_per_period,
    y_rate = y_rate,
    y_type = y_type,
    y_m = y_m,
    redemption = redemption,
    tol = tol,
    check = check
  )

  min(out$price_at_target_yield)
}


#' Callable-bond pricing table at a target minimum yield
#'
#' Computes the price of a callable bond under each possible redemption
#' scenario implied by the call schedule and final maturity.
#'
#' This is a reporting wrapper around the callable-bond pricing logic.
#' It returns one row per scenario and identifies the worst-case scenario
#' for the investor, that is, the scenario producing the smallest price
#' at the target yield.
#'
#' Assumptions:
#' \itemize{
#'   \item Coupons are paid in arrears at regular intervals.
#'   \item \code{years_to_maturity * coupons_per_year} must be an integer.
#'   \item Each \code{call_times * coupons_per_year} must be an integer.
#'   \item Stub periods are not supported.
#'   \item Pricing is performed at a coupon date (dirty price, no accrued interest).
#' }
#'
#' Yield input conventions:
#' \itemize{
#'   \item If \code{y_effective_per_period} is supplied, it takes precedence over
#'         \code{y_rate}, \code{y_type}, and \code{y_m}, and is interpreted as the effective
#'         yield per coupon period.
#'   \item Otherwise, \code{y_rate}, \code{y_type}, and \code{y_m} define an annual yield
#'         specification, which is converted first to annual effective yield and
#'         then to effective yield per coupon period.
#' }
#'
#' @param face Numeric scalar. Face (par) value of the bond.
#' @param coupon_rate Numeric scalar. Annual coupon rate as a proportion.
#' @param years_to_maturity Numeric scalar. Final maturity in years.
#' @param coupons_per_year Positive integer. Number of coupon payments per year.
#' @param call_times Numeric vector of callable times in years.
#' @param call_prices Numeric vector of call prices corresponding to \code{call_times}.
#' @param y_effective_per_period Optional numeric scalar. Effective yield per
#'   coupon period. If supplied, it is used directly.
#' @param y_rate Optional numeric scalar. Annual yield rate value.
#' @param y_type Character string indicating the annual yield type:
#'   \code{"effective"}, \code{"nominal_interest"}, \code{"nominal_discount"}, or \code{"force"}.
#' @param y_m Positive integer. Compounding frequency for nominal annual yields.
#' @param redemption Numeric scalar. Redemption value at final maturity.
#'   If \code{NULL}, defaults to \code{face}.
#' @param tol Numeric scalar. Tolerance used in alignment checks.
#' @param check Logical scalar. If \code{TRUE}, performs input validation.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{scenario_id}{Scenario index.}
#'   \item{scenario_type}{\code{"call"} or \code{"maturity"}.}
#'   \item{scenario_time}{Redemption time in years for the scenario.}
#'   \item{redemption_value}{Call price or final redemption value.}
#'   \item{price_at_target_yield}{Price consistent with the target yield under that scenario.}
#'   \item{is_worst_case}{Logical flag indicating the worst-case scenario(s).}
#'   \item{yield_per_period}{Effective yield per coupon period used in pricing.}
#'   \item{yield_effective_annual}{Equivalent annual effective yield.}
#'   \item{coupons_per_year}{Coupon frequency.}
#'   \item{face}{Face value of the bond.}
#'   \item{coupon_rate}{Annual coupon rate.}
#'   \item{years_to_maturity}{Final maturity in years.}
#' }
#'
#' @details
#' This function evaluates the callable bond under:
#' \itemize{
#'   \item each call date and associated call price, and
#'   \item final maturity and final redemption value.
#' }
#'
#' For each scenario, the price consistent with the target yield is computed.
#' The column \code{is_worst_case} marks the scenario(s) producing the smallest
#' price, which corresponds to the least favorable redemption scenario for
#' the investor.
#'
#' @seealso \code{\link{bond_callable_price}}, \code{\link{bond_price}}, \code{\link{bond_cash_flows}},
#'   \code{\link{bond_book_value}}, \code{\link{bond_ytm}}
#'
#' @family bonds
#'
#' @examples
#' bond_callable_price_tbl(
#'   face = 100,
#'   coupon_rate = 0.08,
#'   years_to_maturity = 10,
#'   coupons_per_year = 2,
#'   call_times = c(5, 7),
#'   call_prices = c(105, 102),
#'   y_rate = 0.06,
#'   y_type = "effective"
#' )
#'
#' bond_callable_price_tbl(
#'   face = 1000,
#'   coupon_rate = 0.05,
#'   years_to_maturity = 12,
#'   coupons_per_year = 2,
#'   call_times = c(4, 8),
#'   call_prices = c(1030, 1015),
#'   y_effective_per_period = 0.028
#' )
#'
#' @export
bond_callable_price_tbl <- function(
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
  .bond_callable_price_scenarios(
    face = face,
    coupon_rate = coupon_rate,
    years_to_maturity = years_to_maturity,
    coupons_per_year = coupons_per_year,
    call_times = call_times,
    call_prices = call_prices,
    y_effective_per_period = y_effective_per_period,
    y_rate = y_rate,
    y_type = y_type,
    y_m = y_m,
    redemption = redemption,
    tol = tol,
    check = check
  )
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
    # Core bond validations (shared utility)
    .validate_bond_core(face, coupon_rate, years_to_maturity,
                        coupons_per_year, y_m, redemption, tol)

    # Callable bonds require strictly positive maturity
    if (years_to_maturity <= 0) {
      stop("`years_to_maturity` must be strictly positive for callable bonds.", call. = FALSE)
    }

    # Call schedule validations
    if (missing(call_times) || missing(call_prices)) {
      stop("`call_times` and `call_prices` must be provided.", call. = FALSE)
    }

    if (!is.numeric(call_times) || any(is.na(call_times)) || any(!is.finite(call_times))) {
      stop("`call_times` must be a finite numeric vector.", call. = FALSE)
    }

    if (!is.numeric(call_prices) || any(is.na(call_prices)) || any(!is.finite(call_prices))) {
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

  # --- Resolve yield (shared utility) ---
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
    dplyr::arrange(scenario_time, scenario_id)

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
