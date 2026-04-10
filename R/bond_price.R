#' Price of a level coupon bond from its yield
#'
#' Computes the dirty price of a level coupon bond at time 0 from its yield.
#'
#' Assumptions:
#' \itemize{
#'   \item Coupons are paid in arrears at regular intervals.
#'   \item \code{years_to_maturity * coupons_per_year} must be an integer.
#'   \item Stub periods are not supported.
#'   \item No accrued interest is considered; the price is evaluated at a coupon date.
#' }
#'
#' The yield may be supplied in either of two ways:
#' \itemize{
#'   \item directly as an effective yield per coupon period through
#'         \code{y_effective_per_period}, or
#'   \item as an annual rate specification through \code{y_rate}, \code{y_type}, and \code{y_m}.
#' }
#'
#' If an annual rate specification is supplied, it is first converted to the
#' equivalent annual effective yield and then to the effective yield per
#' coupon period.
#'
#' @param face Numeric scalar. Face (par) value of the bond.
#' @param coupon_rate Numeric scalar. Annual coupon rate as a proportion.
#' @param years_to_maturity Numeric scalar. Time to maturity in years.
#' @param coupons_per_year Positive integer. Number of coupon payments per year.
#' @param y_effective_per_period Optional numeric scalar. Effective yield per
#'   coupon period. If supplied, it is used directly.
#' @param y_rate Optional numeric scalar. Annual yield rate value.
#' @param y_type Character string indicating the annual yield type:
#'   \code{"effective"}, \code{"nominal_interest"}, \code{"nominal_discount"}, or \code{"force"}.
#' @param y_m Positive integer. Compounding frequency for nominal annual yields.
#' @param redemption Numeric scalar. Redemption value at maturity.
#'   If \code{NULL}, defaults to \code{face}.
#' @param tol Numeric scalar. Tolerance used when checking alignment of maturity
#'   with coupon periods.
#' @param check Logical scalar. If \code{TRUE}, performs basic input validation.
#'
#' @return Numeric scalar: dirty price of the bond at time 0.
#'
#' @details
#' Let \eqn{j} be the effective yield per coupon period, \eqn{m} the number of
#' coupon payments per year, and \eqn{N = T \times m} the total number of
#' periods. With coupon per period \eqn{C = F r / m} and discount factor
#' \eqn{v = 1/(1+j)}, the price is:
#' \deqn{P = C \cdot a_{\overline{N|}|j} + R \cdot v^N = \sum_{k=1}^{N} C_k \, v^k}{P = C * a_N||j + R * v^N = sum_(k=1)^(N) C_k v^k}
#' where the sum runs over all cash flows (coupons and redemption) indexed by
#' coupon period \eqn{k}.
#'
#' @seealso \code{\link{bond_ytm}}, \code{\link{bond_cash_flows}}, \code{\link{bond_duration}},
#'   \code{\link{bond_convexity}}, \code{\link{bond_book_value}}, \code{\link{bond_callable_price}}
#'
#' @family bonds
#'
#' @examples
#' # 5-year annual coupon bond, yield given as annual effective
#' bond_price(
#'   face = 100,
#'   coupon_rate = 0.08,
#'   years_to_maturity = 5,
#'   coupons_per_year = 1,
#'   y_rate = 0.06,
#'   y_type = "effective"
#' )
#'
#' # 10-year semiannual bond, yield given directly per coupon period
#' bond_price(
#'   face = 1000,
#'   coupon_rate = 0.05,
#'   years_to_maturity = 10,
#'   coupons_per_year = 2,
#'   y_effective_per_period = 0.03
#' )
#'
#' # Semiannual coupons, nominal annual yield convertible quarterly
#' bond_price(
#'   face = 100,
#'   coupon_rate = 0.08,
#'   years_to_maturity = 5,
#'   coupons_per_year = 2,
#'   y_rate = 0.06,
#'   y_type = "nominal_interest",
#'   y_m = 4
#' )
#'
#' @export
bond_price <- function(
    face,
    coupon_rate,
    years_to_maturity,
    coupons_per_year = 1L,
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
    .validate_bond_core(face, coupon_rate, years_to_maturity,
                        coupons_per_year, y_m, redemption, tol)
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
  N <- as.integer(round(N_raw))

  if (N == 0L) {
    return(redemption)
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

  cf_tbl <- bond_cash_flows(
    face = face,
    coupon_rate = coupon_rate,
    years_to_maturity = years_to_maturity,
    coupons_per_year = coupons_per_year,
    redemption = redemption,
    tol = tol,
    check = FALSE
  )

  # Discount by coupon period (not by time in years)
  cf_periods <- as.integer(round(cf_tbl$time * coupons_per_year))
  sum(cf_tbl$cash_flow / (1 + ip)^cf_periods)
}
