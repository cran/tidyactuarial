#' Price of a level coupon bond from its yield
#'
#' Computes the dirty price of a level coupon bond at time 0 from its yield,
#' using compact actuarial notation.
#'
#' Assumptions:
#' \itemize{
#'   \item Coupons are paid in arrears at regular intervals.
#'   \item \code{n * k} must be an integer.
#'   \item Stub periods are not supported.
#'   \item No accrued interest is considered; the price is evaluated at a coupon date.
#' }
#'
#' The yield may be supplied in either of two ways:
#' \itemize{
#'   \item directly as an effective yield per coupon period through
#'         \code{y_effective_per_period}, or
#'   \item as an annual rate specification through \code{y}, \code{y_type}, and
#'         \code{y_m}.
#' }
#'
#' If an annual rate specification is supplied, it is first converted to the
#' equivalent annual effective yield and then to the effective yield per coupon
#' period.
#'
#' @param face Numeric scalar. Face value of the bond.
#' @param c Numeric scalar. Annual coupon rate as a proportion.
#' @param n Numeric scalar. Time to maturity in years.
#' @param k Positive integer. Number of coupon payments per year.
#' @param y_effective_per_period Optional numeric scalar. Effective yield per
#'   coupon period. If supplied, it is used directly.
#' @param y Optional numeric scalar. Annual yield rate value.
#' @param y_type Character string indicating the annual yield type:
#'   \code{"effective"}, \code{"nominal_interest"}, \code{"nominal_discount"},
#'   or \code{"force"}.
#' @param y_m Positive integer. Conversion frequency for nominal annual yields.
#' @param R Numeric scalar. Redemption value at maturity. If \code{NULL},
#'   defaults to \code{face}.
#' @param tol Numeric scalar. Tolerance used when checking alignment of maturity
#'   with coupon periods.
#' @param check Logical scalar. If \code{TRUE}, performs basic input validation.
#'
#' @return Numeric scalar: dirty price of the bond at time 0.
#'
#' @details
#' This function follows the compact bond notation used in
#' \code{tidyactuarial}: \code{face} is the face value, \code{c} is the annual
#' coupon rate, \code{n} is the time to maturity, \code{k} is the coupon
#' frequency, \code{y} is the annual yield input, and \code{R} is the redemption
#' value.
#'
#' Let \eqn{j} be the effective yield per coupon period, \eqn{k} the number of
#' coupon payments per year, and let \eqn{N = nk} be the total number of coupon
#' periods. With coupon per period \eqn{C = face \cdot c/k} and discount factor
#' \eqn{v = 1/(1+j)}, the price is:
#' \deqn{P = \sum_{r=1}^{N} C_r v^r,}
#' where the sum runs over all coupon and redemption cash flows indexed by
#' coupon period \eqn{r}.
#'
#' @seealso \code{\link{bond_ytm}}, \code{\link{bond_cash_flows}},
#'   \code{\link{bond_duration}}, \code{\link{bond_convexity}},
#'   \code{\link{bond_book_value}}, \code{\link{bond_callable_price}}
#'
#' @family bonds
#'
#' @examples
#' # 5-year annual coupon bond, yield given as annual effective
#' bond_price(
#'   face = 100,
#'   c = 0.08,
#'   n = 5,
#'   k = 1,
#'   y = 0.06,
#'   y_type = "effective"
#' )
#'
#' # 10-year semiannual bond, yield given directly per coupon period
#' bond_price(
#'   face = 1000,
#'   c = 0.05,
#'   n = 10,
#'   k = 2,
#'   y_effective_per_period = 0.03
#' )
#'
#' # Semiannual coupons, nominal annual yield convertible quarterly
#' bond_price(
#'   face = 100,
#'   c = 0.08,
#'   n = 5,
#'   k = 2,
#'   y = 0.06,
#'   y_type = "nominal_interest",
#'   y_m = 4
#' )
#'
#' @export
bond_price <- function(
    face,
    c,
    n,
    k = 1L,
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
    .validate_bond_core(face, c, n, k, y_m, R, tol)
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
  N <- as.integer(round(N_raw))

  if (N == 0L) {
    return(R)
  }

  # --- Resolve yield ---
  yield <- .resolve_bond_yield(
    y_effective_per_period = y_effective_per_period,
    y_rate = y,
    y_type = y_type,
    y_m = y_m,
    coupons_per_year = k
  )
  ip <- yield$ip

  cf_tbl <- bond_cash_flows(
    face = face,
    c = c,
    n = n,
    k = k,
    R = R,
    tol = tol,
    check = FALSE
  )

  # Discount by coupon period, not by time in years.
  cf_periods_raw <- cf_tbl$t * k
  if (any(abs(cf_periods_raw - round(cf_periods_raw)) > tol)) {
    stop("Internal error: coupon times are not aligned with coupon periods.", call. = FALSE)
  }

  cf_periods <- as.integer(round(cf_periods_raw))

  sum(cf_tbl$cf / (1 + ip)^cf_periods)
}
