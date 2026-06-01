#' Macaulay and modified duration of a level coupon bond under a flat yield
#'
#' Computes Macaulay duration and modified-duration measures for a level coupon
#' bond valued under a flat yield-to-maturity assumption, using compact
#' actuarial notation.
#'
#' Assumptions:
#' \itemize{
#'   \item Coupons are paid in arrears at regular intervals.
#'   \item \code{n * k} must be an integer.
#'   \item Stub periods are not supported.
#'   \item Valuation is at a coupon date with no accrued interest.
#'   \item A single flat yield is used to discount all cash flows.
#' }
#'
#' Yield input conventions:
#' \itemize{
#'   \item If \code{y_effective_per_period} is supplied, it is interpreted as the
#'         effective yield per coupon period.
#'   \item Otherwise, \code{y}, \code{y_type}, and \code{y_m} define an annual
#'         yield specification, which is converted first to annual effective
#'         yield and then to effective yield per coupon period.
#' }
#'
#' @param face Numeric scalar. Face value of the bond.
#' @param c Numeric scalar. Annual coupon rate as a proportion.
#' @param n Numeric scalar. Time to maturity in years.
#' @param k Positive integer. Number of coupon payments per year.
#' @param y_effective_per_period Optional numeric scalar. Effective yield per
#'   coupon period.
#' @param y Optional numeric scalar. Annual yield rate value.
#' @param y_type Character string indicating the annual yield type:
#'   \code{"effective"}, \code{"nominal_interest"}, \code{"nominal_discount"},
#'   or \code{"force"}.
#' @param y_m Positive integer. Conversion frequency for nominal annual yields.
#' @param R Numeric scalar. Redemption value at maturity. If \code{NULL},
#'   defaults to \code{face}.
#' @param tol Numeric scalar. Tolerance used to check maturity alignment.
#' @param check Logical scalar. If \code{TRUE}, performs input validation.
#'
#' @return A one-row tibble with:
#' \describe{
#'   \item{price}{Dirty price at the given yield.}
#'   \item{macaulay_duration_periods}{Macaulay duration in coupon periods.}
#'   \item{macaulay_duration_years}{Macaulay duration in years.}
#'   \item{modified_duration_periods_j}{Modified duration with respect to the
#'   effective yield per coupon period, expressed in coupon periods.}
#'   \item{modified_duration_years_i}{Modified duration with respect to the
#'   annual effective yield, expressed in years.}
#'   \item{yield_per_period}{Effective yield per coupon period.}
#'   \item{yield_effective_annual}{Annual effective yield.}
#'   \item{k}{Coupon frequency.}
#'   \item{n_periods}{Total number of coupon periods.}
#' }
#'
#' @details
#' This function follows the compact bond notation used in
#' \code{tidyactuarial}: \code{face} is the face value, \code{c} is the annual
#' coupon rate, \code{n} is the time to maturity, \code{k} is the coupon
#' frequency, \code{y} is the annual yield input, and \code{R} is the redemption
#' value.
#'
#' Let \eqn{j} be the effective yield per coupon period, \eqn{k} the number of
#' coupon payments per year, and let cash flows \eqn{C_r} occur at coupon
#' periods \eqn{r = 1, \dots, N}. With \eqn{v = 1/(1+j)} and
#' \eqn{P = \sum_r C_r v^r}:
#'
#' Macaulay duration in coupon periods:
#' \deqn{D_p =
#' \frac{\sum_{r=1}^{N} r C_r v^r}{P}.}
#'
#' Macaulay duration in years is \eqn{D = D_p / k}.
#'
#' Modified duration with respect to \eqn{j}, in coupon periods, is
#' \eqn{D^*_j = D_p / (1 + j)}.
#'
#' Modified duration with respect to the annual effective rate \eqn{i}, in
#' years, is \eqn{D^*_i = D / (1 + i)}, where
#' \eqn{i = (1+j)^k - 1}.
#'
#' Modified duration measures the first-order sensitivity of the bond price to
#' yield changes. Together with \code{\link{bond_convexity}}, it forms the
#' second-order Taylor approximation of price changes.
#'
#' @seealso \code{\link{bond_convexity}}, \code{\link{bond_price}},
#'   \code{\link{bond_cash_flows}}, \code{\link{bond_book_value}},
#'   \code{\link{bond_ytm}}
#'
#' @family bonds
#'
#' @examples
#' bond_duration(
#'   face = 100,
#'   c = 0.08,
#'   n = 5,
#'   k = 2,
#'   y = 0.06,
#'   y_type = "effective"
#' )
#'
#' bond_duration(
#'   face = 1000,
#'   c = 0.05,
#'   n = 10,
#'   k = 2,
#'   y_effective_per_period = 0.03
#' )
#'
#' @export
bond_duration <- function(
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
  if (!is.logical(check) || length(check) != 1L || is.na(check)) {
    stop("`check` must be TRUE or FALSE.", call. = FALSE)
  }

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
    return(tibble::tibble(
      price = R,
      macaulay_duration_periods = 0,
      macaulay_duration_years = 0,
      modified_duration_periods_j = 0,
      modified_duration_years_i = 0,
      yield_per_period = NA_real_,
      yield_effective_annual = NA_real_,
      k = k,
      n_periods = 0L
    ))
  }

  # --- Resolve yield ---
  yield <- .resolve_bond_yield(
    y_effective_per_period = y_effective_per_period,
    y_rate = y,
    y_type = y_type,
    y_m = y_m,
    coupons_per_year = k
  )
  j <- yield$ip
  i_annual <- yield$i_annual

  # Extra robustness: validate derived yield
  if (!is.finite(i_annual) || i_annual <= -1) {
    stop("Derived annual effective yield is invalid (<= -1).", call. = FALSE)
  }
  if (!is.finite(j) || j <= -1) {
    stop("Derived yield per period is invalid (<= -1).", call. = FALSE)
  }

  cf_tbl <- bond_cash_flows(
    face = face,
    c = c,
    n = n,
    k = k,
    R = R,
    tol = tol,
    check = FALSE
  )

  # Discount by coupon period.
  period_raw <- cf_tbl$t * k
  if (any(abs(period_raw - round(period_raw)) > tol)) {
    stop("Internal error: coupon times are not aligned with coupon periods.", call. = FALSE)
  }
  period <- as.integer(round(period_raw))

  v <- 1 / (1 + j)
  pv <- cf_tbl$cf * v^period
  price <- sum(pv)

  if (!is.finite(price) || abs(price) <= tol) {
    stop(
      "Bond price is zero or numerically indistinguishable from zero; duration is undefined.",
      call. = FALSE
    )
  }

  # Macaulay duration in periods and years
  macaulay_duration_periods <- sum(period * pv) / price
  macaulay_duration_years <- macaulay_duration_periods / k

  # Modified duration
  modified_duration_periods_j <- macaulay_duration_periods / (1 + j)
  modified_duration_years_i <- macaulay_duration_years / (1 + i_annual)

  tibble::tibble(
    price = price,
    macaulay_duration_periods = macaulay_duration_periods,
    macaulay_duration_years = macaulay_duration_years,
    modified_duration_periods_j = modified_duration_periods_j,
    modified_duration_years_i = modified_duration_years_i,
    yield_per_period = j,
    yield_effective_annual = i_annual,
    k = k,
    n_periods = N
  )
}
