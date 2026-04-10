#' Yield to maturity of a level coupon bond
#'
#' Computes the yield to maturity (YTM) of a level coupon bond given its
#' observed dirty price at time 0.
#'
#' The YTM is solved first as the effective yield per coupon period and then
#' reported together with common annual equivalents.
#'
#' Assumptions:
#' \itemize{
#'   \item Coupons are paid in arrears at regular intervals.
#'   \item Price is observed at a coupon date (no accrued interest).
#'   \item \code{years_to_maturity * coupons_per_year} must be an integer.
#'   \item Stub periods are not supported.
#' }
#'
#' @param price Numeric scalar. Observed dirty price of the bond at time 0.
#' @param face Numeric scalar. Face (par) value of the bond.
#' @param coupon_rate Numeric scalar. Annual coupon rate as a proportion.
#' @param years_to_maturity Numeric scalar. Time to maturity in years.
#'   Must be strictly positive.
#' @param coupons_per_year Positive integer. Number of coupon payments per year.
#' @param redemption Numeric scalar. Redemption value at maturity.
#'   If \code{NULL}, defaults to \code{face}.
#' @param interval Optional numeric vector of length 2 giving a bracket for the
#'   effective yield per coupon period.
#' @param tol Numeric scalar. Tolerance passed to \code{\link[stats]{uniroot}}.
#' @param maxiter Positive integer. Maximum number of iterations passed to
#'   \code{\link[stats]{uniroot}}.
#' @param check Logical scalar. If \code{TRUE}, performs basic input checks.
#'
#' @return A one-row tibble with columns:
#' \describe{
#'   \item{price}{Input dirty price.}
#'   \item{i_period}{Effective yield per coupon period.}
#'   \item{j_nominal}{Nominal annual yield convertible \code{coupons_per_year}
#'     times per year (= \code{coupons_per_year * i_period}). When
#'     \code{coupons_per_year = 2}, this is the bond-equivalent yield.}
#'   \item{i_effective_annual}{Annual effective yield.}
#' }
#'
#' @details
#' The effective yield per coupon period \eqn{j} is the solution to
#' \deqn{P = \sum_{k=1}^{N} C_k \, (1+j)^{-k}}{P = sum_(k=1)^(N) C_k (1+j)^-k}
#' where \eqn{P} is the observed price and \eqn{C_k} are the bond's cash
#' flows (coupons and redemption) at coupon periods \eqn{k = 1, \dots, N}.
#'
#' The root is found numerically using \code{\link[stats]{uniroot}}. If no \code{interval}
#' is supplied, the function automatically brackets the root starting from
#' \eqn{(-0.999999, \, 0.10)} and progressively widens the upper bound until
#' a sign change is detected.
#'
#' From the per-period yield, the annual equivalents are:
#' \deqn{j^{(m)} = m \cdot j, \qquad i = (1 + j)^m - 1.}{j^(m) = m * j, i = (1 + j)^m - 1.}
#'
#' @seealso \code{\link{bond_price}}, \code{\link{bond_cash_flows}}, \code{\link{bond_duration}},
#'   \code{\link{bond_convexity}}, \code{\link{bond_callable_price}}
#'
#' @family bonds
#'
#' @examples
#' bond_ytm(
#'   price = 100,
#'   face = 100,
#'   coupon_rate = 0.06,
#'   years_to_maturity = 5,
#'   coupons_per_year = 1
#' )
#'
#' bond_ytm(
#'   price = 950,
#'   face = 1000,
#'   coupon_rate = 0.05,
#'   years_to_maturity = 10,
#'   coupons_per_year = 2
#' )
#'
#' @export
bond_ytm <- function(
    price,
    face,
    coupon_rate,
    years_to_maturity,
    coupons_per_year = 1L,
    redemption = NULL,
    interval = NULL,
    tol = 1e-12,
    maxiter = 1000,
    check = TRUE
) {
  if (is.null(redemption)) {
    redemption <- face
  }

  if (isTRUE(check)) {
    if (!is.numeric(price) || length(price) != 1L || is.na(price) ||
        !is.finite(price) || price <= 0) {
      stop("`price` must be a single finite positive number.", call. = FALSE)
    }

    if (!is.numeric(face) || length(face) != 1L || is.na(face) ||
        !is.finite(face) || face < 0) {
      stop("`face` must be a single finite nonnegative number.", call. = FALSE)
    }

    if (!is.numeric(coupon_rate) || length(coupon_rate) != 1L || is.na(coupon_rate) ||
        !is.finite(coupon_rate) || coupon_rate < 0) {
      stop("`coupon_rate` must be a single finite nonnegative number.", call. = FALSE)
    }

    if (!is.numeric(years_to_maturity) || length(years_to_maturity) != 1L ||
        is.na(years_to_maturity) || !is.finite(years_to_maturity) ||
        years_to_maturity <= 0) {
      stop("`years_to_maturity` must be a single finite positive number.", call. = FALSE)
    }

    if (!is.numeric(coupons_per_year) || length(coupons_per_year) != 1L ||
        is.na(coupons_per_year) || !is.finite(coupons_per_year) ||
        coupons_per_year <= 0 ||
        abs(coupons_per_year - round(coupons_per_year)) > 1e-10) {
      stop("`coupons_per_year` must be a positive integer.", call. = FALSE)
    }

    if (!is.numeric(redemption) || length(redemption) != 1L || is.na(redemption) ||
        !is.finite(redemption) || redemption < 0) {
      stop("`redemption` must be a single finite nonnegative number.", call. = FALSE)
    }

    if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
      stop("`tol` must be a single positive numeric value.", call. = FALSE)
    }

    if (!is.numeric(maxiter) || length(maxiter) != 1L || is.na(maxiter) ||
        !is.finite(maxiter) || maxiter <= 0 || maxiter != floor(maxiter)) {
      stop("`maxiter` must be a positive integer.", call. = FALSE)
    }
  }

  m <- as.integer(round(coupons_per_year))
  maxiter <- as.integer(maxiter)

  cf_tbl <- bond_cash_flows(
    face = face,
    coupon_rate = coupon_rate,
    years_to_maturity = years_to_maturity,
    coupons_per_year = m,
    redemption = redemption,
    check = FALSE
  )

  if (nrow(cf_tbl) == 0L || sum(cf_tbl$cash_flow > 0) == 0L) {
    stop("The bond must have at least one positive future cash flow.", call. = FALSE)
  }

  # Objective function: bond_price(ip) - observed price = 0
  f <- function(ip) {
    bond_price(
      face = face,
      coupon_rate = coupon_rate,
      years_to_maturity = years_to_maturity,
      coupons_per_year = m,
      y_effective_per_period = ip,
      redemption = redemption,
      check = FALSE
    ) - price
  }

  # --- Bracket the root automatically if no interval supplied ---
  if (is.null(interval)) {
    lo <- -0.999999
    hi <- 0.10

    flo <- f(lo)
    fhi <- f(hi)

    it <- 0L
    while (is.finite(flo) && is.finite(fhi) && flo * fhi > 0 && it < 100L) {
      hi <- 2 * hi + 0.05
      fhi <- f(hi)
      it <- it + 1L
    }

    if (!(is.finite(flo) && is.finite(fhi) && flo * fhi <= 0)) {
      stop("Could not bracket the YTM. Provide a wider `interval`.", call. = FALSE)
    }

    interval <- c(lo, hi)

  } else {
    if (!is.numeric(interval) || length(interval) != 2L || anyNA(interval) ||
        interval[1] >= interval[2]) {
      stop(
        "`interval` must be a numeric vector of length 2 with interval[1] < interval[2].",
        call. = FALSE
      )
    }
    if (interval[1] <= -1) {
      stop("The lower bound of `interval` must be greater than -1.", call. = FALSE)
    }
  }

  root <- stats::uniroot(
    f,
    interval = interval,
    tol = tol,
    maxiter = maxiter
  )

  ip <- root$root

  tibble::tibble(
    price = price,
    i_period = ip,
    j_nominal = m * ip,
    i_effective_annual = (1 + ip)^m - 1
  )
}
