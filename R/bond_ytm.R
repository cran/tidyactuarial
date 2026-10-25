#' Yield to maturity of a level coupon bond
#'
#' Computes the yield to maturity (YTM) of a level coupon bond given its
#' observed dirty price at time 0, using compact actuarial notation.
#'
#' The YTM is solved first as the effective yield per coupon period and then
#' reported together with common annual equivalents.
#'
#' Assumptions:
#' \itemize{
#'   \item Coupons are paid in arrears at regular intervals.
#'   \item Price is observed at a coupon date (no accrued interest).
#'   \item \code{n * k} must be an integer.
#'   \item Stub periods are not supported.
#' }
#'
#' @param P Numeric scalar. Observed dirty price of the bond at time 0.
#' @param face Numeric scalar. Face value of the bond.
#' @param c Numeric scalar. Annual coupon rate as a proportion.
#' @param n Numeric scalar. Time to maturity in years. Must be strictly
#'   positive.
#' @param k Positive integer. Number of coupon payments per year.
#' @param R Numeric scalar. Redemption value at maturity. If \code{NULL},
#'   defaults to \code{face}.
#' @param interval Optional finite numeric vector of length 2 giving an increasing
#'   bracket for the effective yield per coupon period, with both bounds above -1.
#' @param tol Finite positive numeric scalar. Yield-coordinate tolerance passed
#'   to \code{\link[stats]{uniroot}}, not a monetary price tolerance.
#' @param maxiter Positive integer. Maximum number of iterations passed to
#'   \code{\link[stats]{uniroot}}.
#' @param check Logical scalar. If \code{TRUE}, performs basic input checks.
#'
#' @return A one-row tibble with columns:
#' \describe{
#'   \item{price}{Input dirty price.}
#'   \item{i_period}{Effective yield per coupon period.}
#'   \item{j_nominal}{Nominal annual yield convertible \code{k} times per year
#'     (= \code{k * i_period}). When \code{k = 2}, this is the bond-equivalent
#'     yield.}
#'   \item{i_effective_annual}{Annual effective yield.}
#'   \item{k}{Coupon frequency.}
#' }
#'
#' @details
#' This function follows the compact bond notation used in
#' \code{tidyactuarial}: \code{P} is the observed dirty price, \code{face} is
#' the face value, \code{c} is the annual coupon rate, \code{n} is the time to
#' maturity, \code{k} is the coupon frequency, and \code{R} is the redemption
#' value.
#'
#' The effective yield per coupon period \eqn{j} is the solution to
#' \deqn{P = \sum_{r=1}^{N} C_r (1+j)^{-r}.}
#'
#' With nonnegative payments and at least one positive future payment, price
#' is strictly decreasing in \eqn{j > -1}, so the yield is unique. The root is
#' found with \code{\link[stats]{uniroot}}, using the log of the discounted
#' cash-flow sum relative to price to avoid monetary overflow or underflow.
#' If no \code{interval} is supplied, the lower bound is the nearest double
#' precision number above -1 and the upper bound starts at 0.10 and is widened
#' until the root is bracketed. A yield outside the representable domain or an
#' exhausted iteration limit causes an error. Tolerance applies to the yield
#' coordinate according to \code{uniroot}; it does not guarantee a given price
#' residual, particularly for yields near -1.
#'
#' From the per-period yield, the annual equivalents are:
#' \deqn{j^{(k)} = k j, \qquad i = (1 + j)^k - 1.}
#'
#' @seealso \code{\link{bond_price}}, \code{\link{bond_cash_flows}},
#'   \code{\link{bond_duration}}, \code{\link{bond_convexity}},
#'   \code{\link{bond_callable_price}}
#'
#' @family bonds
#'
#' @examples
#' bond_ytm(
#'   P = 100,
#'   face = 100,
#'   c = 0.06,
#'   n = 5,
#'   k = 1
#' )
#'
#' bond_ytm(
#'   P = 950,
#'   face = 1000,
#'   c = 0.05,
#'   n = 10,
#'   k = 2
#' )
#'
#' @export
bond_ytm <- function(
    P,
    face,
    c,
    n,
    k = 1L,
    R = NULL,
    interval = NULL,
    tol = 1e-12,
    maxiter = 1000,
    check = TRUE
) {
  if (is.null(R)) {
    R <- face
  }

  if (isTRUE(check)) {
    if (!is.numeric(P) || length(P) != 1L || is.na(P) ||
        !is.finite(P) || P <= 0) {
      stop("`P` must be a single finite positive number.", call. = FALSE)
    }

    if (!is.numeric(face) || length(face) != 1L || is.na(face) ||
        !is.finite(face) || face < 0) {
      stop("`face` must be a single finite nonnegative number.", call. = FALSE)
    }

    if (missing(c) || !is.numeric(c) || length(c) != 1L || is.na(c) ||
        !is.finite(c) || c < 0) {
      stop("`c` must be a single finite nonnegative number.", call. = FALSE)
    }

    if (missing(n) || !is.numeric(n) || length(n) != 1L ||
        is.na(n) || !is.finite(n) || n <= 0) {
      stop("`n` must be a single finite positive number.", call. = FALSE)
    }

    if (!is.numeric(k) || length(k) != 1L ||
        is.na(k) || !is.finite(k) ||
        k <= 0 || k > .Machine$integer.max ||
        abs(k - round(k)) > 1e-10) {
      stop("`k` must be a positive integer.", call. = FALSE)
    }

    if (!is.numeric(R) || length(R) != 1L || is.na(R) ||
        !is.finite(R) || R < 0) {
      stop("`R` must be a single finite nonnegative number.", call. = FALSE)
    }

    if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) ||
        !is.finite(tol) || tol <= 0) {
      stop("`tol` must be a single finite positive numeric value.", call. = FALSE)
    }

    if (!is.numeric(maxiter) || length(maxiter) != 1L || is.na(maxiter) ||
        !is.finite(maxiter) || maxiter <= 0 || maxiter > .Machine$integer.max ||
        maxiter != floor(maxiter)) {
      stop("`maxiter` must be a positive integer.", call. = FALSE)
    }
  }

  k <- as.integer(round(k))
  maxiter <- as.integer(maxiter)

  cf_tbl <- bond_cash_flows(
    face = face,
    c = c,
    n = n,
    k = k,
    R = R,
    check = FALSE
  )

  positive <- cf_tbl$cf > 0 & cf_tbl$period > 0
  if (!any(positive)) {
    stop("The bond must have at least one positive future cash flow.", call. = FALSE)
  }

  # log(sum(C_r / P * exp(-r * log1p(ip)))) has the same unique root.
  # Preserve accurate ordinary ratios; use separate logs only if a ratio loses range.
  ratio <- cf_tbl$cf[positive] / P
  log_ratio <- log(ratio)
  recover <- !is.finite(log_ratio) | ratio < .Machine$double.xmin
  log_ratio[recover] <- log(cf_tbl$cf[positive][recover]) - log(P)
  periods <- cf_tbl$period[positive]
  f <- function(ip) {
    terms <- log_ratio - periods * log1p(ip)
    anchor <- max(terms)
    anchor + log(sum(exp(terms - anchor)))
  }

  # --- Bracket the root automatically if no interval supplied ---
  if (is.null(interval)) {
    lo <- -1 + .Machine$double.eps / 2
    hi <- 0.10

    flo <- f(lo)
    fhi <- f(hi)

    while (fhi > 0 && hi < .Machine$double.xmax) {
      hi <- if (hi > .Machine$double.xmax / 2) .Machine$double.xmax else 2 * hi + 0.05
      fhi <- f(hi)
    }

    if (flo < 0 || fhi > 0) {
      stop("Could not bracket the YTM within the representable yield domain.", call. = FALSE)
    }

    interval <- c(lo, hi)

  } else {
    if (!is.numeric(interval) || length(interval) != 2L || anyNA(interval) ||
        any(!is.finite(interval)) ||
        interval[1] >= interval[2]) {
      stop(
        "`interval` must be a finite numeric vector of length 2 with interval[1] < interval[2].",
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
    maxiter = maxiter,
    check.conv = TRUE
  )

  ip <- root$root

  tibble::tibble(
    price = P,
    i_period = ip,
    j_nominal = k * ip,
    i_effective_annual = expm1(k * log1p(ip)),
    k = k
  )
}
