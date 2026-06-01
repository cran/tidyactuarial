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
#' The root is found numerically using \code{\link[stats]{uniroot}}. If no
#' \code{interval} is supplied, the function automatically brackets the root
#' starting from \eqn{(-0.999999, 0.10)} and progressively widens the upper
#' bound until a sign change is detected.
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
        k <= 0 ||
        abs(k - round(k)) > 1e-10) {
      stop("`k` must be a positive integer.", call. = FALSE)
    }

    if (!is.numeric(R) || length(R) != 1L || is.na(R) ||
        !is.finite(R) || R < 0) {
      stop("`R` must be a single finite nonnegative number.", call. = FALSE)
    }

    if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
      stop("`tol` must be a single positive numeric value.", call. = FALSE)
    }

    if (!is.numeric(maxiter) || length(maxiter) != 1L || is.na(maxiter) ||
        !is.finite(maxiter) || maxiter <= 0 || maxiter != floor(maxiter)) {
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

  if (nrow(cf_tbl) == 0L || sum(cf_tbl$cf > 0) == 0L) {
    stop("The bond must have at least one positive future cash flow.", call. = FALSE)
  }

  # Objective function: bond_price(ip) - observed price = 0
  f <- function(ip) {
    bond_price(
      face = face,
      c = c,
      n = n,
      k = k,
      y_effective_per_period = ip,
      R = R,
      check = FALSE
    ) - P
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
    price = P,
    i_period = ip,
    j_nominal = k * ip,
    i_effective_annual = (1 + ip)^k - 1,
    k = k
  )
}
