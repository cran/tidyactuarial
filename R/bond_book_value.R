#' Book value of a level coupon bond at a coupon date
#'
#' Computes the book value of a level coupon bond at one or more coupon dates,
#' under a specified yield basis, using compact actuarial notation.
#'
#' The book value is interpreted prospectively: at a valuation time that lies
#' on the coupon grid, it equals the present value at that time of all
#' remaining future coupons and the final redemption amount, discounted at the
#' bond's yield basis.
#'
#' This function interprets \code{t} as a time immediately after any coupon due
#' at that date has been paid. Therefore:
#' \itemize{
#'   \item at \code{t = 0}, the book value equals the bond price,
#'   \item at \code{t = n}, the book value is 0.
#' }
#'
#' Assumptions:
#' \itemize{
#'   \item Coupons are paid in arrears at regular intervals.
#'   \item \code{n * k} must be an integer.
#'   \item Each \code{t * k} must be an integer.
#'   \item Stub periods are not supported.
#'   \item Valuation is performed at coupon dates; no accrued interest is included.
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
#' @param n Numeric scalar. Final maturity in years.
#' @param t Numeric vector. Valuation time(s) in years, measured from issue.
#'   Each value must lie between \code{0} and \code{n} and must align with
#'   coupon dates.
#' @param k Positive integer. Number of coupon payments per year.
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
#' @param tidy Logical scalar. If \code{FALSE}, returns a numeric vector of book
#'   values. If \code{TRUE}, returns a tibble with intermediate quantities.
#'
#' @return
#' If \code{tidy = FALSE}, a numeric vector of book values, one for each
#' \code{t}.
#'
#' If \code{tidy = TRUE}, a tibble with valuation times, valuation periods, book
#' values, yield information, and bond inputs.
#'
#' @details
#' This function follows the compact bond notation used in
#' \code{tidyactuarial}: \code{P} is price, \code{face} is the face value,
#' \code{c} is the annual coupon rate, \code{n} is maturity, \code{k} is coupon
#' frequency, \code{y} is the yield, \code{R} is redemption value, and \code{t}
#' is valuation time.
#'
#' Let the valuation time correspond to coupon period \eqn{s}, with total
#' maturity period count \eqn{N}. If the remaining future cash flows are
#' \eqn{C_{s+1}, \dots, C_N}, and \eqn{i_p} is the effective yield per coupon
#' period, then the book value at time \eqn{s} is
#' \deqn{BV_s = \sum_{r=s+1}^{N} C_r (1+i_p)^{-(r-s)}.}
#'
#' This is the prospective book value on the bond's yield basis.
#'
#' @seealso \code{\link{bond_price}}, \code{\link{bond_cash_flows}},
#'   \code{\link{bond_duration}}, \code{\link{bond_convexity}}
#'
#' @family bonds
#'
#' @examples
#' # Book value at time 0 equals price
#' bond_book_value(
#'   face = 100,
#'   c = 0.08,
#'   n = 5,
#'   t = 0,
#'   k = 2,
#'   y = 0.06,
#'   y_type = "effective"
#' )
#'
#' # Book value at several coupon dates
#' bond_book_value(
#'   face = 100,
#'   c = 0.08,
#'   n = 5,
#'   t = c(0, 1, 2, 3, 4, 5),
#'   k = 1,
#'   y = 0.06,
#'   y_type = "effective"
#' )
#'
#' # Tidy output
#' bond_book_value(
#'   face = 100,
#'   c = 0.08,
#'   n = 5,
#'   t = c(0, 1, 2, 3, 4, 5),
#'   k = 1,
#'   y = 0.06,
#'   y_type = "effective",
#'   tidy = TRUE
#' )
#'
#' # Yield given directly per coupon period
#' bond_book_value(
#'   face = 1000,
#'   c = 0.05,
#'   n = 10,
#'   t = c(0, 2, 4, 6),
#'   k = 2,
#'   y_effective_per_period = 0.03
#' )
#'
#' @export
bond_book_value <- function(
    face,
    c,
    n,
    t,
    k = 1L,
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

  out <- .bond_book_value_schedule(
    face = face,
    c = c,
    n = n,
    t = t,
    k = k,
    y_effective_per_period = y_effective_per_period,
    y = y,
    y_type = y_type,
    y_m = y_m,
    R = R,
    tol = tol,
    check = check
  )

  if (!tidy) {
    return(out$book_value)
  }

  out
}


# Internal helper: builds book-value table at requested valuation times
.bond_book_value_schedule <- function(
    face,
    c,
    n,
    t,
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
    .validate_bond_core(
      face = face,
      c = c,
      n = n,
      k = k,
      y_m = y_m,
      R = R,
      tol = tol
    )

    if (missing(t)) {
      stop("`t` must be provided.", call. = FALSE)
    }

    if (!is.numeric(t) ||
        any(is.na(t)) ||
        any(!is.finite(t))) {
      stop("`t` must be a finite numeric vector.", call. = FALSE)
    }

    if (any(t < 0) || any(t > n)) {
      stop(
        "`t` must lie between 0 and `n`.",
        call. = FALSE
      )
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

  N <- as.integer(round(N_raw))

  val_N_raw <- t * k

  if (any(abs(val_N_raw - round(val_N_raw)) > tol)) {
    stop(
      "Each `t * k` must be an integer ",
      "(valuation dates must align with coupon dates).",
      call. = FALSE
    )
  }

  val_N <- as.integer(round(val_N_raw))

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

  cf_tbl <- bond_cash_flows(
    face = face,
    c = c,
    n = n,
    k = k,
    R = R,
    tol = tol,
    check = FALSE
  )

  cf_period <- as.integer(round(cf_tbl$t * k))

  book_value <- vapply(
    seq_along(t),
    function(j) {
      s <- val_N[[j]]

      if (s >= N) {
        return(0)
      }

      idx <- cf_period > s

      if (!any(idx)) {
        return(0)
      }

      rel_periods <- cf_period[idx] - s

      sum(cf_tbl$cf[idx] / (1 + ip)^rel_periods)
    },
    numeric(1L)
  )

  tibble::tibble(
    t = as.numeric(t),
    valuation_period = val_N,
    book_value = book_value,
    yield_per_period = ip,
    yield_effective_annual = i_annual,
    k = k,
    face = face,
    c = c,
    n = n,
    R = R
  )
}


# ============================================================
# Shared internal utilities for bond functions
# ============================================================

#' Validate common bond scalar inputs
#' @noRd
.validate_bond_core <- function(
    face,
    c = NULL,
    n = NULL,
    k = 1L,
    y_m = 1L,
    R = NULL,
    tol = 1e-10,
    coupon_rate = NULL,
    years_to_maturity = NULL,
    coupons_per_year = NULL,
    redemption = NULL
) {
  # Transitional internal compatibility for bond functions not yet migrated.
  if (is.null(c) && !is.null(coupon_rate)) {
    c <- coupon_rate
  }

  if (is.null(n) && !is.null(years_to_maturity)) {
    n <- years_to_maturity
  }

  if (!missing(k) && is.null(coupons_per_year)) {
    coupons_per_year <- k
  } else if (!is.null(coupons_per_year)) {
    k <- coupons_per_year
  }

  if (is.null(R) && !is.null(redemption)) {
    R <- redemption
  }

  if (is.null(R)) {
    R <- face
  }

  if (!is.numeric(face) ||
      length(face) != 1L ||
      is.na(face) ||
      !is.finite(face) ||
      face < 0) {
    stop("`face` must be a single finite nonnegative number.", call. = FALSE)
  }

  if (!is.numeric(c) ||
      length(c) != 1L ||
      is.na(c) ||
      !is.finite(c) ||
      c < 0) {
    stop("`c` must be a single finite nonnegative number.", call. = FALSE)
  }

  if (!is.numeric(n) ||
      length(n) != 1L ||
      is.na(n) ||
      !is.finite(n) ||
      n < 0) {
    stop("`n` must be a single finite nonnegative number.", call. = FALSE)
  }

  if (!is.numeric(k) ||
      length(k) != 1L ||
      is.na(k) ||
      !is.finite(k) ||
      k <= 0 ||
      abs(k - round(k)) > tol) {
    stop("`k` must be a positive integer.", call. = FALSE)
  }

  if (!is.numeric(y_m) ||
      length(y_m) != 1L ||
      is.na(y_m) ||
      !is.finite(y_m) ||
      y_m <= 0 ||
      abs(y_m - round(y_m)) > tol) {
    stop("`y_m` must be a positive integer.", call. = FALSE)
  }

  if (!is.numeric(R) ||
      length(R) != 1L ||
      is.na(R) ||
      !is.finite(R) ||
      R < 0) {
    stop("`R` must be a single finite nonnegative number.", call. = FALSE)
  }

  if (!is.numeric(tol) ||
      length(tol) != 1L ||
      is.na(tol) ||
      tol <= 0) {
    stop("`tol` must be a single positive numeric value.", call. = FALSE)
  }

  invisible(TRUE)
}


#' Resolve bond yield from either per-period or annual specification
#' @return A list with \code{ip} (effective yield per period) and
#'   \code{i_annual}.
#' @noRd
.resolve_bond_yield <- function(
    y_effective_per_period = NULL,
    y = NULL,
    y_type = "effective",
    y_m = 1L,
    k = 1L,
    y_rate = NULL,
    coupons_per_year = NULL
) {
  # Transitional internal compatibility for bond functions not yet migrated.
  if (is.null(y) && !is.null(y_rate)) {
    y <- y_rate
  }

  if (!is.null(coupons_per_year)) {
    k <- coupons_per_year
  }

  if (!is.null(y_effective_per_period)) {
    if (!is.numeric(y_effective_per_period) ||
        length(y_effective_per_period) != 1L ||
        is.na(y_effective_per_period) ||
        !is.finite(y_effective_per_period) ||
        y_effective_per_period <= -1) {
      stop(
        "`y_effective_per_period` must be a single finite numeric value greater than -1.",
        call. = FALSE
      )
    }

    ip <- y_effective_per_period
    i_annual <- (1 + ip)^k - 1
  } else {
    if (is.null(y)) {
      stop(
        "Provide either `y_effective_per_period` or `y`.",
        call. = FALSE
      )
    }

    i_annual <- standardize_interest(
      type = y_type,
      rate = y,
      m = y_m
    )

    ip <- (1 + i_annual)^(1 / k) - 1
  }

  list(ip = ip, i_annual = i_annual)
}
