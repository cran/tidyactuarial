#' Book value of a level coupon bond at a coupon date
#'
#' Computes the book value of a level coupon bond at one or more coupon dates,
#' under a specified yield basis.
#'
#' The book value is interpreted prospectively: at a valuation time that lies
#' on the coupon grid, it equals the present value at that time of all
#' remaining future coupons and the final redemption amount, discounted at the
#' bond's yield basis.
#'
#' This function interprets \code{valuation_time} as a time immediately after any
#' coupon due at that date has been paid. Therefore:
#' \itemize{
#'   \item at \code{valuation_time = 0}, the book value equals the bond price,
#'   \item at \code{valuation_time = years_to_maturity}, the book value is 0.
#' }
#'
#' Assumptions:
#' \itemize{
#'   \item Coupons are paid in arrears at regular intervals.
#'   \item \code{years_to_maturity * coupons_per_year} must be an integer.
#'   \item Each \code{valuation_time * coupons_per_year} must be an integer.
#'   \item Stub periods are not supported.
#'   \item Valuation is performed at coupon dates; no accrued interest is included.
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
#' @param valuation_time Numeric vector. Valuation time(s) in years, measured
#'   from issue. Each value must lie between \code{0} and \code{years_to_maturity}
#'   and must align with coupon dates.
#' @param coupons_per_year Positive integer. Number of coupon payments per year.
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
#' @param output Character string. Use \code{"value"} to return a numeric vector
#'   of book values, or \code{"table"} to return a tibble with intermediate
#'   quantities.
#'
#' @return
#' If \code{output = "value"}, a numeric vector of book values, one for each
#' \code{valuation_time}.
#'
#' If \code{output = "table"}, a tibble with valuation times, valuation periods,
#' book values, yield information, and bond inputs.
#'
#' @details
#' Let the valuation time correspond to coupon period \eqn{k}, with total
#' maturity period count \eqn{N}. If the remaining future cash flows are
#' \eqn{C_{k+1}, \dots, C_N}, and \eqn{i_p} is the effective yield per coupon
#' period, then the book value at time \eqn{k} is
#' \deqn{BV_k = \sum_{j=k+1}^{N} C_j (1+i_p)^{-(j-k)}}{BV_k = sum_{j=k+1}^{N} C_j (1+i_p)^(-(j-k))}
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
#'   coupon_rate = 0.08,
#'   years_to_maturity = 5,
#'   valuation_time = 0,
#'   coupons_per_year = 2,
#'   yield_rate = 0.06,
#'   yield_rate_type = "effective"
#' )
#'
#' # Book value at several coupon dates
#' bond_book_value(
#'   face = 100,
#'   coupon_rate = 0.08,
#'   years_to_maturity = 5,
#'   valuation_time = c(0, 1, 2, 3, 4, 5),
#'   coupons_per_year = 1,
#'   yield_rate = 0.06,
#'   yield_rate_type = "effective"
#' )
#'
#' # Table output
#' bond_book_value(
#'   face = 100,
#'   coupon_rate = 0.08,
#'   years_to_maturity = 5,
#'   valuation_time = c(0, 1, 2, 3, 4, 5),
#'   coupons_per_year = 1,
#'   yield_rate = 0.06,
#'   yield_rate_type = "effective",
#'   output = "table"
#' )
#'
#' # Yield given directly per coupon period
#' bond_book_value(
#'   face = 1000,
#'   coupon_rate = 0.05,
#'   years_to_maturity = 10,
#'   valuation_time = c(0, 2, 4, 6),
#'   coupons_per_year = 2,
#'   yield_effective_per_period = 0.03
#' )
#'
#' @export
bond_book_value <- function(
    face,
    coupon_rate,
    years_to_maturity,
    valuation_time,
    coupons_per_year = 1L,
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

  out <- .bond_book_value_schedule(
    face = face,
    coupon_rate = coupon_rate,
    years_to_maturity = years_to_maturity,
    valuation_time = valuation_time,
    coupons_per_year = coupons_per_year,
    y_effective_per_period = yield_effective_per_period,
    y_rate = yield_rate,
    y_type = yield_rate_type,
    y_m = yield_m,
    redemption = redemption,
    tol = tol,
    check = check
  )

  if (output == "value") {
    return(out$book_value)
  }

  out
}


# Internal helper: builds book-value table at requested valuation times
.bond_book_value_schedule <- function(
    face,
    coupon_rate,
    years_to_maturity,
    valuation_time,
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
    .validate_bond_core(
      face = face,
      coupon_rate = coupon_rate,
      years_to_maturity = years_to_maturity,
      coupons_per_year = coupons_per_year,
      y_m = y_m,
      redemption = redemption,
      tol = tol
    )

    if (missing(valuation_time)) {
      stop("`valuation_time` must be provided.", call. = FALSE)
    }

    if (!is.numeric(valuation_time) ||
        any(is.na(valuation_time)) ||
        any(!is.finite(valuation_time))) {
      stop("`valuation_time` must be a finite numeric vector.", call. = FALSE)
    }

    if (any(valuation_time < 0) || any(valuation_time > years_to_maturity)) {
      stop(
        "`valuation_time` must lie between 0 and `years_to_maturity`.",
        call. = FALSE
      )
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

  N <- as.integer(round(N_raw))

  val_N_raw <- valuation_time * coupons_per_year

  if (any(abs(val_N_raw - round(val_N_raw)) > tol)) {
    stop(
      "Each `valuation_time * coupons_per_year` must be an integer ",
      "(valuation dates must align with coupon dates).",
      call. = FALSE
    )
  }

  val_N <- as.integer(round(val_N_raw))

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

  cf_tbl <- bond_cash_flows(
    face = face,
    coupon_rate = coupon_rate,
    years_to_maturity = years_to_maturity,
    coupons_per_year = coupons_per_year,
    redemption = redemption,
    tol = tol,
    check = FALSE
  )

  cf_period <- as.integer(round(cf_tbl$time * coupons_per_year))

  book_value <- vapply(
    seq_along(valuation_time),
    function(j) {
      k <- val_N[[j]]

      if (k >= N) {
        return(0)
      }

      idx <- cf_period > k

      if (!any(idx)) {
        return(0)
      }

      rel_periods <- cf_period[idx] - k

      sum(cf_tbl$cash_flow[idx] / (1 + ip)^rel_periods)
    },
    numeric(1L)
  )

  tibble::tibble(
    valuation_time = as.numeric(valuation_time),
    valuation_period = val_N,
    book_value = book_value,
    yield_per_period = ip,
    yield_effective_annual = i_annual,
    coupons_per_year = coupons_per_year,
    face = face,
    coupon_rate = coupon_rate,
    years_to_maturity = years_to_maturity,
    redemption = redemption
  )
}


# ============================================================
# Shared internal utilities for bond functions
# ============================================================

#' Validate common bond scalar inputs
#' @noRd
.validate_bond_core <- function(
    face,
    coupon_rate,
    years_to_maturity,
    coupons_per_year,
    y_m,
    redemption,
    tol
) {
  if (!is.numeric(face) ||
      length(face) != 1L ||
      is.na(face) ||
      !is.finite(face) ||
      face < 0) {
    stop("`face` must be a single finite nonnegative number.", call. = FALSE)
  }

  if (!is.numeric(coupon_rate) ||
      length(coupon_rate) != 1L ||
      is.na(coupon_rate) ||
      !is.finite(coupon_rate) ||
      coupon_rate < 0) {
    stop("`coupon_rate` must be a single finite nonnegative number.", call. = FALSE)
  }

  if (!is.numeric(years_to_maturity) ||
      length(years_to_maturity) != 1L ||
      is.na(years_to_maturity) ||
      !is.finite(years_to_maturity) ||
      years_to_maturity < 0) {
    stop("`years_to_maturity` must be a single finite nonnegative number.", call. = FALSE)
  }

  if (!is.numeric(coupons_per_year) ||
      length(coupons_per_year) != 1L ||
      is.na(coupons_per_year) ||
      !is.finite(coupons_per_year) ||
      coupons_per_year <= 0 ||
      abs(coupons_per_year - round(coupons_per_year)) > tol) {
    stop("`coupons_per_year` must be a positive integer.", call. = FALSE)
  }

  if (!is.numeric(y_m) ||
      length(y_m) != 1L ||
      is.na(y_m) ||
      !is.finite(y_m) ||
      y_m <= 0 ||
      abs(y_m - round(y_m)) > tol) {
    stop("`y_m` must be a positive integer.", call. = FALSE)
  }

  if (!is.numeric(redemption) ||
      length(redemption) != 1L ||
      is.na(redemption) ||
      !is.finite(redemption) ||
      redemption < 0) {
    stop("`redemption` must be a single finite nonnegative number.", call. = FALSE)
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
#' @return A list with \code{ip} (effective yield per period) and \code{i_annual}.
#' @noRd
.resolve_bond_yield <- function(
    y_effective_per_period = NULL,
    y_rate = NULL,
    y_type = "effective",
    y_m = 1L,
    coupons_per_year = 1L
) {
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
    i_annual <- (1 + ip)^coupons_per_year - 1
  } else {
    if (is.null(y_rate)) {
      stop(
        "Provide either `yield_effective_per_period` or `yield_rate`.",
        call. = FALSE
      )
    }

    i_annual <- standardize_interest(
      type = y_type,
      rate = y_rate,
      m = y_m
    )

    ip <- (1 + i_annual)^(1 / coupons_per_year) - 1
  }

  list(ip = ip, i_annual = i_annual)
}
