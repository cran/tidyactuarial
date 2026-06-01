#' Cash flow structure of a level coupon bond
#'
#' Builds the cash-flow schedule of a level coupon bond with constant coupon
#' rate and a single redemption payment at maturity, using compact actuarial
#' notation.
#'
#' @param face Numeric scalar. Face value of the bond.
#' @param c Numeric scalar. Annual coupon rate.
#' @param n Numeric scalar. Years to maturity.
#' @param k Positive integer. Number of coupon payments per year.
#' @param R Numeric scalar. Redemption value. If \code{NULL}, it is set equal
#'   to \code{face}.
#' @param tol Numeric scalar. Tolerance.
#' @param check Logical scalar. Input validation.
#'
#' @return A tibble with the bond cash-flow schedule. The main actuarial
#'   columns are \code{t} for payment time and \code{cf} for cash flow.
#'
#' @details
#' This function follows the compact bond notation used in
#' \code{tidyactuarial}: \code{face} is the face value, \code{c} is the annual
#' coupon rate, \code{n} is the term to maturity, \code{k} is the number of
#' coupon payments per year, and \code{R} is the redemption value.
#'
#' Stub periods are not supported; therefore, \code{n * k} must be an integer.
#'
#' @examples
#' bond_cash_flows(
#'   face = 1000,
#'   c = 0.05,
#'   n = 10,
#'   k = 2,
#'   R = 1000
#' )
#'
#' @export
bond_cash_flows <- function(
    face,
    c,
    n,
    k = 1L,
    R = NULL,
    tol = 1e-10,
    check = TRUE
) {
  if (is.null(R)) {
    R <- face
  }

  if (isTRUE(check)) {
    if (!is.numeric(face) || length(face) != 1L || is.na(face) ||
        !is.finite(face) || face < 0) {
      stop("`face` must be a single finite nonnegative number.", call. = FALSE)
    }

    if (missing(c) || !is.numeric(c) || length(c) != 1L || is.na(c) ||
        !is.finite(c) || c < 0) {
      stop("`c` must be a single finite nonnegative number.", call. = FALSE)
    }

    if (missing(n) || !is.numeric(n) || length(n) != 1L ||
        is.na(n) || !is.finite(n) || n < 0) {
      stop("`n` must be a single finite nonnegative number.", call. = FALSE)
    }

    if (!is.numeric(k) || length(k) != 1L ||
        is.na(k) || !is.finite(k) ||
        k <= 0 ||
        abs(k - round(k)) > tol) {
      stop("`k` must be a positive integer.", call. = FALSE)
    }

    if (!is.numeric(R) || length(R) != 1L || is.na(R) ||
        !is.finite(R) || R < 0) {
      stop("`R` must be a single finite nonnegative number.", call. = FALSE)
    }

    if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
      stop("`tol` must be a single positive numeric value.", call. = FALSE)
    }
  }

  k <- as.integer(round(k))

  N_raw <- n * k
  if (abs(N_raw - round(N_raw)) > tol) {
    stop(
      "`n * k` must be an integer ",
      "(stub periods are not supported).",
      call. = FALSE
    )
  }
  N <- as.integer(round(N_raw))

  # Maturity at time 0: immediate redemption only
  if (N == 0L) {
    return(
      tibble::tibble(
        cashflow_id = 1L,
        period = 0L,
        t = 0,
        cf = R,
        type = "redemption"
      )
    )
  }

  coupon_per_period <- face * c / k
  period <- seq_len(N)
  t <- period / k

  coupon_tbl <- tibble::tibble(
    cashflow_id = seq_len(N),
    period = period,
    t = t,
    cf = rep(coupon_per_period, N),
    type = rep("coupon", N)
  )

  redemption_tbl <- tibble::tibble(
    cashflow_id = N + 1L,
    period = N,
    t = N / k,
    cf = R,
    type = "redemption"
  )

  dplyr::bind_rows(coupon_tbl, redemption_tbl) |>
    dplyr::arrange(t, cashflow_id)
}
