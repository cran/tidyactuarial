#' Cash flow structure of a level coupon bond
#'
#' Builds the cash flow schedule of a level coupon bond with constant coupon
#' rate and a single redemption payment at maturity.
#'
#' @param face Numeric scalar. Face value.
#' @param coupon_rate Numeric scalar. Annual coupon rate.
#' @param years_to_maturity Numeric scalar. Years to maturity.
#' @param coupons_per_year Positive integer. Payments per year.
#' @param redemption Numeric scalar. Redemption value.
#' @param tol Numeric scalar. Tolerance.
#' @param check Logical scalar. Input validation.
#'
#' @return A tibble with the bond schedule.
#'
#' @export
bond_cash_flows <- function(
    face,
    coupon_rate,
    years_to_maturity,
    coupons_per_year = 1L,
    redemption = NULL,
    tol = 1e-10,
    check = TRUE
) {
  if (is.null(redemption)) {
    redemption <- face
  }

  if (isTRUE(check)) {
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
        years_to_maturity < 0) {
      stop("`years_to_maturity` must be a single finite nonnegative number.", call. = FALSE)
    }

    if (!is.numeric(coupons_per_year) || length(coupons_per_year) != 1L ||
        is.na(coupons_per_year) || !is.finite(coupons_per_year) ||
        coupons_per_year <= 0 ||
        abs(coupons_per_year - round(coupons_per_year)) > tol) {
      stop("`coupons_per_year` must be a positive integer.", call. = FALSE)
    }

    if (!is.numeric(redemption) || length(redemption) != 1L || is.na(redemption) ||
        !is.finite(redemption) || redemption < 0) {
      stop("`redemption` must be a single finite nonnegative number.", call. = FALSE)
    }

    if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
      stop("`tol` must be a single positive numeric value.", call. = FALSE)
    }
  }

  k <- as.integer(round(coupons_per_year))

  N_raw <- years_to_maturity * k
  if (abs(N_raw - round(N_raw)) > tol) {
    stop(
      "`years_to_maturity * coupons_per_year` must be an integer ",
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
        time = 0,
        cash_flow = redemption,
        type = "redemption"
      )
    )
  }

  coupon_per_period <- face * coupon_rate / k
  period <- seq_len(N)
  time <- period / k

  coupon_tbl <- tibble::tibble(
    cashflow_id = seq_len(N),
    period = period,
    time = time,
    cash_flow = rep(coupon_per_period, N),
    type = rep("coupon", N)
  )

  redemption_tbl <- tibble::tibble(
    cashflow_id = N + 1L,
    period = N,
    time = N / k,
    cash_flow = redemption,
    type = "redemption"
  )

  dplyr::bind_rows(coupon_tbl, redemption_tbl) |>
    dplyr::arrange(time, cashflow_id)
}
