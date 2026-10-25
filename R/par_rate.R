#' Par coupon rate implied by a spot yield curve
#'
#' @md
#' @description
#' Computes the level annual coupon rate that makes a bond's time-zero dirty
#' price equal to its face value, using the spot discount factors at each
#' coupon date. Unlike [bond_ytm()], it does not solve for a flat yield.
#'
#' @details
#' For a bond with face value \eqn{F > 0}, redemption \eqn{R}, term \eqn{n},
#' and \eqn{k} coupon payments per year, the price equation is
#' \deqn{P = (F c/k) \sum_{j=1}^{nk} v(j/k) + R v(n).}
#' Setting the dirty price \eqn{P = F} and solving for the nominal annual
#' coupon rate gives
#' \deqn{c_{par} = k [F - R v(n)] / [F \sum_{j=1}^{nk} v(j/k)].}
#' When redemption equals face value, this reduces to
#' \deqn{c_{par} = k [1 - v(n)] / \sum_{j=1}^{nk} v(j/k).}
#'
#' Spot inputs use annual effective rates by default, with
#' \eqn{v(t) = (1 + i_t)^{-t}}. Other supported interest conventions are
#' converted using [discount_factor_spot()]. The calculation uses only
#' observed nodes: every coupon date must match one curve maturity; there is
#' no interpolation or extrapolation. Small floating-point discrepancies
#' within `tol` are accepted as identical dates.
#'
#' Input can be (1) numeric vectors `t` and `i`, (2) a data frame with
#' numeric `t` and `i` columns, such as [bootstrap_spot_curve()] output, or
#' (3) a tibble with `t` and `i` list-columns, such as [yield_curve()] output.
#' In the third form each row describes one curve, and the input columns are
#' retained. `n`, `face`, `k`, and `R` may be scalars or have one value per
#' curve. For a single curve, `n = NULL` uses the largest observed maturity.
#' A coupon rate can be negative when the implied financial terms require it;
#' the function does not impose a contractual zero floor.
#'
#' This model assumes regular coupons paid in arrears, no stub periods,
#' no accrued interest, and a deterministic redemption at maturity.
#' Calculation separates the dimensionless rate from monetary amounts and uses
#' the original log spot growth at maturity to retain small coupon rates.
#' Final amounts and the positive coupon discount sum must be representable.
#' `repricing_error` reports the floating-point price discrepancy; extreme
#' cancellation, especially with negative coupons, can still limit precision.
#'
#' @param .data Optional data frame containing a spot curve as numeric columns,
#'   or multiple curves as list-columns. For pipe-friendly usage.
#' @param t Numeric vector of curve maturities (years) when `.data = NULL`.
#' @param i Numeric vector of observed spot rates when `.data = NULL`.
#' @param n Bond maturity in years. Defaults to the maximum maturity of each
#'   curve; may also be a scalar or one value per curve.
#' @param face Strictly positive face value; scalar or one value per curve.
#' @param k Positive integer coupon frequency per year; scalar or one value
#'   per curve. Coupon dates are `seq_len(n * k) / k`.
#' @param R Redemption amount, scalar or one value per curve; defaults to
#'   `face`.
#' @param i_type Spot rate type (`"effective"`, `"nominal_interest"`,
#'   `"nominal_discount"`, or `"force"`). Passed to [discount_factor_spot()].
#'   Must have length one or one value per curve node.
#' @param m Conversion frequency for nominal spot rates, passed to
#'   [discount_factor_spot()]. Must have length one or one per curve node.
#' @param col_t Name of the spot maturity column in `.data`.
#' @param col_i Name of the spot rate column in `.data`.
#' @param tol Positive numerical tolerance for aligning coupon dates and
#'   curve maturities. Does not authorize interpolation.
#'
#' @return One tibble row per curve, containing the supplied curve columns
#'   and `n`, `face`, `k`, `R`, `annuity_factor`, `v_n`, `c_par`,
#'   `coupon_per_period`, `repriced_price`, and `repricing_error`.
#'   The `c_par` value is the nominal annual coupon rate as a proportion;
#'   `coupon_per_period` is its monetary coupon amount. `repricing_error`
#'   is the reconstructed dirty price minus the face value.
#'
#' @seealso [bootstrap_spot_curve()], [yield_curve()], [bond_cash_flows()],
#'   [bond_price()], [discount_factor_spot()]
#' @family interest
#'
#' @examples
#' # An annual-pay bond whose spot curve was obtained from market prices.
#' curve <- tibble::tibble(
#'   t = c(1, 2, 3),
#'   i = c(0.04, 0.045, 0.05)
#' )
#' curve |> par_rate(face = 100, k = 1)
#'
#' # A two-year bond with semiannual coupons: spot nodes every six months.
#' par_rate(
#'   t = c(0.5, 1, 1.5, 2),
#'   i = rep(0.04, 4),
#'   n = 2,
#'   k = 2
#' )
#'
#' # A one-row curve returned by yield_curve() is also accepted.
#' yield_curve(t = 1:3, i = c(0.04, 0.045, 0.05)) |>
#'   par_rate()
#'
#' @export
par_rate <- function(
    .data = NULL,
    t = NULL,
    i = NULL,
    n = NULL,
    face = 100,
    k = 1L,
    R = NULL,
    i_type = "effective",
    m = 1L,
    col_t = "t",
    col_i = "i",
    tol = 1e-10
) {
  abort <- function(message) stop(message, call. = FALSE)

  if (!is.numeric(tol) || length(tol) != 1L ||
      is.na(tol) || !is.finite(tol) || tol <= 0) {
    abort("`tol` must be a finite positive number.")
  }
  if (!is.character(col_t) || length(col_t) != 1L ||
      is.na(col_t) || !nzchar(col_t) ||
      !is.character(col_i) || length(col_i) != 1L ||
      is.na(col_i) || !nzchar(col_i) || identical(col_t, col_i)) {
    abort("`col_t` and `col_i` must be distinct nonempty column names.")
  }
  if (!is.character(i_type) || !length(i_type) || anyNA(i_type)) {
    abort("`i_type` must be a nonmissing character vector.")
  }
  if (!is.numeric(m) || !length(m) || anyNA(m) ||
      any(!is.finite(m)) || any(m < 1) || any(m != floor(m))) {
    abort("`m` must contain positive integer conversion frequencies.")
  }

  if (is.null(.data)) {
    if (!is.numeric(t) || !is.numeric(i)) {
      abort("Supply numeric `t` and `i` vectors or a `.data` curve.")
    }
    data_in <- tibble::tibble(t = list(t), i = list(i))
    t_name <- "t"
    i_name <- "i"
  } else {
    if (!is.data.frame(.data) || !nrow(.data)) {
      abort("`.data` must be a nonempty data frame.")
    }
    if (!is.null(t) || !is.null(i)) {
      abort("Do not supply `t` or `i` separately when `.data` is provided.")
    }
    if (!all(c(col_t, col_i) %in% names(.data))) {
      abort("The spot maturity and rate columns were not found in `.data`.")
    }
    t_raw <- .data[[col_t]]
    i_raw <- .data[[col_i]]
    if (is.numeric(t_raw) && is.numeric(i_raw)) {
      # One node per row: collapse to one curve (as in bootstrap_spot_curve()).
      data_in <- tibble::tibble(t = list(t_raw), i = list(i_raw))
      t_name <- "t"
      i_name <- "i"
    } else if (is.list(t_raw) && is.list(i_raw)) {
      # One curve per row: preserve curve identifiers and other metadata.
      data_in <- tibble::as_tibble(.data)
      t_name <- col_t
      i_name <- col_i
    } else {
      abort("`t` and `i` must both be numeric columns or both be list-columns.")
    }
  }

  n_curves <- nrow(data_in)
  recycle_curve <- function(x, name, default = NULL) {
    if (is.null(x)) x <- default
    if (!is.numeric(x) || !(length(x) %in% c(1L, n_curves)) ||
        anyNA(x) || any(!is.finite(x))) {
      abort(sprintf("`%s` must contain one finite numeric value per curve, or a scalar.", name))
    }
    rep_len(x, n_curves)
  }
  faces <- recycle_curve(face, "face")
  frequencies <- recycle_curve(k, "k")
  redemptions <- recycle_curve(R, "R", faces)
  maturities <- if (is.null(n)) NULL else recycle_curve(n, "n")

  if (any(faces <= 0)) abort("`face` must be strictly positive.")
  if (any(redemptions < 0)) abort("`R` must be nonnegative.")
  if (any(frequencies < 1 | frequencies != floor(frequencies)) ||
      any(frequencies > .Machine$integer.max)) {
    abort("`k` must contain positive integers.")
  }

  extra <- c("n", "face", "k", "R", "annuity_factor", "v_n", "c_par",
             "coupon_per_period", "repriced_price", "repricing_error")
  conflicting <- intersect(extra, names(data_in))
  if (length(conflicting)) {
    abort(paste0("Output column(s) already exist in `.data`: ",
                 paste(conflicting, collapse = ", "), "."))
  }

  out <- vector("list", n_curves)
  for (row in seq_len(n_curves)) {
    tt <- data_in[[t_name]][[row]]
    ii <- data_in[[i_name]][[row]]
    if (!is.numeric(tt) || !is.null(dim(tt)) ||
        !is.numeric(ii) || !is.null(dim(ii)) ||
        !length(tt) || length(tt) != length(ii) ||
        anyNA(tt) || any(!is.finite(tt)) || any(tt <= 0) ||
        anyNA(ii) || any(!is.finite(ii))) {
      abort(sprintf("Curve %d must have finite numeric time/rate vectors of equal nonzero length; times must be positive.", row))
    }
    idx <- order(tt)
    times <- tt[idx]
    rates <- ii[idx]
    if (length(times) > 1L) {
      gaps <- diff(times)
      adjacent_scale <- pmax(1, abs(times[-1L]), abs(times[-length(times)]))
      if (any(gaps <= tol * adjacent_scale)) {
        abort(sprintf("Curve %d has duplicate or indistinguishable maturity nodes.", row))
      }
    }
    maturity <- if (is.null(maturities)) max(times) else maturities[[row]]
    freq <- frequencies[[row]]
    if (maturity <= 0) abort("`n` must be strictly positive.")
    payment_count <- maturity * freq
    if (!is.finite(payment_count) || payment_count < 1 ||
        abs(payment_count - round(payment_count)) > tol ||
        payment_count > 1e6) {
      abort(sprintf("Curve %d requires `n * k` to be a positive integer (at most 1e6).", row))
    }
    coupon_dates <- seq_len(as.integer(round(payment_count))) / freq
    matched <- vapply(coupon_dates, function(coupon_date) {
      valid <- which(abs(times - coupon_date) <=
                       tol * pmax(1, abs(times), abs(coupon_date)))
      if (length(valid) != 1L) {
        abort(sprintf("Curve %d has no unique discount-factor node at coupon time %.12g; no interpolation is performed.",
                      row, coupon_date))
      }
      valid[[1L]]
    }, integer(1L))

    # Align rate conventions with the sorted spot nodes, when nodewise input
    # has been supplied. Scalar conventions apply to every maturity.
    typ <- if (length(i_type) == length(tt)) i_type[idx] else i_type
    frequency <- if (length(m) == length(tt)) m[idx] else m
    if (!(length(typ) %in% c(1L, length(tt))) ||
        !(length(frequency) %in% c(1L, length(tt)))) {
      abort(sprintf("`i_type` and `m` must be scalar or have one value per node in curve %d.", row))
    }
    discounts <- discount_factor_spot(
      t = times, i = rates, i_type = typ, m = frequency
    )
    if (anyNA(discounts) || any(!is.finite(discounts)) || any(discounts <= 0)) {
      abort(sprintf("Curve %d implies invalid or nonfinite discount factors.", row))
    }

    coupon_v <- discounts[matched]
    annuity <- sum(coupon_v)
    v_terminal <- tail(coupon_v, 1L)
    nominal <- faces[[row]]
    redemption <- redemptions[[row]]
    if (!is.finite(annuity) || annuity <= 0) {
      abort(sprintf("Curve %d has an invalid coupon discount sum.", row))
    }
    # Work per unit of face and retain 1 - exp(z) for small rates and spreads.
    growth <- .interest_basis(i_type = typ, i = rates, m = frequency)$growth
    log_v_terminal <- -times[tail(matched, 1L)] * growth[tail(matched, 1L)]
    ratio <- redemption / nominal
    log_ratio <- if (abs(redemption - nominal) <= nominal / 2) {
      log1p((redemption - nominal) / nominal)
    } else if (is.finite(ratio) && ratio >= .Machine$double.xmin) {
      log(ratio)
    } else {
      log(redemption) - log(nominal)
    }
    z <- log_ratio + log_v_terminal
    balance_sign <- -sign(z)
    log_balance <- if (z > 0) z + log(-expm1(-z)) else log(-expm1(z))
    rate <- .cash_flow_value(balance_sign, log(freq) + log_balance - log(annuity))
    coupon <- .cash_flow_value(balance_sign * nominal, log_balance - log(annuity))
    repriced <- .cash_flow_sum(c(coupon, redemption), c(log(annuity), log_v_terminal))
    if (!is.finite(rate) || !is.finite(coupon) || !is.finite(repriced)) {
      abort(sprintf("Curve %d implies a nonfinite par coupon rate or price.", row))
    }
    out[[row]] <- c(
      n = maturity, face = nominal, k = freq, R = redemption,
      annuity_factor = annuity, v_n = v_terminal, c_par = rate,
      coupon_per_period = coupon, repriced_price = repriced,
      repricing_error = repriced - nominal
    )
  }

  result <- tibble::as_tibble(do.call(rbind, out))
  for (name in names(result)) data_in[[name]] <- result[[name]]
  data_in
}
