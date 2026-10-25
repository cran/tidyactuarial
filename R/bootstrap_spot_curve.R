#' Bootstrap a spot yield curve from instrument prices
#'
#' @md
#' @description
#' Recovers annual effective spot rates and discount factors sequentially
#' from observed dirty prices and dated cash flows. Each instrument supplies
#' one new maturity; all its earlier positive cash flows must occur at
#' maturities already bootstrapped.
#'
#' @details
#' For instrument \eqn{j}, with maturity \eqn{t_j}, dirty price \eqn{P_j},
#' and cash flows \eqn{CF_{j,k}}, the price identity is
#' \deqn{P_j = \sum_{k < j} CF_{j,k} v(t_k) + CF_{j,j} v(t_j).}
#' Hence the new discount factor is recovered exactly as
#' \deqn{v(t_j) = \frac{P_j - \sum_{k < j} CF_{j,k} v(t_k)}{CF_{j,j}}.}
#' The annual effective spot rate is
#' \deqn{i_{t_j} = v(t_j)^{-1/t_j} - 1.}
#'
#' Each input row represents one instrument; \code{t} and \code{cf} are
#' list-columns. Repeated payment times within an instrument are combined,
#' allowing coupon and redemption flows produced by [bond_cash_flows()].
#' Instruments are processed by increasing final positive-payment time, not
#' input row order. Zero cash flows are ignored. Duplicate final maturities,
#' missing earlier discount-factor nodes, zero or negative residual prices,
#' and nonpositive terminal cash flows are rejected.
#'
#' This is exact triangular bootstrapping: it does **not** interpolate,
#' extrapolate, fit a curve, or estimate a least-squares solution. An annual
#' effective spot rate can be negative; its discount factor must be positive.
#' Input prices and cash flows must use the same monetary units.
#' Both the positive finite discount factor and its annual effective rate above
#' -1 must be representable. Rounded input prices can make a small terminal
#' residual ill-conditioned; no tolerance silently replaces a nonpositive
#' residual. `repricing_error` reports the recomputed floating-point discrepancy,
#' not an independent accuracy guarantee or a convergence test.
#'
#' @param .data Nonempty data frame with one row per priced instrument.
#' @param col_instrument Character scalar naming the unique instrument-ID column.
#' @param col_price Character scalar naming the observed dirty-price column.
#' @param col_t Character scalar naming the list-column of positive payment
#'   times in years.
#' @param col_cf Character scalar naming the list-column of nonnegative cash
#'   flows corresponding to \code{col_t}.
#'
#' @return A tibble ordered by maturity, with columns \code{instrument},
#'   \code{t}, \code{price}, \code{known_pv}, \code{terminal_cf}, \code{v},
#'   \code{i}, \code{repriced_price}, and \code{repricing_error}. The columns
#'   \code{v} and \code{i} contain the new maturity's discount factor and
#'   annual effective spot rate. \code{repricing_error} is the model price
#'   minus the observed dirty price, in monetary units.
#'
#' @seealso [yield_curve()], [discount_factor_spot()], [forward_rate()],
#'   [bond_cash_flows()]
#' @family interest
#'
#' @examples
#' # Recover a known 1-, 2-, and 3-year annual effective spot curve.
#' t <- 1:3
#' i <- c(0.04, 0.045, 0.05)
#' v <- (1 + i)^(-t)
#' instruments <- tibble::tibble(
#'   instrument = c("Zero 1Y", "Coupon 2Y", "Coupon 3Y"),
#'   price = c(100 * v[1],
#'             5 * v[1] + 105 * v[2],
#'             6 * v[1] + 6 * v[2] + 106 * v[3]),
#'   t = list(1, c(1, 2), c(1, 2, 3)),
#'   cf = list(100, c(5, 105), c(6, 6, 106))
#' )
#' curve <- instruments |> bootstrap_spot_curve()
#' curve
#' yield_curve(t = curve$t, i = curve$i)
#'
#' @export
bootstrap_spot_curve <- function(
    .data,
    col_instrument = "instrument",
    col_price = "price",
    col_t = "t",
    col_cf = "cf"
) {
  if (!is.data.frame(.data) || nrow(.data) == 0L) {
    stop("`.data` must be a nonempty data frame with one row per instrument.",
         call. = FALSE)
  }

  col_names <- list(
    col_instrument = col_instrument,
    col_price = col_price,
    col_t = col_t,
    col_cf = col_cf
  )
  for (arg in names(col_names)) {
    value <- col_names[[arg]]
    if (!is.character(value) || length(value) != 1L ||
        is.na(value) || !nzchar(value)) {
      stop(sprintf("`%s` must be a single nonempty column name.", arg),
           call. = FALSE)
    }
  }
  required <- unlist(col_names, use.names = FALSE)
  if (anyDuplicated(required)) {
    stop("Input column names must be distinct.", call. = FALSE)
  }
  absent <- setdiff(required, names(.data))
  if (length(absent)) {
    stop(paste0("Missing required column(s): ",
                paste(absent, collapse = ", "), "."), call. = FALSE)
  }

  id <- .data[[col_instrument]]
  if (is.factor(id)) id <- as.character(id)
  if (!is.character(id) || anyNA(id) || any(!nzchar(id)) ||
      anyDuplicated(id)) {
    stop("Instrument IDs must be unique, nonmissing, nonempty strings.",
         call. = FALSE)
  }

  price <- .data[[col_price]]
  if (!is.numeric(price) || length(price) != nrow(.data) ||
      anyNA(price) || any(!is.finite(price)) || any(price <= 0)) {
    stop("Prices must be finite positive numeric values.", call. = FALSE)
  }

  t_col <- .data[[col_t]]
  cf_col <- .data[[col_cf]]
  if (!is.list(t_col) || !is.list(cf_col)) {
    stop("Time and cash-flow columns must be list-columns of numeric vectors.",
         call. = FALSE)
  }

  n <- nrow(.data)
  # Validated, positive cash flows for each instrument.
  schedules <- vector("list", n)
  maturities <- numeric(n)
  for (j in seq_len(n)) {
    tt <- t_col[[j]]
    cash <- cf_col[[j]]
    if (!is.numeric(tt) || !is.null(dim(tt)) ||
        !is.numeric(cash) || !is.null(dim(cash)) ||
        !length(tt) || length(tt) != length(cash) ||
        anyNA(tt) || any(!is.finite(tt)) || any(tt <= 0) ||
        anyNA(cash) || any(!is.finite(cash)) || any(cash < 0)) {
      stop(sprintf("Invalid payment times or cash flows for instrument `%s`.",
                   id[[j]]), call. = FALSE)
    }

    # Sum coupon and redemption paid at the same date, and discard zeros.
    positive <- cash > 0
    if (!any(positive)) {
      stop(sprintf("Instrument `%s` must have a positive cash flow.", id[[j]]),
           call. = FALSE)
    }
    tt <- tt[positive]
    cash <- cash[positive]
    order_t <- order(tt)
    tt <- tt[order_t]
    cash <- cash[order_t]
    times <- unique(tt)
    amounts <- vapply(times, function(s) sum(cash[tt == s]), numeric(1L))
    if (any(!is.finite(amounts))) {
      stop(sprintf("Cash-flow totals overflow for instrument `%s`.", id[[j]]),
           call. = FALSE)
    }
    schedules[[j]] <- list(t = times, cf = amounts)
    maturities[[j]] <- tail(times, 1L)
  }

  if (anyDuplicated(maturities)) {
    stop("Each instrument must introduce a distinct final maturity.",
         call. = FALSE)
  }

  ord <- order(maturities)
  result <- vector("list", n)
  calibrated_t <- numeric(0L)
  calibrated_v <- numeric(0L)

  for (step in seq_len(n)) {
    j <- ord[[step]]
    payment_times <- schedules[[j]]$t
    payment_amounts <- schedules[[j]]$cf
    maturity <- tail(payment_times, 1L)
    terminal_cf <- tail(payment_amounts, 1L)

    earlier <- seq_len(length(payment_times) - 1L)
    if (length(earlier)) {
      known_idx <- match(payment_times[earlier], calibrated_t)
      if (anyNA(known_idx)) {
        unknown <- payment_times[earlier][is.na(known_idx)]
        stop(sprintf(
          "Instrument `%s` has an earlier cash flow at t = %s without a bootstrapped discount factor. No interpolation is performed.",
          id[[j]], paste(unknown, collapse = ", ")
        ), call. = FALSE)
      }
      known_pv <- sum(payment_amounts[earlier] * calibrated_v[known_idx])
    } else {
      known_pv <- 0
    }

    # Keep the market price outside tibble's data mask; `price` is also
    # an output column and must not shadow the original price vector.
    market_price <- price[[j]]
    residual <- market_price - known_pv
    if (!is.finite(residual) || residual <= 0) {
      stop(sprintf(
        "Instrument `%s` implies a nonpositive or nonfinite residual present value.",
        id[[j]]
      ), call. = FALSE)
    }
    v <- residual / terminal_cf
    i <- expm1(-log(v) / maturity)
    if (!is.finite(v) || v <= 0 || !is.finite(i) || i <= -1) {
      stop(sprintf("Instrument `%s` implies an invalid spot discount factor or rate.",
                   id[[j]]), call. = FALSE)
    }

    repriced_price <- known_pv + terminal_cf * v
    result[[step]] <- tibble::tibble(
      instrument = id[[j]],
      t = maturity,
      price = market_price,
      known_pv = known_pv,
      terminal_cf = terminal_cf,
      v = v,
      i = i,
      repriced_price = repriced_price,
      repricing_error = repriced_price - market_price
    )
    calibrated_t <- c(calibrated_t, maturity)
    calibrated_v <- c(calibrated_v, v)
  }

  dplyr::bind_rows(result)
}
