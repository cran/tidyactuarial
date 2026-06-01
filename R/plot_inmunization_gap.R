#' Plot immunization performance under interest-rate shifts
#'
#' Computes and plots the difference between the present value of liabilities
#' and the present value of an immunized asset portfolio under small interest
#' rate changes. This allows visual evaluation of duration or
#' duration-convexity immunization quality.
#'
#' @param L Numeric vector of liability payments.
#' @param t Numeric vector of times of each liability payment.
#' @param asset_cashflows A list where each element is a list with components
#'   \code{$cf} and \code{$t}, defining each asset's cash flow.
#' @param w Numeric vector of portfolio weights or units. Must have the same
#'   length as \code{asset_cashflows}.
#' @param i Base interest-rate input.
#' @param i_type Character string indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal rates. Ignored
#'   for \code{i_type = "effective"} and \code{i_type = "force"}.
#' @param delta A numeric value defining the range of annual effective rates:
#'   from \code{i_effective - delta} to \code{i_effective + delta}, where
#'   \code{i_effective} is the annual effective rate equivalent to
#'   \code{i}, \code{i_type}, and \code{m}.
#' @param n_grid Number of rate values to evaluate.
#'
#' @return A \code{ggplot2} object showing the PV difference curve
#'   \eqn{PV_A(i) - PV_L(i)}{PV_A(i) - PV_L(i)} and a zero reference line.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{L} denotes liabilities, \code{t} denotes payment
#' times, \code{w} denotes asset weights, \code{cf} denotes cash flows,
#' \code{i} denotes the interest-rate input, \code{i_type} denotes the
#' interest-rate type, and \code{m} denotes the conversion frequency for
#' nominal rates.
#'
#' Let \eqn{v(i) = 1/(1+i)}{v(i) = 1/(1+i)}. For a liability stream
#' \eqn{L_k}{L_k} at time \eqn{t_k}{t_k}:
#' \deqn{PV_L(i) = \sum_k L_k \, v(i)^{t_k}}{PV_L(i) = sum(L_k * v(i)^t_k)}
#'
#' For a portfolio of assets with weights \eqn{w_j}{w_j}:
#' \deqn{PV_A(i) = \sum_j w_j \, PV_j(i)}{PV_A(i) = sum(w_j * PV_j(i))}
#'
#' The curve \eqn{\Delta(i) = PV_A(i) - PV_L(i)}{Delta(i) = PV_A(i) - PV_L(i)}
#' illustrates immunization robustness. Under perfect duration immunization,
#' this curve is tangent to zero at the base rate and non-negative nearby if
#' the convexity condition is also met.
#'
#' @seealso \code{\link{immunize_duration}},
#'   \code{\link{immunize_duration_convexity}},
#'   \code{\link{bond_duration}}, \code{\link{bond_convexity}}
#'
#' @family immunization
#'
#' @examples
#' # Two-asset duration immunization gap
#' plot_immunization_gap(
#'   L = c(5000, 8000),
#'   t = c(3, 7),
#'   asset_cashflows = list(
#'     list(cf = c(0, 0, 100), t = c(1, 2, 3)),
#'     list(cf = c(0, 0, 0, 0, 0, 0, 200), t = 1:7)
#'   ),
#'   w = c(5, 2.5),
#'   i = 0.05,
#'   delta = 0.02
#' )
#'
#' @export
plot_immunization_gap <- function(
    L,
    t,
    asset_cashflows,
    w,
    i,
    i_type = "effective",
    m = 1L,
    delta = 0.01,
    n_grid = 200L
) {
  if (missing(L) || !is.numeric(L)) {
    stop("`L` must be a numeric vector of liabilities.", call. = FALSE)
  }

  if (missing(t) || !is.numeric(t)) {
    stop("`t` must be a numeric vector of liability times.", call. = FALSE)
  }

  if (length(L) != length(t)) {
    stop("`L` and `t` must have the same length.", call. = FALSE)
  }

  if (length(L) == 0L) {
    stop("`L` and `t` must have positive length.", call. = FALSE)
  }

  if (any(is.na(L)) || any(!is.finite(L))) {
    stop("`L` must contain only finite numeric values.", call. = FALSE)
  }

  if (any(is.na(t)) || any(!is.finite(t)) || any(t < 0)) {
    stop("`t` must contain only finite values greater than or equal to 0.", call. = FALSE)
  }

  if (!is.list(asset_cashflows) || length(asset_cashflows) == 0L) {
    stop("`asset_cashflows` must be a non-empty list.", call. = FALSE)
  }

  if (missing(w) || !is.numeric(w)) {
    stop("`w` must be a numeric vector of asset weights.", call. = FALSE)
  }

  if (length(asset_cashflows) != length(w)) {
    stop("`asset_cashflows` and `w` must have the same length.", call. = FALSE)
  }

  if (any(is.na(w)) || any(!is.finite(w))) {
    stop("`w` must contain only finite numeric values.", call. = FALSE)
  }

  if (missing(i) || !is.numeric(i) || length(i) != 1L ||
      is.na(i) || !is.finite(i)) {
    stop("`i` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.character(i_type) || length(i_type) != 1L || is.na(i_type)) {
    stop("`i_type` must be a single character string.", call. = FALSE)
  }

  valid_i_type <- c(
    "effective",
    "nominal_interest",
    "nominal_discount",
    "force"
  )

  if (!i_type %in% valid_i_type) {
    stop(
      "`i_type` must be one of: ",
      paste(sprintf("'%s'", valid_i_type), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!is.numeric(m) || length(m) != 1L || is.na(m) ||
      !is.finite(m) || m < 1 || abs(m - round(m)) > 1e-10) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }

  m <- as.integer(round(m))

  if (!is.numeric(delta) || length(delta) != 1L || is.na(delta) ||
      !is.finite(delta) || delta <= 0) {
    stop("`delta` must be a single finite positive number.", call. = FALSE)
  }

  if (!is.numeric(n_grid) || length(n_grid) != 1L || is.na(n_grid) ||
      !is.finite(n_grid) || n_grid < 2 || abs(n_grid - round(n_grid)) > 1e-10) {
    stop("`n_grid` must be a single integer greater than or equal to 2.", call. = FALSE)
  }

  n_grid <- as.integer(round(n_grid))

  i_effective <- standardize_interest(
    type = i_type,
    rate = i,
    m = m
  )

  if (!is.numeric(i_effective) ||
      length(i_effective) != 1L ||
      is.na(i_effective) ||
      !is.finite(i_effective) ||
      i_effective <= -1) {
    stop(
      "The standardized annual effective interest rate must be greater than -1.",
      call. = FALSE
    )
  }

  lower <- max(-0.999999, i_effective - delta)
  upper <- i_effective + delta

  if (upper <= lower) {
    stop("The interest-rate grid is invalid. Reduce `delta`.", call. = FALSE)
  }

  i_grid <- seq(lower, upper, length.out = n_grid)

  pv_liab <- function(rate) {
    v <- 1 / (1 + rate)
    sum(L * v^t)
  }

  pv_asset_j <- function(asset_cf, rate) {
    if (!is.list(asset_cf) ||
        !all(c("cf", "t") %in% names(asset_cf))) {
      stop("Each element of `asset_cashflows` must contain components `cf` and `t`.",
           call. = FALSE)
    }

    cf_j <- asset_cf$cf
    t_j <- asset_cf$t

    if (!is.numeric(cf_j) || !is.numeric(t_j) || length(cf_j) != length(t_j)) {
      stop("Each asset cash flow must have numeric `cf` and `t` of the same length.",
           call. = FALSE)
    }

    if (any(is.na(cf_j)) || any(!is.finite(cf_j)) ||
        any(is.na(t_j)) || any(!is.finite(t_j)) || any(t_j < 0)) {
      stop("Asset cash-flow components `cf` and `t` must be finite, with `t >= 0`.",
           call. = FALSE)
    }

    v <- 1 / (1 + rate)
    sum(cf_j * v^t_j)
  }

  pv_portfolio <- function(rate) {
    sum(vapply(seq_along(asset_cashflows), function(j) {
      w[j] * pv_asset_j(asset_cashflows[[j]], rate)
    }, numeric(1L)))
  }

  PV_L <- vapply(i_grid, pv_liab, numeric(1L))
  PV_A <- vapply(i_grid, pv_portfolio, numeric(1L))
  GAP  <- PV_A - PV_L

  df <- data.frame(i = i_grid, gap = GAP)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$i, y = .data$gap)) +
    ggplot2::geom_hline(yintercept = 0, color = "red", linewidth = 0.8) +
    ggplot2::geom_line(color = "blue", linewidth = 1.1) +
    ggplot2::labs(
      title = "Immunization Gap: PV_assets(i) - PV_liabilities(i)",
      x = "Annual effective interest rate i",
      y = "PV difference"
    ) +
    ggplot2::theme_minimal()
}
