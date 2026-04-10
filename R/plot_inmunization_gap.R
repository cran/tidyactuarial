#' Plot immunization performance under interest rate shifts
#'
#' Computes and plots the difference between the present value of liabilities
#' and the present value of an immunized asset portfolio under small interest
#' rate changes. This allows visual evaluation of duration or duration-convexity
#' immunization quality.
#'
#' @param liabilities Numeric vector of liability payments.
#' @param t_liabilities Numeric vector of times (periods) of each liability
#'   payment.
#' @param asset_cashflows A list where each element is a list with components
#'   \code{$payment} and \code{$time}, defining each asset's cash flow.
#' @param weights Numeric vector of portfolio weights (amount invested in each
#'   asset). Must have same length as \code{asset_cashflows}.
#' @param i0 Base effective interest rate per period.
#' @param delta A numeric value defining the range of rates: from
#'   \code{i0 - delta} to \code{i0 + delta}.
#' @param n_grid Number of rate values to evaluate.
#'
#' @return A \code{ggplot2} object showing the PV difference curve
#'   \eqn{PV_A(i) - PV_L(i)}{PV_A(i) - PV_L(i)} and a zero reference line.
#'
#' @details
#' Let \eqn{v(i) = 1/(1+i)}{v(i) = 1/(1+i)}. For a liability stream
#' \eqn{L_k}{L_k} at time \eqn{t_k}{t_k}:
#' \deqn{PV_L(i) = \sum_k L_k \, v(i)^{t_k}}{PV_L(i) = sum(L_k * v(i)^t_k)}
#'
#' For a portfolio of assets with weights \eqn{w_j}{w_j}:
#' \deqn{PV_A(i) = \sum_j w_j \, PV_j(i)}{PV_A(i) = sum(w_j * PV_j(i))}
#'
#' The curve \eqn{\Delta(i) = PV_A(i) - PV_L(i)}{Delta(i) = PV_A(i) - PV_L(i)}
#' illustrates immunization robustness. Under perfect duration immunization,
#' this curve is tangent to zero at \eqn{i = i_0}{i = i0} and non-negative
#' nearby if the convexity condition is also met.
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
#'   liabilities = c(5000, 8000),
#'   t_liabilities = c(3, 7),
#'   asset_cashflows = list(
#'     list(payment = c(0, 0, 100), time = c(1, 2, 3)),
#'     list(payment = c(0, 0, 0, 0, 0, 0, 200), time = 1:7)
#'   ),
#'   weights = c(5, 2.5),
#'   i0 = 0.05,
#'   delta = 0.02
#' )
#'
#' @export
plot_immunization_gap <- function(
    liabilities,
    t_liabilities,
    asset_cashflows,
    weights,
    i0,
    delta = 0.01,
    n_grid = 200L
) {
  if (length(asset_cashflows) != length(weights)) {
    stop("'asset_cashflows' and 'weights' must have same length.", call. = FALSE)
  }

  i_grid <- seq(i0 - delta, i0 + delta, length.out = n_grid)

  pv_liab <- function(i) {
    v <- 1 / (1 + i)
    sum(liabilities * v^t_liabilities)
  }

  pv_asset_j <- function(cf, i) {
    v <- 1 / (1 + i)
    sum(cf$payment * v^cf$time)
  }

  pv_portfolio <- function(i) {
    sum(vapply(seq_along(asset_cashflows), function(j) {
      weights[j] * pv_asset_j(asset_cashflows[[j]], i)
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
      x = "Interest rate i",
      y = "PV difference"
    ) +
    ggplot2::theme_minimal()
}
