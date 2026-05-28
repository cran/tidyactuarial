#' Gross (expense-loaded) premium from net premium
#'
#' Adjusts a net premium using a simple expense structure
#' \eqn{(\alpha, \beta, \gamma)} to obtain the gross or commercial premium
#' through the extended equivalence principle.
#'
#' The function is designed to work with the detailed output of
#' \code{\link{premium_x}} using \code{output = "table"}. It can also be used
#' with any one-row tibble containing the columns \code{premium} and
#' \code{apv_premiums}.
#'
#' @param prem A one-row data frame or tibble containing at least:
#'   \itemize{
#'     \item \code{premium}: net premium per payment.
#'     \item \code{apv_premiums}: APV of the premium annuity.
#'   }
#' @param alpha Numeric scalar greater than or equal to 0. Initial acquisition
#'   expense as a multiple of one gross premium payment. The initial expense is
#'   \eqn{\alpha G}, paid once at issue.
#' @param beta Numeric scalar in \eqn{[0,1)}. Proportional collection expense
#'   as a fraction of each gross premium payment.
#' @param gamma Numeric scalar greater than or equal to 0. Fixed maintenance
#'   expense per premium payment period, in monetary units.
#' @param output Character string. Use \code{"value"} to return a numeric gross
#'   premium, or \code{"table"} to return a one-row tibble with the expense
#'   breakdown.
#'
#' @details
#' The extended equivalence principle equates the APV of gross premiums with
#' the APV of benefits plus expenses:
#' \deqn{
#' G \ddot{a} =
#' P_{\text{net}}\ddot{a}
#' + \alpha G
#' + \beta G \ddot{a}
#' + \gamma \ddot{a}.
#' }
#'
#' Solving for the gross premium gives:
#' \deqn{
#' G =
#' \frac{P_{\text{net}} + \gamma}
#' {(1-\beta) - \alpha / \ddot{a}}.
#' }
#'
#' In this function, \eqn{\ddot{a}} is supplied through the
#' \code{apv_premiums} column of \code{prem}.
#'
#' @return
#' If \code{output = "value"}, a numeric gross premium per payment.
#'
#' If \code{output = "table"}, a one-row tibble with columns
#' \code{gross_premium}, \code{net_premium}, \code{alpha}, \code{beta},
#' \code{gamma}, \code{loading_pct}, and \code{apv_premiums}.
#'
#' @seealso \code{\link{premium_x}} for single-life net premiums,
#'   \code{\link{premium_xy}} for two-life net premiums,
#'   \code{\link{annuity_x}} for building custom expense APVs.
#'
#' @family life-contingencies
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:66,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000, 86000)
#' )
#'
#' # Full workflow: net premium -> gross premium
#' net <- premium_x(
#'   mortality_table = lt,
#'   age = 60,
#'   rate = 0.05,
#'   insurance_type = "whole",
#'   benefit = 100000,
#'   output = "table"
#' )
#'
#' premium_gross(net, alpha = 0.5, beta = 0.05, gamma = 50)
#'
#' # Finan-style expense structure:
#' # 10% of each premium plus fixed expenses of 275 per payment period
#' premium_gross(net, alpha = 0, beta = 0.10, gamma = 275)
#'
#' # Detailed output with expense breakdown
#' premium_gross(
#'   net,
#'   alpha = 0.5,
#'   beta = 0.05,
#'   gamma = 50,
#'   output = "table"
#' )
#'
#' # No expenses: gross = net
#' premium_gross(net)
#'
#' @export
premium_gross <- function(
    prem,
    alpha = 0,
    beta = 0,
    gamma = 0,
    output = c("value", "table")
) {
  output <- match.arg(output)

  # --- Input checks ---
  if (!inherits(prem, "data.frame") || nrow(prem) != 1L) {
    stop(
      "'prem' must be a one-row tibble from ",
      "premium_x(..., output = \"table\") or another compatible premium table.",
      call. = FALSE
    )
  }

  if (!("premium" %in% names(prem))) {
    stop("Column `premium` was not found in `prem`.", call. = FALSE)
  }

  if (!("apv_premiums" %in% names(prem))) {
    stop("Column `apv_premiums` was not found in `prem`.", call. = FALSE)
  }

  P_net <- prem$premium[[1]]
  apv_premiums <- prem$apv_premiums[[1]]

  if (!is.numeric(P_net) ||
      length(P_net) != 1L ||
      is.na(P_net) ||
      !is.finite(P_net) ||
      P_net <= 0) {
    stop("The net premium must be a single positive finite number.", call. = FALSE)
  }

  if (!is.numeric(apv_premiums) ||
      length(apv_premiums) != 1L ||
      is.na(apv_premiums) ||
      !is.finite(apv_premiums) ||
      apv_premiums <= 0) {
    stop(
      "The APV of premiums must be a single positive finite number.",
      call. = FALSE
    )
  }

  if (!is.numeric(alpha) ||
      length(alpha) != 1L ||
      is.na(alpha) ||
      !is.finite(alpha) ||
      alpha < 0) {
    stop("`alpha` must be a single finite nonnegative number.", call. = FALSE)
  }

  if (!is.numeric(beta) ||
      length(beta) != 1L ||
      is.na(beta) ||
      !is.finite(beta) ||
      beta < 0 ||
      beta >= 1) {
    stop("`beta` must be a single finite number in [0, 1).", call. = FALSE)
  }

  if (!is.numeric(gamma) ||
      length(gamma) != 1L ||
      is.na(gamma) ||
      !is.finite(gamma) ||
      gamma < 0) {
    stop("`gamma` must be a single finite nonnegative number.", call. = FALSE)
  }

  # G = (P_net + gamma) / ((1 - beta) - alpha / apv_premiums)
  denominator <- (1 - beta) - alpha / apv_premiums

  if (!is.finite(denominator) || denominator <= 0) {
    stop(
      "No level gross premium exists: the expense structure exceeds ",
      "the premium annuity capacity.",
      call. = FALSE
    )
  }

  gross_premium <- (P_net + gamma) / denominator

  if (output == "value") {
    return(gross_premium)
  }

  tibble::tibble(
    gross_premium = gross_premium,
    net_premium = P_net,
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    loading_pct = (gross_premium - P_net) / P_net * 100,
    apv_premiums = apv_premiums
  )
}
