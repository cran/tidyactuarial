#' Gross (expense-loaded) premium from net premium
#'
#' Adjusts a net premium using a simple expense structure
#' (\eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}) to obtain the gross
#' (commercial) premium via the extended equivalence principle
#' (Finan, Sections 70--71).
#'
#' Designed to work directly with the output of
#' \code{\link{premium_x}} or \code{\link{premium_xy}} when
#' \code{tidy = TRUE}.
#'
#' @param prem A one-row tibble returned by
#'   \code{premium_x(..., tidy = TRUE)} or
#'   \code{premium_xy(..., tidy = TRUE)}. Must contain columns
#'   \code{premium} (net premium per payment) and \code{apv_premiums}
#'   (APV of the premium annuity).
#' @param alpha Numeric \eqn{\ge 0}. Initial (acquisition) expense as
#'   a multiple of one gross premium payment \eqn{G}. The initial
#'   expense is \eqn{\alpha \cdot G}, paid once at issue.
#'   Default \code{0}.
#' @param beta Numeric in \eqn{[0, 1)}. Proportional (collection)
#'   expense as a fraction of each gross premium payment. Each period
#'   the insurer incurs \eqn{\beta \cdot G}. Default \code{0}.
#' @param gamma Numeric \eqn{\ge 0}. Fixed maintenance expense per
#'   premium payment period (in monetary units). Default \code{0}.
#' @param tidy Logical. If \code{TRUE}, returns a one-row tibble with
#'   the gross premium and expense breakdown.
#'
#' @details
#' The extended equivalence principle (Finan, Section 70) equates the
#' APV of gross premiums with the APV of benefits plus expenses:
#' \deqn{G \cdot \ddot{a} = P_{\text{net}} \cdot \ddot{a} +
#'   \alpha \cdot G + \beta \cdot G \cdot \ddot{a} +
#'   \gamma \cdot \ddot{a}.}
#'
#' Solving for \eqn{G}:
#' \deqn{G = \frac{P_{\text{net}} + \gamma}{(1 - \beta) -
#'   \alpha / \ddot{a}}}
#' where \eqn{\ddot{a}} is the APV of the premium annuity
#' (\code{apv_premiums} column in \code{prem}).
#'
#' The expense structure maps to common actuarial categories
#' (Finan, Section 71):
#' \itemize{
#'   \item \eqn{\alpha}: acquisition costs (agent commission,
#'     underwriting) - higher in year 1
#'   \item \eqn{\beta}: collection costs - proportional to premium
#'   \item \eqn{\gamma}: maintenance costs - fixed per period
#' }
#'
#' For more complex expense structures (different first-year vs.
#' renewal rates, per-policy vs. per-thousand, settlement expenses),
#' the user should build the expense APV manually using
#' \code{\link{annuity_x}} and apply the equivalence principle
#' directly (see Finan, Examples 70.2, 71.1--71.3).
#'
#' @return A numeric gross premium per payment, or a one-row tibble
#'   if \code{tidy = TRUE} with columns \code{gross_premium},
#'   \code{net_premium}, \code{alpha}, \code{beta}, \code{gamma},
#'   \code{loading_pct}.
#'
#' @seealso \code{\link{premium_x}} for single-life net premiums,
#'   \code{\link{premium_xy}} for two-life net premiums,
#'   \code{\link{annuity_x}} for building custom expense APVs.
#'
#' @examples
#' lt <- data.frame(
#'   x  = 60:66,
#'   lx = c(100000, 99000, 97500, 95500, 93000, 90000, 86000)
#' )
#'
#' # Full workflow: net premium -> gross premium
#' net <- premium_x(lt, x = 60, i = 0.05,
#'                  product = "whole", benefit = 100000,
#'                  tidy = TRUE)
#' premium_gross(net, alpha = 0.5, beta = 0.05, gamma = 50)
#'
#' # Finan Example 70.3 style: 10% of premium + $25/yr + $250/yr
#' # beta = 0.10, gamma = 275 (= 25 + 250), alpha = 0
#' premium_gross(net, alpha = 0, beta = 0.10, gamma = 275)
#'
#' # Tidy output with expense breakdown
#' premium_gross(net, alpha = 0.5, beta = 0.05, gamma = 50,
#'               tidy = TRUE)
#'
#' # Two-life workflow
#' net_xy <- premium_xy(lt, x = 60, y = 62, i = 0.05,
#'                      type = "whole", cohort = "first",
#'                      benefit = 100000, tidy = TRUE)
#' premium_gross(net_xy, alpha = 0.3, beta = 0.04, gamma = 30)
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
    tidy = FALSE
) {
  # --- checks ---
  if (!inherits(prem, "data.frame") || nrow(prem) != 1) {
    stop(
      "'prem' must be a one-row tibble from ",
      "premium_x() or premium_xy() with tidy = TRUE."
    )
  }
  if (!("premium" %in% names(prem))) {
    stop("Column 'premium' not found in 'prem'.")
  }
  if (!("apv_premiums" %in% names(prem))) {
    stop("Column 'apv_premiums' not found in 'prem'.")
  }

  P_net <- prem$premium[[1]]
  EY    <- prem$apv_premiums[[1]]

  if (!is.numeric(P_net) || !is.finite(P_net) || P_net <= 0) {
    stop("Net premium must be a positive finite number.")
  }
  if (!is.numeric(EY) || !is.finite(EY) || EY <= 0) {
    stop("APV of premiums must be a positive finite number.")
  }
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      is.na(alpha) || alpha < 0) {
    stop("'alpha' must be a single nonnegative number.")
  }
  if (!is.numeric(beta) || length(beta) != 1L ||
      is.na(beta) || beta < 0 || beta >= 1) {
    stop("'beta' must be a single number in [0, 1).")
  }
  if (!is.numeric(gamma) || length(gamma) != 1L ||
      is.na(gamma) || gamma < 0) {
    stop("'gamma' must be a single nonnegative number.")
  }

  # G = (P_net + gamma) / ((1 - beta) - alpha / EY)
  denom <- (1 - beta) - alpha / EY
  if (denom <= 0) {
    stop(
      "No level gross premium exists: expenses exceed ",
      "the premium annuity capacity."
    )
  }

  G <- (P_net + gamma) / denom

  if (!isTRUE(tidy)) return(G)

  tibble::tibble(
    gross_premium = G,
    net_premium   = P_net,
    alpha         = alpha,
    beta          = beta,
    gamma         = gamma,
    loading_pct   = (G - P_net) / P_net * 100,
    apv_premiums  = EY
  )
}
