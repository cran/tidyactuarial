#' Gross (expense-loaded) premium from a net premium
#'
#' Adjusts a net premium using a simple expense structure
#' \eqn{(\alpha, \beta, \gamma)} to obtain the gross or commercial premium
#' through the extended equivalence principle, using compact actuarial notation.
#'
#' The function is designed to work with the detailed output of
#' \code{\link{premium_x}} or \code{\link{premium_xy}} when those functions
#' return a one-row tibble containing the APV of the premium annuity. It can
#' also be used with any compatible one-row tibble.
#'
#' @param prem A one-row data frame or tibble containing at least:
#'   \itemize{
#'     \item \code{premium} or \code{P}: net premium per payment.
#'     \item \code{apv_premiums} or \code{a_premiums}: APV of the premium
#'     annuity.
#'   }
#' @param alpha Numeric scalar greater than or equal to 0. Initial acquisition
#'   expense as a multiple of one gross premium payment. The initial expense is
#'   \eqn{\alpha G}, paid once at issue.
#' @param beta Numeric scalar in \eqn{[0,1)}. Proportional collection expense
#'   as a fraction of each gross premium payment.
#' @param gamma Numeric scalar greater than or equal to 0. Fixed maintenance
#'   expense per premium payment period, in monetary units.
#' @param tidy Logical scalar. If \code{FALSE}, returns the gross premium as a
#'   numeric value. If \code{TRUE}, returns a one-row tibble with the expense
#'   breakdown.
#' @param ... Transitional compatibility for older calls using
#'   \code{output = "value"} or \code{output = "table"}. This argument will be
#'   removed in a future version.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}. The gross premium is denoted by \code{G}, the net
#' premium by \code{P_net}, and the APV of the premium annuity by
#' \code{a_premiums}.
#'
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
#' \code{apv_premiums} or \code{a_premiums} column of \code{prem}.
#'
#' @return
#' If \code{tidy = FALSE}, a numeric gross premium per payment.
#'
#' If \code{tidy = TRUE}, a one-row tibble with columns \code{G},
#' \code{P_net}, \code{alpha}, \code{beta}, \code{gamma}, \code{loading_pct},
#' and \code{a_premiums}.
#'
#' @seealso \code{\link{premium_x}} for single-life net premiums,
#'   \code{\link{premium_xy}} for two-life net premiums,
#'   \code{\link{annuity_x}} for building custom expense APVs.
#'
#' @family life-contingencies
#'
#' @examples
#' # Compatible one-row premium table
#' net <- tibble::tibble(
#'   premium = 1200,
#'   apv_premiums = 12.5
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
#'   tidy = TRUE
#' )
#'
#' # Also accepts compact actuarial column names
#' net2 <- tibble::tibble(
#'   P = 1200,
#'   a_premiums = 12.5
#' )
#'
#' premium_gross(net2, alpha = 0.5, beta = 0.05, gamma = 50)
#'
#' @export
premium_gross <- function(
    prem,
    alpha = 0,
    beta = 0,
    gamma = 0,
    tidy = FALSE,
    ...
) {
  dots <- list(...)

  if (length(dots) > 0L) {
    allowed_old <- "output"
    bad_dots <- setdiff(names(dots), allowed_old)

    if (length(bad_dots) > 0L) {
      stop(
        "Unused argument(s): ",
        paste(sprintf("`%s`", bad_dots), collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    if (!is.null(dots$output)) {
      if (!identical(tidy, FALSE)) {
        stop("Provide only one of `tidy` or deprecated `output`.", call. = FALSE)
      }

      output <- match.arg(dots$output, c("value", "table"))
      tidy <- identical(output, "table")
    }
  }

  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  # --- Input checks ---
  if (!inherits(prem, "data.frame") || nrow(prem) != 1L) {
    stop(
      "`prem` must be a one-row tibble from ",
      "premium_x(..., tidy = TRUE), premium_xy(..., tidy = TRUE), ",
      "or another compatible premium table.",
      call. = FALSE
    )
  }

  # Accept both old and compact column names while the package is being migrated.
  premium_col <- if ("P" %in% names(prem)) {
    "P"
  } else if ("premium" %in% names(prem)) {
    "premium"
  } else {
    NA_character_
  }

  apv_col <- if ("a_premiums" %in% names(prem)) {
    "a_premiums"
  } else if ("apv_premiums" %in% names(prem)) {
    "apv_premiums"
  } else {
    NA_character_
  }

  if (is.na(premium_col)) {
    stop("Column `P` or `premium` was not found in `prem`.", call. = FALSE)
  }

  if (is.na(apv_col)) {
    stop("Column `a_premiums` or `apv_premiums` was not found in `prem`.", call. = FALSE)
  }

  P_net <- prem[[premium_col]][[1]]
  a_premiums <- prem[[apv_col]][[1]]

  if (!is.numeric(P_net) ||
      length(P_net) != 1L ||
      is.na(P_net) ||
      !is.finite(P_net) ||
      P_net <= 0) {
    stop("The net premium must be a single positive finite number.", call. = FALSE)
  }

  if (!is.numeric(a_premiums) ||
      length(a_premiums) != 1L ||
      is.na(a_premiums) ||
      !is.finite(a_premiums) ||
      a_premiums <= 0) {
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

  # G = (P_net + gamma) / ((1 - beta) - alpha / a_premiums)
  denominator <- (1 - beta) - alpha / a_premiums

  if (!is.finite(denominator) || denominator <= 0) {
    stop(
      "No level gross premium exists: the expense structure exceeds ",
      "the premium annuity capacity.",
      call. = FALSE
    )
  }

  G <- (P_net + gamma) / denominator

  if (!tidy) {
    return(G)
  }

  tibble::tibble(
    G = G,
    P_net = P_net,
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    loading_pct = (G - P_net) / P_net * 100,
    a_premiums = a_premiums
  )
}
