#' Compute Monte Carlo loss random variables for life contingencies
#'
#' Computes simulated actuarial loss random variables from Monte Carlo present
#' values of benefits, premium annuities, and premiums, using compact actuarial
#' notation.
#'
#' This function constructs the simulated loss random variable
#'
#' \deqn{
#'   L = Z - P Y,
#' }
#'
#' where \eqn{Z} is the present value random variable of insurance benefits,
#' \eqn{Y} is the present value random variable of the premium annuity, and
#' \eqn{P} is the premium per payment.
#'
#' @param .data A data frame or tibble containing simulated present values of
#'   benefits and premium annuities.
#' @param col_Z Character string. Name of the column containing the simulated
#'   present value of benefits. Default is \code{"pv_benefit"}.
#' @param col_Y Character string. Name of the column containing the simulated
#'   present value of premium annuities. Default is \code{"pv_annuity"}.
#' @param col_P Character string. Name of the column containing the premium.
#'   Default is \code{"P"}. If this column is not found and \code{col_P = "P"},
#'   the legacy column \code{"premium"} is used when available.
#' @param col_L Character string. Name of the output column containing the
#'   simulated loss random variable. Default is \code{"L"}.
#' @param P Optional numeric scalar. If supplied, this value is used as the
#'   premium instead of reading the premium from \code{col_P}.
#' @param ... Transitional compatibility for older calls using \code{data},
#'   \code{benefit_col}, \code{annuity_col}, \code{premium_col},
#'   \code{loss_col}, and \code{premium}.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{Z} denotes the present value random variable of
#' benefits, \code{Y} denotes the present value random variable of the premium
#' annuity, \code{P} denotes the premium per payment, and \code{L} denotes the
#' actuarial loss random variable at issue.
#'
#' The function does not estimate the premium. It only constructs the simulated
#' loss random variable. The premium may come from a column previously created
#' by \code{\link{mc_premium}}, or it may be supplied directly through the
#' \code{P} argument.
#'
#' The interpretation of \code{P} depends on how the premium annuity present
#' value was constructed. If \code{pv_annuity} was generated with annual
#' payments, then \code{P} corresponds to that annual premium structure. If
#' \code{pv_annuity} was generated using fractional payments in
#' \code{\link{mc_annuity}}, such as \code{payment = 1 / 12} and
#' \code{k = 12}, then \code{P} is applied to that same fractional payment
#' pattern.
#'
#' Thus, \code{mc_loss()} can be used without modification for annual premiums,
#' monthly premiums, quarterly premiums, semiannual premiums, and multiple-life
#' premium structures, provided that \code{col_Y} contains the appropriate
#' simulated premium annuity present value.
#'
#' The resulting loss column can be used to estimate quantities such as expected
#' loss, variance, standard deviation, probability of positive loss, loss
#' quantiles, empirical value-at-risk measures, and sensitivities across
#' actuarial assumptions.
#'
#' @return A tibble with the original columns and one additional column
#' containing the simulated loss random variable. The name of this column is
#' controlled by \code{col_L}. For transition, if \code{col_L != "loss"} and
#' the input does not already contain a column named \code{loss}, a legacy
#' column \code{loss} is also added with the same value.
#'
#' @seealso
#' \code{\link{simulate_lifetime}}, \code{\link{simulate_lifetimes}},
#' \code{\link{mc_multilife_status}}, \code{\link{mc_insurance}},
#' \code{\link{mc_annuity}}, \code{\link{mc_premium}},
#' \code{\link{mc_reserve}}, \code{\link{summary_mc}}
#'
#' @references
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A.,
#' and Nesbitt, C. J. (1997). \emph{Actuarial Mathematics}. Second Edition.
#' Society of Actuaries.
#'
#' @family monte-carlo
#'
#' @examples
#' # Example 1: direct use with simulated present values
#' sim_values <- tibble::tibble(
#'   sim_id = 1:5,
#'   pv_benefit = c(0.80, 0.75, 0.95, 0.60, 0.70),
#'   pv_annuity = c(8.0, 7.5, 9.0, 6.0, 7.0),
#'   P = 0.10
#' )
#'
#' sim_values |>
#'   mc_loss()
#'
#' # Example 2: using a premium supplied directly
#' sim_values |>
#'   mc_loss(P = 0.12)
#'
#' # Example 3: annual whole life loss
#' lt <- tibble::tibble(
#'   x = 40:100,
#'   qx = seq(0.002, 1, length.out = 61)
#' )
#'
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     seed = 123
#'   ) |>
#'   mc_insurance(
#'     i = 0.05,
#'     type = "whole",
#'     benefit = 1
#'   ) |>
#'   mc_annuity(
#'     i = 0.05,
#'     type = "whole",
#'     payment = 1,
#'     k = 1,
#'     timing = "due"
#'   ) |>
#'   mc_premium() |>
#'   mc_loss()
#'
#' # Example 4: monthly whole life loss
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     frac = "udd",
#'     seed = 123
#'   ) |>
#'   mc_insurance(
#'     i = 0.05,
#'     type = "whole",
#'     benefit = 1
#'   ) |>
#'   mc_annuity(
#'     i = 0.05,
#'     type = "whole",
#'     payment = 1 / 12,
#'     k = 12,
#'     timing = "due"
#'   ) |>
#'   mc_premium() |>
#'   mc_loss()
#'
#' # Example 5: summarising simulated losses
#' lt |>
#'   simulate_lifetime(
#'     x = 45,
#'     n_sim = 25,
#'     frac = "udd",
#'     seed = 123
#'   ) |>
#'   mc_insurance(
#'     i = 0.04,
#'     type = "term",
#'     n = 20,
#'     benefit = 100000
#'   ) |>
#'   mc_annuity(
#'     i = 0.04,
#'     type = "temporary",
#'     n = 20,
#'     payment = 1 / 12,
#'     k = 12,
#'     timing = "due"
#'   ) |>
#'   mc_premium() |>
#'   mc_loss() |>
#'   summary_mc(value_col = "L")
#'
#' # Example 6: joint-life annual loss
#' lt |>
#'   simulate_lifetimes(
#'     x = c(60, 58),
#'     n_sim = 25,
#'     seed = 123
#'   ) |>
#'   mc_multilife_status(status = "joint") |>
#'   mc_insurance(
#'     i = 0.04,
#'     type = "whole",
#'     benefit = 100000,
#'     col_K = "K_status",
#'     col_T = "T_status"
#'   ) |>
#'   mc_annuity(
#'     i = 0.04,
#'     type = "whole",
#'     payment = 1,
#'     k = 1,
#'     timing = "due",
#'     col_K = "K_status",
#'     col_T = "T_status"
#'   ) |>
#'   mc_premium() |>
#'   mc_loss()
#'
#' # Transitional compatibility with old column arguments
#' sim_values |>
#'   mc_loss(
#'     benefit_col = "pv_benefit",
#'     annuity_col = "pv_annuity",
#'     premium_col = "P",
#'     loss_col = "loss"
#'   )
#'
#' @export
mc_loss <- function(
    .data = NULL,
    col_Z = "pv_benefit",
    col_Y = "pv_annuity",
    col_P = "P",
    col_L = "L",
    P = NULL,
    ...
) {
  dots <- list(...)

  # Use exact name matching for deprecated arguments.
  # Do not use `dots$premium`, because `$` partially matches `premium_col`.
  dot_has <- function(nm) {
    nm %in% names(dots)
  }

  dot_get <- function(nm) {
    dots[[nm]]
  }

  # -------------------------------------------------------------------------
  # Transitional compatibility with the previous public API
  # -------------------------------------------------------------------------

  allowed_old <- c(
    "data",
    "benefit_col",
    "annuity_col",
    "premium_col",
    "loss_col",
    "premium"
  )

  bad_dots <- setdiff(names(dots), allowed_old)

  if (length(bad_dots) > 0L) {
    stop(
      "Unused argument(s): ",
      paste(sprintf("`%s`", bad_dots), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (dot_has("data")) {
    if (!is.null(.data)) {
      stop("Provide only one of `.data` or deprecated `data`.", call. = FALSE)
    }

    .data <- dot_get("data")
  }

  if (dot_has("benefit_col")) {
    if (!identical(col_Z, "pv_benefit")) {
      stop("Provide only one of `col_Z` or deprecated `benefit_col`.",
           call. = FALSE)
    }

    col_Z <- dot_get("benefit_col")
  }

  if (dot_has("annuity_col")) {
    if (!identical(col_Y, "pv_annuity")) {
      stop("Provide only one of `col_Y` or deprecated `annuity_col`.",
           call. = FALSE)
    }

    col_Y <- dot_get("annuity_col")
  }

  if (dot_has("premium_col")) {
    if (!identical(col_P, "P")) {
      stop("Provide only one of `col_P` or deprecated `premium_col`.",
           call. = FALSE)
    }

    col_P <- dot_get("premium_col")
  }

  if (dot_has("loss_col")) {
    if (!identical(col_L, "L")) {
      stop("Provide only one of `col_L` or deprecated `loss_col`.",
           call. = FALSE)
    }

    col_L <- dot_get("loss_col")
  }

  if (dot_has("premium")) {
    if (!is.null(P)) {
      stop("Provide only one of `P` or deprecated `premium`.", call. = FALSE)
    }

    P <- dot_get("premium")
  }

  # -------------------------------------------------------------------------
  # Validation
  # -------------------------------------------------------------------------

  if (!is.data.frame(.data)) {
    stop("`.data` must be a data frame or tibble.", call. = FALSE)
  }

  .mc_assert_numeric_column(.data, col_Z, "col_Z")
  .mc_assert_numeric_column(.data, col_Y, "col_Y")
  .mc_assert_character_scalar(col_P, "col_P")
  .mc_assert_character_scalar(col_L, "col_L")

  if (col_L %in% c(col_Z, col_Y, col_P)) {
    stop(
      "`col_L` must be different from `col_Z`, `col_Y`, and `col_P`.",
      call. = FALSE
    )
  }

  Z <- .data[[col_Z]]
  Y <- .data[[col_Y]]

  if (!is.null(P)) {
    .mc_assert_numeric_scalar(P, "P")
    P_value <- rep(P, length(Z))
  } else {
    if (!col_P %in% names(.data) &&
        identical(col_P, "P") &&
        "premium" %in% names(.data)) {
      col_P <- "premium"
    }

    .mc_assert_numeric_column(.data, col_P, "col_P")
    P_value <- .data[[col_P]]
  }

  L <- Z - P_value * Y

  out <- .data |>
    dplyr::mutate(
      "{col_P}" := P_value,
      "{col_L}" := L
    )

  if (!identical(col_L, "loss") && !"loss" %in% names(out)) {
    out <- out |>
      dplyr::mutate(
        loss = L
      )
  }

  out
}

