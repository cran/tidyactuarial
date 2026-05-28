#' Compute Monte Carlo net premiums for life contingencies
#'
#' Computes simulated net premiums from Monte Carlo present values of insurance
#' benefits and premium annuities.
#'
#' This function applies the actuarial equivalence principle to simulated
#' present value random variables. If \eqn{Z} denotes the present value random
#' variable of the benefit and \eqn{Y} denotes the present value random variable
#' of the premium annuity, the net premium is estimated as
#'
#' \deqn{
#'   \hat{P} = \frac{\bar{Z}}{\bar{Y}}.
#' }
#'
#' Equivalently, this is the Monte Carlo estimator of
#'
#' \deqn{
#'   P = \frac{E[Z]}{E[Y]}.
#' }
#'
#' @param data A data frame or tibble containing simulated present values.
#'   Usually this object is obtained after applying [mc_insurance()] and
#'   [mc_annuity()] to the same simulated lifetime sample.
#' @param benefit_col Character string. Name of the column containing the
#'   simulated present value of the insurance benefit. Default is
#'   `"pv_benefit"`.
#' @param annuity_col Character string. Name of the column containing the
#'   simulated present value of the premium annuity. Default is `"pv_annuity"`.
#' @param premium_col Character string. Name of the output column containing
#'   the simulated net premium. Default is `"premium"`.
#' @param by Optional character vector with grouping columns. If supplied, the
#'   premium is computed separately within each group. If `by = NULL` and
#'   `data` is already grouped with [dplyr::group_by()], the current grouping
#'   structure is used.
#' @param na_rm Logical. Should missing values be removed when computing
#'   simulated means? Default is `TRUE`.
#'
#' @details
#' The function does not simulate lifetimes and does not calculate present
#' values directly. It only computes the Monte Carlo net premium from columns
#' that already contain simulated present values.
#'
#' In a typical workflow, [simulate_lifetime()] generates simulated values of
#' \eqn{K_x} and possibly \eqn{T_x}; [mc_insurance()] creates the simulated
#' benefit present value \eqn{Z}; [mc_annuity()] creates the simulated premium
#' annuity present value \eqn{Y}; and `mc_premium()` estimates the net level
#' premium.
#'
#' The estimated premium is attached to every row of the input data. This is
#' intentional: it makes it easy to construct the simulated loss random variable
#' with [mc_loss()], for example
#'
#' \deqn{
#'   L = Z - \hat{P}Y.
#' }
#'
#' The function is also valid for premiums payable more than once per year.
#' In that case, the payment frequency is not specified in `mc_premium()`;
#' it is already embedded in the simulated premium annuity present value
#' supplied through `annuity_col`.
#'
#' For example, if [mc_annuity()] was called with
#' `payment = 1 / 12` and `payments_per_year = 12`, then `pv_annuity`
#' represents the present value of a monthly premium stream whose total annual
#' amount is 1. The premium estimated by `mc_premium()` is then consistent with
#' that monthly payment structure and corresponds to a Monte Carlo estimate of
#' a premium such as \eqn{P_x^{(12)}}.
#'
#' If [mc_annuity()] was called with `payment = 1` and
#' `payments_per_year = 12`, then `pv_annuity` represents a stream of payments
#' of 1 each month, and the resulting premium should be interpreted relative to
#' that payment pattern.
#'
#' This function computes net premiums only. It does not include expenses,
#' safety loadings, profit margins, taxes, surrender charges, commissions, or
#' other practical pricing adjustments.
#'
#' @return A tibble with the original columns and one additional column
#' containing the simulated net premium. The name of this column is controlled
#' by `premium_col`.
#'
#' @seealso
#' [simulate_lifetime()], [simulate_lifetimes()], [mc_insurance()],
#' [mc_annuity()], [mc_loss()], [mc_reserve()], [summary_mc()]
#'
#' @references
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A.,
#' and Nesbitt, C. J. (1997). *Actuarial Mathematics*. Second Edition.
#' Society of Actuaries.
#'
#' @examples
#' # Example 1: direct use with simulated present values
#' sim_values <- tibble::tibble(
#'   sim_id = 1:6,
#'   pv_benefit = c(0.82, 0.74, 0.61, 0.95, 0.70, 0.88),
#'   pv_annuity = c(8.2, 7.5, 6.1, 9.0, 7.2, 8.8)
#' )
#'
#' sim_values |>
#'   mc_premium()
#'
#' # Example 2: grouped premiums by age
#' sim_by_age <- tibble::tibble(
#'   sim_id = rep(1:6, times = 2),
#'   age = rep(c(40, 50), each = 6),
#'   pv_benefit = c(
#'     0.82, 0.74, 0.61, 0.95, 0.70, 0.88,
#'     0.91, 0.86, 0.79, 0.98, 0.83, 0.94
#'   ),
#'   pv_annuity = c(
#'     8.2, 7.5, 6.1, 9.0, 7.2, 8.8,
#'     6.8, 6.4, 5.9, 7.1, 6.2, 6.7
#'   )
#' )
#'
#' sim_by_age |>
#'   mc_premium(by = "age")
#'
#' # Example 3: using dplyr grouping
#' sim_by_age |>
#'   dplyr::group_by(age) |>
#'   mc_premium()
#'
#' # Example 4: annual whole life net premium
#' life_table <- tibble::tibble(
#'   age = 40:100,
#'   qx = seq(0.002, 1, length.out = 61)
#' )
#'
#' life_table |>
#'   simulate_lifetime(age = 40, n_sim = 500, seed = 123) |>
#'   mc_insurance(
#'     rate = 0.05,
#'     insurance = "whole_life",
#'     benefit = 1
#'   ) |>
#'   mc_annuity(
#'     rate = 0.05,
#'     annuity = "whole_life",
#'     payment = 1,
#'     payments_per_year = 1,
#'     timing = "due"
#'   ) |>
#'   mc_premium()
#'
#' # Example 5: monthly whole life net premium
#' life_table |>
#'   simulate_lifetime(
#'     age = 40,
#'     n_sim = 500,
#'     fractional = "udd",
#'     seed = 123
#'   ) |>
#'   mc_insurance(
#'     rate = 0.05,
#'     insurance = "whole_life",
#'     benefit = 1
#'   ) |>
#'   mc_annuity(
#'     rate = 0.05,
#'     annuity = "whole_life",
#'     payment = 1 / 12,
#'     payments_per_year = 12,
#'     timing = "due"
#'   ) |>
#'   mc_premium()
#'
#' # Example 6: term insurance with monthly premium annuity
#' life_table |>
#'   simulate_lifetime(
#'     age = 45,
#'     n_sim = 500,
#'     fractional = "udd",
#'     seed = 123
#'   ) |>
#'   mc_insurance(
#'     rate = 0.04,
#'     insurance = "term",
#'     term = 20,
#'     benefit = 100000
#'   ) |>
#'   mc_annuity(
#'     rate = 0.04,
#'     annuity = "temporary",
#'     term = 20,
#'     payment = 1 / 12,
#'     payments_per_year = 12,
#'     timing = "due"
#'   ) |>
#'   mc_premium()
#'
#' # Example 7: joint-life monthly net premium
#' life_table |>
#'   simulate_lifetimes(
#'     ages = c(60, 58),
#'     n_sim = 500,
#'     fractional = "udd",
#'     seed = 123
#'   ) |>
#'   mc_multilife_status(status = "joint_life") |>
#'   mc_insurance(
#'     rate = 0.04,
#'     insurance = "whole_life",
#'     benefit = 100000,
#'     k_col = "K_status",
#'     tx_col = "T_status"
#'   ) |>
#'   mc_annuity(
#'     rate = 0.04,
#'     annuity = "whole_life",
#'     payment = 1 / 12,
#'     payments_per_year = 12,
#'     timing = "due",
#'     k_col = "K_status",
#'     tx_col = "T_status"
#'   ) |>
#'   mc_premium()
#'
#' @export
mc_premium <- function(data,
                       benefit_col = "pv_benefit",
                       annuity_col = "pv_annuity",
                       premium_col = "premium",
                       by = NULL,
                       na_rm = TRUE) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or tibble.", call. = FALSE)
  }

  .mc_assert_numeric_column(data, benefit_col, "benefit_col")
  .mc_assert_numeric_column(data, annuity_col, "annuity_col")
  .mc_assert_character_scalar(premium_col, "premium_col")

  if (!is.logical(na_rm) || length(na_rm) != 1 || is.na(na_rm)) {
    stop("`na_rm` must be a logical scalar.", call. = FALSE)
  }

  if (!is.null(by)) {
    if (!is.character(by) || anyNA(by)) {
      stop("`by` must be NULL or a character vector.", call. = FALSE)
    }

    if (!all(by %in% names(data))) {
      stop("All columns supplied in `by` must exist in `data`.", call. = FALSE)
    }
  }

  active_by <- by

  if (is.null(active_by)) {
    active_by <- dplyr::group_vars(data)
  }

  compute_premium <- function(benefit, annuity) {
    denominator <- mean(annuity, na.rm = na_rm)

    if (is.na(denominator) || denominator <= 0) {
      stop(
        "The mean simulated annuity present value must be positive.",
        call. = FALSE
      )
    }

    numerator <- mean(benefit, na.rm = na_rm)

    if (is.na(numerator)) {
      stop(
        "The mean simulated benefit present value could not be computed.",
        call. = FALSE
      )
    }

    numerator / denominator
  }

  if (length(active_by) == 0) {
    premium_value <- compute_premium(
      benefit = data[[benefit_col]],
      annuity = data[[annuity_col]]
    )

    return(
      data |>
        dplyr::ungroup() |>
        dplyr::mutate(
          "{premium_col}" := premium_value
        )
    )
  }

  data |>
    dplyr::ungroup() |>
    dplyr::group_by(dplyr::across(dplyr::all_of(active_by))) |>
    dplyr::mutate(
      "{premium_col}" := compute_premium(
        benefit = .data[[benefit_col]],
        annuity = .data[[annuity_col]]
      )
    ) |>
    dplyr::ungroup()
}
