#' Compute Monte Carlo loss random variables for life contingencies
#'
#' Computes simulated actuarial loss random variables from Monte Carlo present
#' values of benefits, premium annuities, and premiums.
#'
#' This function constructs the simulated loss random variable
#'
#' \deqn{
#'   L = Z - P Y,
#' }
#'
#' where \eqn{Z} is the present value random variable of the insurance benefit,
#' \eqn{Y} is the present value random variable of the premium annuity, and
#' \eqn{P} is the premium.
#'
#' @param data A data frame or tibble containing simulated present values of
#'   benefits and premium annuities.
#' @param benefit_col Character string. Name of the column containing the
#'   simulated present value of benefits. Default is `"pv_benefit"`.
#' @param annuity_col Character string. Name of the column containing the
#'   simulated present value of premium annuities. Default is `"pv_annuity"`.
#' @param premium_col Character string. Name of the column containing the
#'   premium. Default is `"premium"`.
#' @param loss_col Character string. Name of the output column containing the
#'   simulated loss random variable. Default is `"loss"`.
#' @param premium Optional numeric scalar. If supplied, this value is used as
#'   the premium instead of reading the premium from `premium_col`.
#'
#' @details
#' In actuarial notation, the loss at issue is commonly written as
#'
#' \deqn{
#'   L_0 = Z - P Y.
#' }
#'
#' The random variable \eqn{Z} represents the present value of future benefits,
#' while \eqn{Y} represents the present value of future premium payments.
#'
#' This function does not estimate the premium. It only constructs the simulated
#' loss random variable. The premium may come from a column previously created
#' by [mc_premium()], or it may be supplied directly through the `premium`
#' argument.
#'
#' The interpretation of `premium` depends on how the premium annuity present
#' value was constructed. If `pv_annuity` was generated with annual payments,
#' then `premium` corresponds to that annual premium structure. If `pv_annuity`
#' was generated using fractional payments in [mc_annuity()], such as
#' `payment = 1 / 12` and `payments_per_year = 12`, then the premium is applied
#' to that same fractional payment pattern.
#'
#' Thus, `mc_loss()` can be used without modification for annual premiums,
#' monthly premiums, quarterly premiums, semiannual premiums, and multiple-life
#' premium structures, provided that `annuity_col` contains the appropriate
#' simulated premium annuity present value.
#'
#' The resulting loss column can be used to estimate quantities such as:
#'
#' * expected loss;
#' * variance and standard deviation of loss;
#' * probability of positive loss;
#' * loss quantiles;
#' * empirical value-at-risk type measures;
#' * sensitivity of loss distributions across ages, benefits, products, payment
#'   frequencies, interest rates, or multiple-life statuses.
#'
#' @return A tibble with the original columns and one additional column
#' containing the simulated loss random variable. The name of this column is
#' controlled by `loss_col`.
#'
#' @seealso
#' [simulate_lifetime()], [simulate_lifetimes()], [mc_multilife_status()],
#' [mc_insurance()], [mc_annuity()], [mc_premium()], [mc_reserve()],
#' [summary_mc()]
#'
#' @references
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A.,
#' and Nesbitt, C. J. (1997). *Actuarial Mathematics*. Second Edition.
#' Society of Actuaries.
#'
#' @examples
#' # Example 1: direct use with simulated present values
#' sim_values <- tibble::tibble(
#'   sim_id = 1:5,
#'   pv_benefit = c(0.80, 0.75, 0.95, 0.60, 0.70),
#'   pv_annuity = c(8.0, 7.5, 9.0, 6.0, 7.0),
#'   premium = 0.10
#' )
#'
#' sim_values |>
#'   mc_loss()
#'
#' # Example 2: using a premium supplied directly
#' sim_values |>
#'   mc_loss(premium = 0.12)
#'
#' # Example 3: annual whole life loss
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
#'   mc_premium() |>
#'   mc_loss()
#'
#' # Example 4: monthly whole life loss
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
#'   mc_premium() |>
#'   mc_loss()
#'
#' # Example 5: summarising simulated losses
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
#'   mc_premium() |>
#'   mc_loss() |>
#'   summary_mc(value_col = "loss")
#'
#' # Example 6: joint-life loss with monthly premium annuity
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
#'   mc_premium() |>
#'   mc_loss()
#'
#' @export
mc_loss <- function(data,
                    benefit_col = "pv_benefit",
                    annuity_col = "pv_annuity",
                    premium_col = "premium",
                    loss_col = "loss",
                    premium = NULL) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or tibble.", call. = FALSE)
  }

  .mc_assert_numeric_column(data, benefit_col, "benefit_col")
  .mc_assert_numeric_column(data, annuity_col, "annuity_col")
  .mc_assert_character_scalar(premium_col, "premium_col")
  .mc_assert_character_scalar(loss_col, "loss_col")

  benefit_value <- data[[benefit_col]]
  annuity_value <- data[[annuity_col]]

  if (!is.null(premium)) {
    .mc_assert_numeric_scalar(premium, "premium")
    premium_value <- premium
  } else {
    .mc_assert_numeric_column(data, premium_col, "premium_col")
    premium_value <- data[[premium_col]]
  }

  loss_value <- benefit_value - premium_value * annuity_value

  data |>
    dplyr::mutate(
      "{loss_col}" := loss_value
    )
}
