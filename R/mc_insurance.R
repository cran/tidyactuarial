#' Compute simulated present values for life insurance benefits
#'
#' Computes Monte Carlo simulated present values of life insurance benefits
#' from simulated future lifetimes.
#'
#' This function is designed to be used after [simulate_lifetime()]. It takes
#' simulated values of the curtate future lifetime \eqn{K_x}, and optionally
#' the complete future lifetime \eqn{T_x}, and evaluates the present value
#' random variable associated with classical life insurance benefits.
#'
#' @param data A data frame or tibble containing simulated future lifetimes,
#'   typically returned by [simulate_lifetime()] or by [mc_multilife_status()]
#'   when working with multiple-life statuses.
#' @param rate Numeric scalar. Interest rate used for discounting.
#' @param benefit Numeric scalar. Benefit amount payable under the insurance.
#'   Default is `1`.
#' @param insurance Character string specifying the type of insurance.
#'   Available options are `"whole_life"`, `"term"`, `"deferred"`,
#'   `"deferred_term"`, `"pure_endowment"`, and `"endowment"`.
#' @param term Numeric scalar. Term of the insurance in years. Required for
#'   `"term"`, `"deferred_term"`, `"pure_endowment"`, and `"endowment"`.
#' @param deferral_years Numeric scalar. Deferral period in years. Default is
#'   `0`. Required to be positive for `"deferred"` and `"deferred_term"`
#'   insurance.
#' @param payment_timing Character string specifying when death benefits are
#'   paid. Available options are `"end_of_year"` and `"moment_of_death"`.
#'   Default is `"end_of_year"`.
#' @param interest_type Character string specifying the interest rate convention.
#'   Available options are `"effective"`, `"nominal"`, and `"force"`.
#'   Default is `"effective"`.
#' @param m Numeric scalar. Number of interest conversion periods per year when
#'   `interest_type = "nominal"`. Default is `1`. This argument controls the
#'   interest-rate conversion frequency only. It does not represent benefit
#'   frequency or premium payment frequency.
#' @param k_col Character string. Name of the column containing simulated
#'   curtate future lifetimes. Default is `"Kx"`.
#' @param tx_col Character string. Name of the column containing simulated
#'   complete future lifetimes. Required when
#'   `payment_timing = "moment_of_death"`. Default is `"Tx"`.
#' @param benefit_col Character string. Name of the output column containing
#'   simulated present values of benefits. Default is `"pv_benefit"`.
#'
#' @details
#' Let \eqn{K_x} denote the curtate future lifetime of a life aged \eqn{x},
#' let \eqn{T_x} denote the complete future lifetime, and let \eqn{v} denote
#' the annual discount factor.
#'
#' If `payment_timing = "end_of_year"`, death benefits are discounted using
#' \eqn{K_x + 1}. If `payment_timing = "moment_of_death"`, death benefits are
#' discounted using \eqn{T_x}.
#'
#' The arguments `interest_type` and `m` determine how the supplied rate is
#' converted into an annual effective rate before discounting. When
#' `interest_type = "effective"`, `rate` is interpreted as an annual effective
#' interest rate. When `interest_type = "nominal"`, `rate` is interpreted as a
#' nominal annual interest rate convertible `m` times per year. When
#' `interest_type = "force"`, `rate` is interpreted as a constant force of
#' interest.
#'
#' The following insurance types are supported:
#'
#' * `"whole_life"`: benefit is paid whenever death occurs.
#'
#' * `"term"`: benefit is paid if death occurs within `term` years.
#'
#' * `"deferred"`: benefit is paid if death occurs after the deferral period.
#'
#' * `"deferred_term"`: benefit is paid if death occurs after the deferral
#'   period and within the following `term` years.
#'
#' * `"pure_endowment"`: benefit is paid at time `term` if the life survives
#'   to that time.
#'
#' * `"endowment"`: death benefit is paid if death occurs within `term` years;
#'   otherwise, a survival benefit is paid at time `term`.
#'
#' The function returns simulated present values, not only their expected
#' values. Therefore the resulting column can be summarized with
#' [summary_mc()], plotted with `ggplot2`, or used to construct premiums,
#' losses, and reserves.
#'
#' @return A tibble with the original simulation columns and additional columns:
#'
#' * `rate`: original rate supplied.
#' * `interest_type`: interest rate convention.
#' * `m`: interest conversion frequency.
#' * `effective_rate`: equivalent annual effective interest rate.
#' * `discount_factor`: annual discount factor.
#' * `insurance`: insurance type.
#' * `benefit`: benefit amount.
#' * `term`: insurance term, if applicable.
#' * `deferral_years`: deferral period.
#' * `payment_timing`: timing used for death benefits.
#' * `benefit_time`: simulated payment time of the benefit.
#' * `benefit_indicator`: indicator that the benefit is paid.
#' * `pv_benefit`: simulated present value of the benefit, or another name
#'   supplied through `benefit_col`.
#'
#' @seealso
#' [simulate_lifetime()], [simulate_lifetimes()], [mc_multilife_status()],
#' [mc_annuity()], [mc_premium()], [mc_loss()], [mc_reserve()],
#' [summary_mc()]
#'
#' @references
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A.,
#' and Nesbitt, C. J. (1997). *Actuarial Mathematics*. Second Edition.
#' Society of Actuaries.
#'
#' @examples
#' life_table <- tibble::tibble(
#'   age = 40:100,
#'   qx = seq(0.002, 1, length.out = 61)
#' )
#'
#' # Whole life insurance payable at the end of the year of death
#' life_table |>
#'   simulate_lifetime(age = 40, n_sim = 1000, seed = 123) |>
#'   mc_insurance(
#'     rate = 0.05,
#'     insurance = "whole_life",
#'     benefit = 1
#'   )
#'
#' # 20-year term insurance
#' life_table |>
#'   simulate_lifetime(age = 40, n_sim = 1000, seed = 123) |>
#'   mc_insurance(
#'     rate = 0.05,
#'     insurance = "term",
#'     term = 20,
#'     benefit = 100000
#'   )
#'
#' # 10-year deferred whole life insurance
#' life_table |>
#'   simulate_lifetime(age = 40, n_sim = 1000, seed = 123) |>
#'   mc_insurance(
#'     rate = 0.05,
#'     insurance = "deferred",
#'     deferral_years = 10,
#'     benefit = 1
#'   )
#'
#' # 10-year deferred, 20-year term insurance
#' life_table |>
#'   simulate_lifetime(age = 40, n_sim = 1000, seed = 123) |>
#'   mc_insurance(
#'     rate = 0.05,
#'     insurance = "deferred_term",
#'     deferral_years = 10,
#'     term = 20,
#'     benefit = 1
#'   )
#'
#' # Pure endowment
#' life_table |>
#'   simulate_lifetime(age = 40, n_sim = 1000, seed = 123) |>
#'   mc_insurance(
#'     rate = 0.05,
#'     insurance = "pure_endowment",
#'     term = 20,
#'     benefit = 1
#'   )
#'
#' # Endowment insurance payable at the moment of death if death occurs
#' life_table |>
#'   simulate_lifetime(age = 40, n_sim = 1000, fractional = "udd", seed = 123) |>
#'   mc_insurance(
#'     rate = 0.05,
#'     insurance = "endowment",
#'     term = 20,
#'     payment_timing = "moment_of_death",
#'     benefit = 1
#'   )
#'
#' # Nominal rate convertible monthly
#' life_table |>
#'   simulate_lifetime(age = 40, n_sim = 1000, seed = 123) |>
#'   mc_insurance(
#'     rate = 0.06,
#'     interest_type = "nominal",
#'     m = 12,
#'     insurance = "whole_life",
#'     benefit = 1
#'   )
#'
#' # First-death insurance using a multiple-life status
#' life_table |>
#'   simulate_lifetimes(
#'     ages = c(60, 58),
#'     n_sim = 1000,
#'     fractional = "udd",
#'     seed = 123
#'   ) |>
#'   mc_multilife_status(status = "first_death") |>
#'   mc_insurance(
#'     rate = 0.04,
#'     insurance = "whole_life",
#'     benefit = 100000,
#'     k_col = "K_status",
#'     tx_col = "T_status"
#'   )
#'
#' @export
mc_insurance <- function(data,
                         rate,
                         benefit = 1,
                         insurance = c(
                           "whole_life",
                           "term",
                           "deferred",
                           "deferred_term",
                           "pure_endowment",
                           "endowment"
                         ),
                         term = NULL,
                         deferral_years = 0,
                         payment_timing = c("end_of_year", "moment_of_death"),
                         interest_type = c("effective", "nominal", "force"),
                         m = 1,
                         k_col = "Kx",
                         tx_col = "Tx",
                         benefit_col = "pv_benefit") {
  insurance <- match.arg(insurance)
  payment_timing <- match.arg(payment_timing)
  interest_type <- match.arg(interest_type)

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or tibble.", call. = FALSE)
  }

  .mc_assert_numeric_scalar(rate, "rate")
  .mc_assert_numeric_scalar(benefit, "benefit", min = 0)
  .mc_assert_numeric_scalar(m, "m", min = 0, strict_min = TRUE)
  .mc_assert_numeric_column(data, k_col, "k_col")
  .mc_assert_character_scalar(benefit_col, "benefit_col")

  if (payment_timing == "moment_of_death") {
    .mc_assert_numeric_column(data, tx_col, "tx_col")
  }

  if (insurance %in% c("term", "deferred_term", "pure_endowment", "endowment")) {
    .mc_assert_numeric_scalar(term, "term", min = 0, strict_min = TRUE)
  }

  if (insurance %in% c("deferred", "deferred_term")) {
    .mc_assert_numeric_scalar(
      deferral_years,
      "deferral_years",
      min = 0,
      strict_min = TRUE
    )
  } else {
    .mc_assert_numeric_scalar(deferral_years, "deferral_years", min = 0)
  }

  effective_rate <- .mc_effective_rate(
    rate = rate,
    interest_type = interest_type,
    m = m
  )

  v <- 1 / (1 + effective_rate)

  Kx <- data[[k_col]]

  if (any(Kx < 0, na.rm = TRUE)) {
    stop("`k_col` must contain non-negative simulated lifetimes.", call. = FALSE)
  }

  Tx <- NULL

  if (payment_timing == "moment_of_death") {
    Tx <- data[[tx_col]]

    if (any(Tx < 0, na.rm = TRUE)) {
      stop("`tx_col` must contain non-negative simulated lifetimes.", call. = FALSE)
    }

    if (all(is.na(Tx))) {
      stop(
        "`tx_col` contains only missing values. Benefits payable at the moment ",
        "of death require complete future lifetimes.",
        call. = FALSE
      )
    }
  }

  death_time <- if (payment_timing == "end_of_year") {
    Kx + 1
  } else {
    Tx
  }

  benefit_indicator <- rep(0, length(Kx))
  benefit_time <- rep(NA_real_, length(Kx))

  if (insurance == "whole_life") {
    benefit_indicator <- rep(1, length(Kx))
    benefit_time <- death_time
  }

  if (insurance == "term") {
    if (payment_timing == "end_of_year") {
      benefit_indicator <- as.numeric(Kx < term)
    } else {
      benefit_indicator <- as.numeric(Tx <= term)
    }

    benefit_time <- ifelse(benefit_indicator == 1, death_time, NA_real_)
  }

  if (insurance == "deferred") {
    if (payment_timing == "end_of_year") {
      benefit_indicator <- as.numeric(Kx >= deferral_years)
    } else {
      benefit_indicator <- as.numeric(Tx > deferral_years)
    }

    benefit_time <- ifelse(benefit_indicator == 1, death_time, NA_real_)
  }

  if (insurance == "deferred_term") {
    if (payment_timing == "end_of_year") {
      benefit_indicator <- as.numeric(
        Kx >= deferral_years & Kx < deferral_years + term
      )
    } else {
      benefit_indicator <- as.numeric(
        Tx > deferral_years & Tx <= deferral_years + term
      )
    }

    benefit_time <- ifelse(benefit_indicator == 1, death_time, NA_real_)
  }

  if (insurance == "pure_endowment") {
    if (payment_timing == "end_of_year") {
      benefit_indicator <- as.numeric(Kx >= term)
    } else {
      benefit_indicator <- as.numeric(Tx >= term)
    }

    benefit_time <- ifelse(benefit_indicator == 1, term, NA_real_)
  }

  if (insurance == "endowment") {
    if (payment_timing == "end_of_year") {
      death_before_term <- Kx < term
    } else {
      death_before_term <- Tx <= term
    }

    benefit_indicator <- rep(1, length(Kx))
    benefit_time <- ifelse(death_before_term, death_time, term)
  }

  pv_benefit <- ifelse(
    benefit_indicator == 1,
    benefit * v^benefit_time,
    0
  )

  data |>
    dplyr::mutate(
      rate = rate,
      interest_type = interest_type,
      m = m,
      effective_rate = effective_rate,
      discount_factor = v,
      insurance = insurance,
      benefit = benefit,
      term = term,
      deferral_years = deferral_years,
      payment_timing = payment_timing,
      benefit_time = benefit_time,
      benefit_indicator = benefit_indicator,
      "{benefit_col}" := pv_benefit
    )
}
