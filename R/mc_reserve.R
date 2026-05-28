#' Compute Monte Carlo prospective reserves for life contingencies
#'
#' Computes simulated prospective reserve losses at one or more policy
#' durations from simulated future lifetimes.
#'
#' This function recalculates future benefit and future premium present values
#' from each valuation duration. It is not a wrapper around [mc_loss()], because
#' reserves require valuing only the cash flows that remain after the valuation
#' time.
#'
#' For a policy in force at duration \eqn{t}, the simulated prospective loss is
#'
#' \deqn{
#'   L_t = Z_t - P Y_t,
#' }
#'
#' where \eqn{Z_t} is the present value at duration \eqn{t} of future benefits,
#' \eqn{Y_t} is the present value at duration \eqn{t} of future premium
#' payments, and \eqn{P} is the premium.
#'
#' @param data A data frame or tibble containing simulated future lifetimes,
#'   typically returned by [simulate_lifetime()] or [mc_multilife_status()].
#' @param duration Numeric vector. Policy duration or durations at which the
#'   reserve is evaluated. Default is `0`.
#' @param rate Numeric scalar. Interest rate used for discounting.
#' @param premium Optional numeric scalar. Premium used in the reserve loss.
#'   If `NULL`, the function first tries to use `premium_col` from `data`.
#'   If `premium_col` is not found, a Monte Carlo net premium at issue is
#'   estimated internally as \eqn{\bar{Z}_0 / \bar{Y}_0}.
#' @param premium_col Character string. Name of the premium column in `data`.
#'   Default is `"premium"`.
#' @param benefit Numeric scalar. Benefit amount payable under the insurance.
#'   Default is `1`.
#' @param payment Numeric scalar. Amount of each premium annuity payment.
#'   Default is `1`.
#' @param payments_per_year Positive integer-like scalar. Number of premium
#'   payments per year. Default is `1`, corresponding to annual premiums. Use
#'   `12` for monthly premiums, `4` for quarterly premiums, and `2` for
#'   semiannual premiums.
#' @param insurance Character string specifying the insurance type. Available
#'   options are `"whole_life"`, `"term"`, `"deferred"`,
#'   `"deferred_term"`, `"pure_endowment"`, and `"endowment"`.
#' @param annuity Character string specifying the premium annuity type.
#'   Available options are `"whole_life"`, `"temporary"`, `"deferred"`,
#'   `"deferred_temporary"`, `"certain"`, and `"guaranteed"`.
#' @param term Numeric scalar. Contract term in years. Required for insurance
#'   types `"term"`, `"deferred_term"`, `"pure_endowment"`, and `"endowment"`,
#'   and for annuity types `"temporary"`, `"deferred_temporary"`, and
#'   `"certain"`.
#' @param deferral_years Numeric scalar. Deferral period in years. Default is
#'   `0`.
#' @param guarantee_years Numeric scalar. Guaranteed payment period in years.
#'   Required when `annuity = "guaranteed"`.
#' @param payment_timing Character string specifying when death benefits are
#'   paid. Available options are `"end_of_year"` and `"moment_of_death"`.
#'   Default is `"end_of_year"`.
#' @param premium_timing Character string specifying the premium payment timing.
#'   Available options are `"due"` and `"immediate"`. Default is `"due"`.
#' @param reserve_timing Character string specifying whether payments due
#'   exactly at the valuation duration are included. Available options are
#'   `"before_payment"` and `"after_payment"`. Default is `"before_payment"`.
#' @param interest_type Character string specifying the interest rate convention.
#'   Available options are `"effective"`, `"nominal"`, and `"force"`.
#'   Default is `"effective"`.
#' @param m Numeric scalar. Number of interest conversion periods per year when
#'   `interest_type = "nominal"`. Default is `1`. This argument controls the
#'   interest-rate conversion frequency only. It does not represent premium
#'   payment frequency.
#' @param k_col Character string. Name of the column containing simulated
#'   curtate future lifetimes. Default is `"Kx"`.
#' @param tx_col Character string. Name of the column containing simulated
#'   complete future lifetimes. Default is `"Tx"`.
#' @param in_force_basis Character string specifying how the in-force indicator
#'   is evaluated. Available options are `"auto"`, `"complete"`, and
#'   `"curtate"`. Default is `"auto"`.
#' @param not_in_force Character string specifying what to return for scenarios
#'   that are not in force at the valuation duration. Available options are
#'   `"na"` and `"zero"`. Default is `"na"`.
#' @param reserve_col Character string. Name of the output reserve loss column.
#'   Default is `"reserve_loss"`.
#'
#' @details
#' The arguments `m` and `payments_per_year` have different meanings:
#'
#' * `m` is used only when `interest_type = "nominal"` and controls the
#'   frequency of interest conversion.
#' * `payments_per_year` controls how frequently future premium payments are
#'   made.
#'
#' The argument `payment` represents the amount of each premium payment. Thus,
#' for monthly premiums with total annual premium equal to 1, use
#' `payment = 1 / 12` and `payments_per_year = 12`.
#'
#' If `payments_per_year = 1`, future premiums are annual. If
#' `payments_per_year > 1`, future premiums are made at fractional times, and a
#' valid complete future lifetime column supplied through `tx_col` is required.
#'
#' Durations may be integer or fractional. Fractional reserve durations require
#' complete future lifetimes. For example, monthly reserve calculations may use
#' `duration = seq(0, 20, by = 1 / 12)`.
#'
#' If `reserve_timing = "before_payment"`, cash flows occurring exactly at the
#' valuation duration are included. If `reserve_timing = "after_payment"`, cash
#' flows occurring exactly at the valuation duration are excluded.
#'
#' If `not_in_force = "na"`, scenarios that are not in force at a given duration
#' receive `NA` values for future present values and reserve losses. This is
#' useful for estimating reserves conditional on the policy still being in
#' force. If `not_in_force = "zero"`, those scenarios receive zero values,
#' which may be useful for portfolio run-off summaries.
#'
#' This function computes prospective reserves under the simulated model. It
#' does not include expenses, surrender values, taxes, profit loadings, or
#' statutory reserving adjustments.
#'
#' @return A tibble with one row per original simulation and per requested
#' duration. It contains the original simulation columns and additional columns:
#'
#' * `duration`: valuation duration.
#' * `in_force`: logical indicator for survival or in-force status.
#' * `reserve_timing`: timing convention used at valuation.
#' * `not_in_force`: convention used for scenarios not in force.
#' * `rate`: original rate supplied.
#' * `interest_type`: interest rate convention.
#' * `m`: interest conversion frequency.
#' * `effective_rate`: equivalent annual effective interest rate.
#' * `discount_factor`: annual discount factor.
#' * `insurance`: insurance type.
#' * `annuity`: premium annuity type.
#' * `benefit`: benefit amount.
#' * `payment`: amount of each premium payment.
#' * `payments_per_year`: premium payment frequency.
#' * `future_pv_benefit`: present value at duration of future benefits.
#' * `future_pv_premiums`: present value at duration of future premium payments.
#' * `premium`: premium used in the simulated reserve loss.
#' * `reserve_loss`: simulated reserve loss, or another name supplied through
#'   `reserve_col`.
#'
#' @seealso
#' [simulate_lifetime()], [simulate_lifetimes()], [mc_multilife_status()],
#' [mc_insurance()], [mc_annuity()], [mc_premium()], [mc_loss()],
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
#' # Annual prospective reserves for whole life insurance
#' life_table |>
#'   simulate_lifetime(age = 40, n_sim = 1000, seed = 123) |>
#'   mc_reserve(
#'     duration = c(0, 5, 10, 20),
#'     rate = 0.05,
#'     insurance = "whole_life",
#'     annuity = "whole_life",
#'     benefit = 1,
#'     payment = 1,
#'     payments_per_year = 1,
#'     premium_timing = "due"
#'   )
#'
#' # Monthly premium reserve curve
#' life_table |>
#'   simulate_lifetime(
#'     age = 40,
#'     n_sim = 1000,
#'     fractional = "udd",
#'     seed = 123
#'   ) |>
#'   mc_reserve(
#'     duration = seq(0, 10, by = 1),
#'     rate = 0.05,
#'     insurance = "whole_life",
#'     annuity = "whole_life",
#'     benefit = 1,
#'     payment = 1 / 12,
#'     payments_per_year = 12,
#'     premium_timing = "due"
#'   )
#'
#' # Monthly reserves by fractional policy duration
#' life_table |>
#'   simulate_lifetime(
#'     age = 40,
#'     n_sim = 1000,
#'     fractional = "udd",
#'     seed = 123
#'   ) |>
#'   mc_reserve(
#'     duration = seq(0, 5, by = 1 / 12),
#'     rate = 0.05,
#'     insurance = "whole_life",
#'     annuity = "whole_life",
#'     benefit = 1,
#'     payment = 1 / 12,
#'     payments_per_year = 12,
#'     premium_timing = "due"
#'   ) |>
#'   summary_mc(value_col = "reserve_loss", by = "duration")
#'
#' # Joint-life reserve using a multiple-life status
#' life_table |>
#'   simulate_lifetimes(
#'     ages = c(60, 58),
#'     n_sim = 1000,
#'     fractional = "udd",
#'     seed = 123
#'   ) |>
#'   mc_multilife_status(status = "joint_life") |>
#'   mc_reserve(
#'     duration = c(0, 5, 10),
#'     rate = 0.04,
#'     insurance = "whole_life",
#'     annuity = "whole_life",
#'     benefit = 1,
#'     payment = 1,
#'     payments_per_year = 1,
#'     k_col = "K_status",
#'     tx_col = "T_status"
#'   )
#'
#' @export
mc_reserve <- function(data,
                       duration = 0,
                       rate,
                       premium = NULL,
                       premium_col = "premium",
                       benefit = 1,
                       payment = 1,
                       payments_per_year = 1,
                       insurance = c(
                         "whole_life",
                         "term",
                         "deferred",
                         "deferred_term",
                         "pure_endowment",
                         "endowment"
                       ),
                       annuity = c(
                         "whole_life",
                         "temporary",
                         "deferred",
                         "deferred_temporary",
                         "certain",
                         "guaranteed"
                       ),
                       term = NULL,
                       deferral_years = 0,
                       guarantee_years = NULL,
                       payment_timing = c("end_of_year", "moment_of_death"),
                       premium_timing = c("due", "immediate"),
                       reserve_timing = c("before_payment", "after_payment"),
                       interest_type = c("effective", "nominal", "force"),
                       m = 1,
                       k_col = "Kx",
                       tx_col = "Tx",
                       in_force_basis = c("auto", "complete", "curtate"),
                       not_in_force = c("na", "zero"),
                       reserve_col = "reserve_loss") {
  insurance <- match.arg(insurance)
  annuity <- match.arg(annuity)
  payment_timing <- match.arg(payment_timing)
  premium_timing <- match.arg(premium_timing)
  reserve_timing <- match.arg(reserve_timing)
  interest_type <- match.arg(interest_type)
  in_force_basis <- match.arg(in_force_basis)
  not_in_force <- match.arg(not_in_force)

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or tibble.", call. = FALSE)
  }

  if (!is.numeric(duration) ||
      length(duration) == 0 ||
      anyNA(duration) ||
      any(!is.finite(duration)) ||
      any(duration < 0)) {
    stop("`duration` must be a non-negative numeric vector.", call. = FALSE)
  }

  .mc_assert_numeric_scalar(rate, "rate")
  .mc_assert_numeric_scalar(benefit, "benefit", min = 0)
  .mc_assert_numeric_scalar(payment, "payment", min = 0)
  .mc_assert_positive_integer(payments_per_year, "payments_per_year")
  .mc_assert_numeric_scalar(m, "m", min = 0, strict_min = TRUE)
  .mc_assert_numeric_column(data, k_col, "k_col")
  .mc_assert_character_scalar(premium_col, "premium_col")
  .mc_assert_character_scalar(reserve_col, "reserve_col")
  .mc_assert_character_scalar(tx_col, "tx_col")

  if (insurance %in% c("term", "deferred_term", "pure_endowment", "endowment") ||
      annuity %in% c("temporary", "deferred_temporary", "certain")) {
    .mc_assert_numeric_scalar(term, "term", min = 0, strict_min = TRUE)
  }

  if (insurance %in% c("deferred", "deferred_term") ||
      annuity %in% c("deferred", "deferred_temporary")) {
    .mc_assert_numeric_scalar(
      deferral_years,
      "deferral_years",
      min = 0,
      strict_min = TRUE
    )
  } else {
    .mc_assert_numeric_scalar(deferral_years, "deferral_years", min = 0)
  }

  if (annuity == "guaranteed") {
    .mc_assert_numeric_scalar(
      guarantee_years,
      "guarantee_years",
      min = 0,
      strict_min = TRUE
    )
  } else {
    .mc_assert_numeric_scalar(
      guarantee_years,
      "guarantee_years",
      min = 0,
      allow_null = TRUE
    )
  }

  if (!is.null(premium)) {
    .mc_assert_numeric_scalar(premium, "premium")
  }

  Kx <- data[[k_col]]

  if (any(Kx < 0, na.rm = TRUE)) {
    stop("`k_col` must contain non-negative simulated lifetimes.", call. = FALSE)
  }

  has_tx <- tx_col %in% names(data) &&
    is.numeric(data[[tx_col]]) &&
    !all(is.na(data[[tx_col]]))

  if (payment_timing == "moment_of_death" && !has_tx) {
    stop(
      "Benefits payable at the moment of death require a valid complete ",
      "future lifetime column identified by `tx_col`.",
      call. = FALSE
    )
  }

  if (payments_per_year > 1 && annuity != "certain" && !has_tx) {
    stop(
      "Fractional life-contingent premium payments require a valid complete ",
      "future lifetime column identified by `tx_col`.",
      call. = FALSE
    )
  }

  has_fractional_duration <- any(
    abs(duration - round(duration)) > sqrt(.Machine$double.eps)
  )

  if (has_fractional_duration && !has_tx) {
    stop(
      "Fractional reserve durations require a valid complete future lifetime ",
      "column identified by `tx_col`.",
      call. = FALSE
    )
  }

  if (in_force_basis == "auto") {
    in_force_basis <- if (has_tx) "complete" else "curtate"
  }

  if (in_force_basis == "complete") {
    .mc_assert_numeric_column(data, tx_col, "tx_col")
  }

  if (in_force_basis == "curtate" && has_fractional_duration) {
    stop(
      "Fractional reserve durations require `in_force_basis = 'complete'` ",
      "and a valid `tx_col`.",
      call. = FALSE
    )
  }

  effective_rate <- .mc_effective_rate(
    rate = rate,
    interest_type = interest_type,
    m = m
  )

  v <- 1 / (1 + effective_rate)

  Tx <- if (has_tx) {
    data[[tx_col]]
  } else {
    rep(NA_real_, length(Kx))
  }

  death_time <- if (payment_timing == "end_of_year") {
    Kx + 1
  } else {
    Tx
  }

  benefit_info <- .mc_reserve_benefit_info(
    Kx = Kx,
    Tx = Tx,
    death_time = death_time,
    insurance = insurance,
    term = term,
    deferral_years = deferral_years,
    payment_timing = payment_timing
  )

  premium_payment_times <- Map(
    function(k, tx) {
      if (payments_per_year == 1) {
        return(
          .mc_annuity_payment_times_annual(
            k = k,
            annuity = annuity,
            term = term,
            deferral_years = deferral_years,
            guarantee_years = guarantee_years,
            timing = premium_timing
          )
        )
      }

      .mc_annuity_payment_times_fractional(
        tx = tx,
        annuity = annuity,
        term = term,
        deferral_years = deferral_years,
        guarantee_years = guarantee_years,
        timing = premium_timing,
        payment_interval = 1 / payments_per_year
      )
    },
    Kx,
    Tx
  )

  issue_pv_benefit <- ifelse(
    benefit_info$benefit_indicator == 1,
    benefit * v^benefit_info$benefit_time,
    0
  )

  issue_pv_premiums <- vapply(
    premium_payment_times,
    function(times) {
      payment * sum(v^times)
    },
    numeric(1)
  )

  if (!is.null(premium)) {
    premium_used <- rep(premium, length(Kx))
  } else if (premium_col %in% names(data)) {
    if (!is.numeric(data[[premium_col]])) {
      stop("`premium_col` must identify a numeric column.", call. = FALSE)
    }

    premium_used <- data[[premium_col]]
  } else {
    denominator <- mean(issue_pv_premiums)

    if (is.na(denominator) || denominator <= 0) {
      stop(
        "The mean simulated premium annuity present value must be positive ",
        "when estimating the net premium internally.",
        call. = FALSE
      )
    }

    premium_hat <- mean(issue_pv_benefit) / denominator
    premium_used <- rep(premium_hat, length(Kx))
  }

  row_id <- seq_len(nrow(data))
  n_duration <- length(duration)

  expanded <- tibble::as_tibble(data)[
    rep(row_id, each = n_duration),
    ,
    drop = FALSE
  ]

  row_rep <- rep(row_id, each = n_duration)
  duration_rep <- rep(duration, times = length(row_id))

  Kx_rep <- Kx[row_rep]
  Tx_rep <- Tx[row_rep]
  premium_rep <- premium_used[row_rep]

  in_force <- if (in_force_basis == "complete") {
    Tx_rep > duration_rep
  } else {
    Kx_rep >= duration_rep
  }

  benefit_time_rep <- benefit_info$benefit_time[row_rep]
  benefit_indicator_rep <- benefit_info$benefit_indicator[row_rep]

  include_benefit <- if (reserve_timing == "before_payment") {
    benefit_time_rep >= duration_rep
  } else {
    benefit_time_rep > duration_rep
  }

  future_pv_benefit <- ifelse(
    benefit_indicator_rep == 1 & include_benefit,
    benefit * v^(benefit_time_rep - duration_rep),
    0
  )

  future_pv_premiums <- vapply(
    seq_along(row_rep),
    function(i) {
      times <- premium_payment_times[[row_rep[i]]]

      future_times <- if (reserve_timing == "before_payment") {
        times[times >= duration_rep[i]]
      } else {
        times[times > duration_rep[i]]
      }

      payment * sum(v^(future_times - duration_rep[i]))
    },
    numeric(1)
  )

  if (not_in_force == "na") {
    future_pv_benefit[!in_force] <- NA_real_
    future_pv_premiums[!in_force] <- NA_real_
  } else {
    future_pv_benefit[!in_force] <- 0
    future_pv_premiums[!in_force] <- 0
  }

  reserve_loss <- future_pv_benefit - premium_rep * future_pv_premiums

  expanded |>
    dplyr::mutate(
      duration = duration_rep,
      in_force = in_force,
      reserve_timing = reserve_timing,
      not_in_force = not_in_force,
      rate = rate,
      interest_type = interest_type,
      m = m,
      effective_rate = effective_rate,
      discount_factor = v,
      insurance = insurance,
      annuity = annuity,
      benefit = benefit,
      payment = payment,
      payments_per_year = payments_per_year,
      term = term,
      deferral_years = deferral_years,
      guarantee_years = guarantee_years,
      payment_timing = payment_timing,
      premium_timing = premium_timing,
      future_pv_benefit = future_pv_benefit,
      future_pv_premiums = future_pv_premiums,
      premium = premium_rep,
      "{reserve_col}" := reserve_loss
    )
}


#' Internal helper: benefit payment information for reserve calculations
#'
#' Determines whether a benefit is paid and at what time for each simulated
#' scenario.
#'
#' @param Kx Numeric vector of simulated curtate future lifetimes.
#' @param Tx Numeric vector of simulated complete future lifetimes.
#' @param death_time Numeric vector of death benefit payment times.
#' @param insurance Character string specifying the insurance type.
#' @param term Numeric scalar or `NULL`.
#' @param deferral_years Numeric scalar.
#' @param payment_timing Character string. Either `"end_of_year"` or
#'   `"moment_of_death"`.
#'
#' @return A list with `benefit_indicator` and `benefit_time`.
#'
#' @keywords internal
.mc_reserve_benefit_info <- function(Kx,
                                     Tx,
                                     death_time,
                                     insurance,
                                     term,
                                     deferral_years,
                                     payment_timing) {
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

  list(
    benefit_indicator = benefit_indicator,
    benefit_time = benefit_time
  )
}
