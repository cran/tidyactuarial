#' Compute simulated present values for life annuities
#'
#' Computes Monte Carlo simulated present values of life annuity payments from
#' simulated future lifetimes.
#'
#' This function is designed to be used after [simulate_lifetime()]. It takes
#' simulated values of the curtate future lifetime \eqn{K_x}, and when needed
#' the complete future lifetime \eqn{T_x}, and evaluates the present value
#' random variable associated with several classical annuity benefits.
#'
#' @param data A data frame or tibble containing simulated future lifetimes,
#'   typically returned by [simulate_lifetime()] or by [mc_multilife_status()]
#'   when working with multiple-life statuses.
#' @param rate Numeric scalar. Interest rate used for discounting.
#' @param payment Numeric scalar. Amount of each annuity payment. Default is
#'   `1`.
#' @param payments_per_year Positive integer-like scalar. Number of annuity
#'   payments per year. Default is `1`, corresponding to annual payments. For
#'   example, use `payments_per_year = 12` for monthly payments, `4` for
#'   quarterly payments, and `2` for semiannual payments.
#' @param annuity Character string specifying the annuity type. Available
#'   options are `"whole_life"`, `"temporary"`, `"deferred"`,
#'   `"deferred_temporary"`, `"certain"`, and `"guaranteed"`.
#' @param term Numeric scalar. Term of the annuity in years. Required for
#'   `"temporary"`, `"deferred_temporary"`, and `"certain"` annuities.
#' @param deferral_years Numeric scalar. Deferral period in years. Default is
#'   `0`. Required to be positive for `"deferred"` and
#'   `"deferred_temporary"` annuities.
#' @param guarantee_years Numeric scalar. Guaranteed payment period in years.
#'   Required for `"guaranteed"` annuities.
#' @param timing Character string specifying the annuity payment timing.
#'   Available options are `"immediate"` and `"due"`. Default is `"immediate"`.
#' @param interest_type Character string specifying the interest rate convention.
#'   Available options are `"effective"`, `"nominal"`, and `"force"`.
#'   Default is `"effective"`.
#' @param m Numeric scalar. Number of interest conversion periods per year when
#'   `interest_type = "nominal"`. Default is `1`. This argument controls the
#'   interest-rate conversion frequency only. It does not represent annuity
#'   payment frequency.
#' @param k_col Character string. Name of the column containing simulated
#'   curtate future lifetimes. Default is `"Kx"`.
#' @param tx_col Character string. Name of the column containing simulated
#'   complete future lifetimes. Default is `"Tx"`. This column is required when
#'   `payments_per_year > 1`, except for annuities certain.
#' @param annuity_col Character string. Name of the output column containing
#'   simulated present values of annuity payments. Default is `"pv_annuity"`.
#'
#' @details
#' Let \eqn{K_x} denote the curtate future lifetime of a life aged \eqn{x},
#' let \eqn{T_x} denote the complete future lifetime, and let \eqn{v} denote
#' the annual discount factor.
#'
#' The arguments `m` and `payments_per_year` have different meanings:
#'
#' * `m` is used only when `interest_type = "nominal"` and controls the
#'   frequency of interest conversion.
#' * `payments_per_year` controls how frequently annuity payments are made.
#'
#' The argument `payment` represents the amount of each annuity payment. Thus,
#' for a monthly annuity with total annual payment equal to 1, use
#' `payment = 1 / 12` and `payments_per_year = 12`.
#'
#' For annual payments, `payments_per_year = 1`, the function works directly
#' with \eqn{K_x}. For fractional payments, such as monthly, quarterly, or
#' semiannual payments, the function uses \eqn{T_x} to determine whether the
#' life is alive at each fractional payment time.
#'
#' For annual whole life annuity-immediate, the simulated present value is
#'
#' \deqn{
#'   Y = \sum_{j=1}^{K_x} c v^j,
#' }
#'
#' where \eqn{c} is the amount of each payment.
#'
#' For annual whole life annuity-due, the simulated present value is
#'
#' \deqn{
#'   \ddot{Y} = \sum_{j=0}^{K_x} c v^j.
#' }
#'
#' For fractional payments with payment frequency \eqn{m_p}, annuity-immediate
#' payments are made at times \eqn{1/m_p, 2/m_p, \ldots} while the life is
#' alive. Annuity-due payments are made at times
#' \eqn{0, 1/m_p, 2/m_p, \ldots} while the life is alive.
#'
#' The following annuity types are supported:
#'
#' * `"whole_life"`: payments continue while the life is alive.
#' * `"temporary"`: payments continue while the life is alive, but for at most
#'   `term` years.
#' * `"deferred"`: payments begin after `deferral_years` years and continue
#'   while the life is alive.
#' * `"deferred_temporary"`: payments begin after `deferral_years` years and
#'   continue while the life is alive, but for at most `term` years after the
#'   deferral period.
#' * `"certain"`: payments are made for `term` years regardless of survival.
#' * `"guaranteed"`: payments continue while the life is alive, with at least
#'   `guarantee_years` years of payments guaranteed.
#'
#' The function returns simulated present values, not only their expected value.
#' Therefore the resulting column can be summarized with [summary_mc()],
#' plotted with `ggplot2`, or used to construct premiums, losses, and reserves.
#'
#' @return A tibble with the original simulation columns and additional columns:
#'
#' * `rate`: original rate supplied.
#' * `interest_type`: interest rate convention.
#' * `m`: interest conversion frequency.
#' * `effective_rate`: equivalent annual effective interest rate.
#' * `discount_factor`: annual discount factor.
#' * `annuity`: annuity type.
#' * `payment`: amount of each annuity payment.
#' * `payments_per_year`: annuity payment frequency.
#' * `term`: annuity term, if applicable.
#' * `deferral_years`: deferral period.
#' * `guarantee_years`: guaranteed period, if applicable.
#' * `timing`: annuity payment timing.
#' * `n_payments`: number of payments made in the simulated scenario.
#' * `first_payment_time`: first payment time in the simulated scenario.
#' * `last_payment_time`: last payment time in the simulated scenario.
#' * `pv_annuity`: simulated present value of annuity payments, or another
#'   name supplied through `annuity_col`.
#'
#' @seealso
#' [simulate_lifetime()], [simulate_lifetimes()], [mc_multilife_status()],
#' [mc_insurance()], [mc_premium()], [mc_loss()], [mc_reserve()],
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
#' # Annual whole life annuity-due
#' life_table |>
#'   simulate_lifetime(age = 40, n_sim = 1000, seed = 123) |>
#'   mc_annuity(
#'     rate = 0.05,
#'     annuity = "whole_life",
#'     payment = 1,
#'     payments_per_year = 1,
#'     timing = "due"
#'   )
#'
#' # Monthly whole life annuity-due with total annual payment equal to 1
#' life_table |>
#'   simulate_lifetime(
#'     age = 40,
#'     n_sim = 1000,
#'     fractional = "udd",
#'     seed = 123
#'   ) |>
#'   mc_annuity(
#'     rate = 0.05,
#'     annuity = "whole_life",
#'     payment = 1 / 12,
#'     payments_per_year = 12,
#'     timing = "due"
#'   )
#'
#' # Quarterly temporary life annuity-immediate
#' life_table |>
#'   simulate_lifetime(
#'     age = 40,
#'     n_sim = 1000,
#'     fractional = "udd",
#'     seed = 123
#'   ) |>
#'   mc_annuity(
#'     rate = 0.05,
#'     annuity = "temporary",
#'     term = 20,
#'     payment = 1 / 4,
#'     payments_per_year = 4,
#'     timing = "immediate"
#'   )
#'
#' # Monthly joint-life annuity using a multiple-life status
#' life_table |>
#'   simulate_lifetimes(
#'     ages = c(60, 58),
#'     n_sim = 1000,
#'     fractional = "udd",
#'     seed = 123
#'   ) |>
#'   mc_multilife_status(status = "joint_life") |>
#'   mc_annuity(
#'     rate = 0.04,
#'     annuity = "whole_life",
#'     payment = 1 / 12,
#'     payments_per_year = 12,
#'     timing = "due",
#'     k_col = "K_status",
#'     tx_col = "T_status"
#'   )
#'
#' # Nominal rate convertible monthly, with quarterly payments
#' life_table |>
#'   simulate_lifetime(
#'     age = 40,
#'     n_sim = 1000,
#'     fractional = "udd",
#'     seed = 123
#'   ) |>
#'   mc_annuity(
#'     rate = 0.06,
#'     interest_type = "nominal",
#'     m = 12,
#'     annuity = "whole_life",
#'     payment = 1 / 4,
#'     payments_per_year = 4,
#'     timing = "due"
#'   )
#'
#' @export
mc_annuity <- function(data,
                       rate,
                       payment = 1,
                       payments_per_year = 1,
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
                       timing = c("immediate", "due"),
                       interest_type = c("effective", "nominal", "force"),
                       m = 1,
                       k_col = "Kx",
                       tx_col = "Tx",
                       annuity_col = "pv_annuity") {
  annuity <- match.arg(annuity)
  timing <- match.arg(timing)
  interest_type <- match.arg(interest_type)

  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or tibble.", call. = FALSE)
  }

  .mc_assert_numeric_scalar(rate, "rate")
  .mc_assert_numeric_scalar(payment, "payment", min = 0)
  .mc_assert_positive_integer(payments_per_year, "payments_per_year")
  .mc_assert_numeric_scalar(m, "m", min = 0, strict_min = TRUE)
  .mc_assert_numeric_column(data, k_col, "k_col")
  .mc_assert_character_scalar(tx_col, "tx_col")
  .mc_assert_character_scalar(annuity_col, "annuity_col")

  if (annuity %in% c("temporary", "deferred_temporary", "certain")) {
    .mc_assert_numeric_scalar(term, "term", min = 0, strict_min = TRUE)
  }

  if (annuity %in% c("deferred", "deferred_temporary")) {
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

  needs_tx <- payments_per_year > 1 && annuity != "certain"

  if (needs_tx) {
    .mc_assert_numeric_column(data, tx_col, "tx_col")

    if (all(is.na(data[[tx_col]]))) {
      stop(
        "`tx_col` contains only missing values. Fractional life-contingent ",
        "payments require complete future lifetimes.",
        call. = FALSE
      )
    }
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

  Tx <- if (tx_col %in% names(data) && is.numeric(data[[tx_col]])) {
    data[[tx_col]]
  } else {
    rep(NA_real_, length(Kx))
  }

  if (needs_tx && any(Tx < 0, na.rm = TRUE)) {
    stop("`tx_col` must contain non-negative simulated lifetimes.", call. = FALSE)
  }

  payment_interval <- 1 / payments_per_year

  payment_times <- Map(
    function(k, tx) {
      if (payments_per_year == 1) {
        return(
          .mc_annuity_payment_times_annual(
            k = k,
            annuity = annuity,
            term = term,
            deferral_years = deferral_years,
            guarantee_years = guarantee_years,
            timing = timing
          )
        )
      }

      .mc_annuity_payment_times_fractional(
        tx = tx,
        annuity = annuity,
        term = term,
        deferral_years = deferral_years,
        guarantee_years = guarantee_years,
        timing = timing,
        payment_interval = payment_interval
      )
    },
    Kx,
    Tx
  )

  pv_annuity <- vapply(
    payment_times,
    function(times) {
      payment * sum(v^times)
    },
    numeric(1)
  )

  n_payments <- vapply(payment_times, length, integer(1))

  first_payment_time <- vapply(
    payment_times,
    function(times) {
      if (length(times) == 0) {
        return(NA_real_)
      }

      min(times)
    },
    numeric(1)
  )

  last_payment_time <- vapply(
    payment_times,
    function(times) {
      if (length(times) == 0) {
        return(NA_real_)
      }

      max(times)
    },
    numeric(1)
  )

  data |>
    dplyr::mutate(
      rate = rate,
      interest_type = interest_type,
      m = m,
      effective_rate = effective_rate,
      discount_factor = v,
      annuity = annuity,
      payment = payment,
      payments_per_year = payments_per_year,
      term = term,
      deferral_years = deferral_years,
      guarantee_years = guarantee_years,
      timing = timing,
      n_payments = n_payments,
      first_payment_time = first_payment_time,
      last_payment_time = last_payment_time,
      "{annuity_col}" := pv_annuity
    )
}


#' Internal helper: annual annuity payment times
#'
#' @param k Numeric scalar. Simulated curtate future lifetime.
#' @param annuity Character string. Annuity type.
#' @param term Numeric scalar or `NULL`.
#' @param deferral_years Numeric scalar.
#' @param guarantee_years Numeric scalar or `NULL`.
#' @param timing Character string. Either `"immediate"` or `"due"`.
#'
#' @return Numeric vector with annual payment times.
#'
#' @keywords internal
.mc_annuity_payment_times_annual <- function(k,
                                             annuity,
                                             term,
                                             deferral_years,
                                             guarantee_years,
                                             timing) {
  if (annuity == "whole_life") {
    if (timing == "immediate") {
      return(.mc_payment_times(1, k))
    }

    if (timing == "due") {
      return(.mc_payment_times(0, k))
    }
  }

  if (annuity == "temporary") {
    if (timing == "immediate") {
      return(.mc_payment_times(1, min(k, term)))
    }

    if (timing == "due") {
      return(.mc_payment_times(0, min(k, term - 1)))
    }
  }

  if (annuity == "deferred") {
    if (timing == "immediate") {
      return(.mc_payment_times(deferral_years + 1, k))
    }

    if (timing == "due") {
      return(.mc_payment_times(deferral_years, k))
    }
  }

  if (annuity == "deferred_temporary") {
    if (timing == "immediate") {
      return(
        .mc_payment_times(
          deferral_years + 1,
          min(k, deferral_years + term)
        )
      )
    }

    if (timing == "due") {
      return(
        .mc_payment_times(
          deferral_years,
          min(k, deferral_years + term - 1)
        )
      )
    }
  }

  if (annuity == "certain") {
    if (timing == "immediate") {
      return(.mc_payment_times(1, term))
    }

    if (timing == "due") {
      return(.mc_payment_times(0, term - 1))
    }
  }

  if (annuity == "guaranteed") {
    if (timing == "immediate") {
      return(.mc_payment_times(1, max(k, guarantee_years)))
    }

    if (timing == "due") {
      return(.mc_payment_times(0, max(k, guarantee_years - 1)))
    }
  }

  numeric(0)
}


#' Internal helper: fractional annuity payment times
#'
#' @param tx Numeric scalar. Simulated complete future lifetime.
#' @param annuity Character string. Annuity type.
#' @param term Numeric scalar or `NULL`.
#' @param deferral_years Numeric scalar.
#' @param guarantee_years Numeric scalar or `NULL`.
#' @param timing Character string. Either `"immediate"` or `"due"`.
#' @param payment_interval Numeric scalar. Time between payments.
#'
#' @return Numeric vector with fractional payment times.
#'
#' @keywords internal
.mc_annuity_payment_times_fractional <- function(tx,
                                                 annuity,
                                                 term,
                                                 deferral_years,
                                                 guarantee_years,
                                                 timing,
                                                 payment_interval) {
  eps <- sqrt(.Machine$double.eps)

  if (annuity == "certain") {
    if (timing == "immediate") {
      return(.mc_payment_times(payment_interval, term, by = payment_interval))
    }

    if (timing == "due") {
      return(
        .mc_payment_times(
          0,
          term - payment_interval,
          by = payment_interval
        )
      )
    }
  }

  if (is.na(tx) || !is.finite(tx) || tx < 0) {
    return(numeric(0))
  }

  last_alive_time <- if (timing == "immediate") {
    floor(tx / payment_interval + eps) * payment_interval
  } else {
    floor((tx - eps) / payment_interval) * payment_interval
  }

  if (last_alive_time < 0) {
    return(numeric(0))
  }

  if (annuity == "whole_life") {
    if (timing == "immediate") {
      return(
        .mc_payment_times(
          payment_interval,
          last_alive_time,
          by = payment_interval
        )
      )
    }

    if (timing == "due") {
      return(.mc_payment_times(0, last_alive_time, by = payment_interval))
    }
  }

  if (annuity == "temporary") {
    if (timing == "immediate") {
      last_time <- min(last_alive_time, term)

      return(
        .mc_payment_times(
          payment_interval,
          last_time,
          by = payment_interval
        )
      )
    }

    if (timing == "due") {
      last_time <- min(last_alive_time, term - payment_interval)

      return(.mc_payment_times(0, last_time, by = payment_interval))
    }
  }

  if (annuity == "deferred") {
    if (timing == "immediate") {
      first_time <- deferral_years + payment_interval

      return(.mc_payment_times(first_time, last_alive_time, by = payment_interval))
    }

    if (timing == "due") {
      first_time <- deferral_years

      return(.mc_payment_times(first_time, last_alive_time, by = payment_interval))
    }
  }

  if (annuity == "deferred_temporary") {
    if (timing == "immediate") {
      first_time <- deferral_years + payment_interval
      last_time <- min(last_alive_time, deferral_years + term)

      return(.mc_payment_times(first_time, last_time, by = payment_interval))
    }

    if (timing == "due") {
      first_time <- deferral_years
      last_time <- min(last_alive_time, deferral_years + term - payment_interval)

      return(.mc_payment_times(first_time, last_time, by = payment_interval))
    }
  }

  if (annuity == "guaranteed") {
    if (timing == "immediate") {
      guaranteed_last_time <- guarantee_years
      last_time <- max(last_alive_time, guaranteed_last_time)

      return(
        .mc_payment_times(
          payment_interval,
          last_time,
          by = payment_interval
        )
      )
    }

    if (timing == "due") {
      guaranteed_last_time <- guarantee_years - payment_interval
      last_time <- max(last_alive_time, guaranteed_last_time)

      return(.mc_payment_times(0, last_time, by = payment_interval))
    }
  }

  numeric(0)
}
