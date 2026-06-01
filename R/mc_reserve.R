#' Compute Monte Carlo prospective reserves for life contingencies
#'
#' Computes simulated prospective reserve losses at one or more policy
#' durations from simulated future lifetimes, using compact actuarial notation.
#'
#' This function recalculates future benefit and future premium present values
#' from each valuation duration. It is not a wrapper around
#' \code{\link{mc_loss}}, because reserves require valuing only the cash flows
#' that remain after the valuation time.
#'
#' For a policy in force at duration \eqn{t}, the simulated prospective loss is
#'
#' \deqn{
#'   L_t = Z_t - P Y_t,
#' }
#'
#' where \eqn{Z_t} is the present value at duration \eqn{t} of future benefits,
#' \eqn{Y_t} is the present value at duration \eqn{t} of future premium
#' payments, and \eqn{P} is the premium per payment.
#'
#' @param .data A data frame or tibble containing simulated future lifetimes,
#'   typically returned by \code{\link{simulate_lifetime}} or by
#'   \code{\link{mc_multilife_status}}.
#' @param t Numeric vector. Policy duration or durations at which the reserve
#'   is evaluated. Default is \code{0}.
#' @param i Numeric scalar. Interest-rate input used for discounting.
#' @param P Optional numeric scalar. Premium used in the reserve loss. If
#'   \code{NULL}, the function first tries to use \code{col_P} from
#'   \code{.data}. If \code{col_P} is not found, a Monte Carlo net premium at
#'   issue is estimated internally as \eqn{\bar{Z}_0 / \bar{Y}_0}.
#' @param col_P Character string. Name of the premium column in \code{.data}.
#'   Default is \code{"P"}. If this column is not found and \code{col_P = "P"},
#'   the legacy column \code{"premium"} is used when available.
#' @param benefit Numeric scalar. Benefit amount payable under the insurance.
#'   Default is \code{1}.
#' @param payment Numeric scalar. Amount of each premium annuity payment.
#'   Default is \code{1}.
#' @param k Positive integer. Number of premium payments per year. Default is
#'   \code{1}, corresponding to annual premiums.
#' @param type Character string specifying the insurance type. Canonical
#'   options are \code{"whole"}, \code{"term"}, \code{"deferred"},
#'   \code{"deferred_term"}, \code{"pure_endowment"}, and
#'   \code{"endowment"}. Transitional aliases \code{"whole_life"} and
#'   \code{"deferred_temporary"} are also accepted.
#' @param annuity_type Character string specifying the premium annuity type.
#'   Canonical options are \code{"whole"}, \code{"temporary"},
#'   \code{"deferred"}, \code{"deferred_temporary"}, \code{"certain"}, and
#'   \code{"guaranteed"}. The transitional alias \code{"whole_life"} is also
#'   accepted.
#' @param n Numeric scalar. Contract term in years. Required for insurance
#'   types \code{"term"}, \code{"deferred_term"}, \code{"pure_endowment"},
#'   and \code{"endowment"}, and for annuity types \code{"temporary"},
#'   \code{"deferred_temporary"}, and \code{"certain"}.
#' @param h Numeric scalar. Deferral period in years. Default is \code{0}.
#' @param n_guar Numeric scalar. Guaranteed payment period in years. Required
#'   when \code{annuity_type = "guaranteed"}.
#' @param timing Character string specifying when death benefits are paid.
#'   Available options are \code{"end_of_year"} and
#'   \code{"moment_of_death"}. Default is \code{"end_of_year"}.
#' @param premium_timing Character string specifying the premium payment timing.
#'   Available options are \code{"due"} and \code{"immediate"}. Default is
#'   \code{"due"}.
#' @param reserve_timing Character string specifying whether payments due
#'   exactly at the valuation duration are included. Available options are
#'   \code{"before_payment"} and \code{"after_payment"}.
#' @param i_type Character string specifying the interest-rate convention.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}. The transitional value
#'   \code{"nominal"} is accepted and treated as \code{"nominal_interest"}.
#' @param m Positive integer. Number of interest conversion periods per year
#'   for nominal annual rates. Default is \code{1}. This argument controls the
#'   interest-rate conversion frequency only. It does not represent premium
#'   payment frequency.
#' @param col_K Character string. Name of the column containing simulated
#'   curtate future lifetimes. Default is \code{"Kx"}.
#' @param col_T Character string. Name of the column containing simulated
#'   complete future lifetimes. Default is \code{"Tx"}.
#' @param in_force_basis Character string specifying how the in-force indicator
#'   is evaluated. Available options are \code{"auto"}, \code{"complete"},
#'   and \code{"curtate"}.
#' @param not_in_force Character string specifying what to return for scenarios
#'   that are not in force at the valuation duration. Available options are
#'   \code{"na"} and \code{"zero"}.
#' @param col_L Character string. Name of the output reserve loss column.
#'   Default is \code{"L_t"}.
#' @param ... Transitional compatibility for older calls using \code{data},
#'   \code{duration}, \code{rate}, \code{premium},
#'   \code{premium_col}, \code{payments_per_year}, \code{insurance},
#'   \code{annuity}, \code{term}, \code{deferral_years},
#'   \code{guarantee_years}, \code{payment_timing},
#'   \code{interest_type}, \code{k_col}, \code{tx_col}, and
#'   \code{reserve_col}.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{t} is the valuation duration, \code{i} is the
#' interest-rate input, \code{i_type} is the interest-rate convention, \code{m}
#' is the nominal conversion frequency, \code{k} is the premium payment
#' frequency, \code{n} is the contract term, \code{h} is the deferral period,
#' \code{P} is the premium per payment, and \code{L_t} is the simulated
#' prospective loss at duration \code{t}.
#'
#' The arguments \code{m} and \code{k} have deliberately different meanings:
#' \itemize{
#'   \item \code{m} is used only for nominal interest-rate conversion.
#'   \item \code{k} controls how frequently future premium payments are made.
#' }
#'
#' The argument \code{payment} represents the amount of each premium annuity
#' payment used to construct \eqn{Y_t}. Thus, for monthly premiums with total
#' annual premium equal to 1, use \code{payment = 1 / 12} and \code{k = 12}.
#'
#' If \code{k = 1}, future premiums are annual. If \code{k > 1}, future
#' premiums are made at fractional times, and a valid complete future lifetime
#' column supplied through \code{col_T} is required for life-contingent premium
#' annuities.
#'
#' Durations may be integer or fractional. Fractional reserve durations require
#' complete future lifetimes. For example, monthly reserve calculations may use
#' \code{t = seq(0, 20, by = 1 / 12)}.
#'
#' If \code{reserve_timing = "before_payment"}, cash flows occurring exactly at
#' the valuation duration are included. If
#' \code{reserve_timing = "after_payment"}, cash flows occurring exactly at the
#' valuation duration are excluded.
#'
#' If \code{not_in_force = "na"}, scenarios that are not in force at a given
#' duration receive \code{NA} values for future present values and reserve
#' losses. This is useful for estimating reserves conditional on the policy
#' still being in force. If \code{not_in_force = "zero"}, those scenarios
#' receive zero values, which may be useful for portfolio run-off summaries.
#'
#' This function computes prospective reserves under the simulated model. It
#' does not include expenses, surrender values, taxes, profit loadings, or
#' statutory reserving adjustments.
#'
#' @return A tibble with one row per original simulation and per requested
#' duration. It contains the original simulation columns and additional columns
#' including \code{t}, \code{in_force}, \code{i}, \code{i_type}, \code{m},
#' \code{i_effective}, \code{v}, \code{type}, \code{annuity_type},
#' \code{benefit}, \code{payment}, \code{k}, \code{Z_t}, \code{Y_t},
#' \code{P}, and \code{L_t} or another name supplied through \code{col_L}.
#'
#' For transition, the output also includes legacy columns such as
#' \code{duration}, \code{rate}, \code{interest_type}, \code{effective_rate},
#' \code{discount_factor}, \code{insurance}, \code{annuity},
#' \code{payments_per_year}, \code{term}, \code{deferral_years},
#' \code{guarantee_years}, \code{future_pv_benefit},
#' \code{future_pv_premiums}, \code{premium}, and \code{reserve_loss}.
#'
#' @seealso
#' \code{\link{simulate_lifetime}}, \code{\link{simulate_lifetimes}},
#' \code{\link{mc_multilife_status}}, \code{\link{mc_insurance}},
#' \code{\link{mc_annuity}}, \code{\link{mc_premium}}, \code{\link{mc_loss}},
#' \code{\link{summary_mc}}
#'
#' @references
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A.,
#' and Nesbitt, C. J. (1997). \emph{Actuarial Mathematics}. Second Edition.
#' Society of Actuaries.
#'
#' @family monte-carlo
#'
#' @examples
#' lt <- tibble::tibble(
#'   x = 40:100,
#'   qx = seq(0.002, 1, length.out = 61)
#' )
#'
#' # Annual prospective reserves for whole-life insurance
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     seed = 123
#'   ) |>
#'   mc_reserve(
#'     t = c(0, 5, 10),
#'     i = 0.05,
#'     type = "whole",
#'     annuity_type = "whole",
#'     benefit = 1,
#'     payment = 1,
#'     k = 1,
#'     premium_timing = "due"
#'   )
#'
#' # Monthly premium reserve curve at annual valuation durations
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     frac = "udd",
#'     seed = 123
#'   ) |>
#'   mc_reserve(
#'     t = c(0, 1, 2, 3),
#'     i = 0.05,
#'     type = "whole",
#'     annuity_type = "whole",
#'     benefit = 1,
#'     payment = 1 / 12,
#'     k = 12,
#'     premium_timing = "due"
#'   )
#'
#' # Fractional reserve durations
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     frac = "udd",
#'     seed = 123
#'   ) |>
#'   mc_reserve(
#'     t = seq(0, 1, by = 1 / 4),
#'     i = 0.05,
#'     type = "whole",
#'     annuity_type = "whole",
#'     benefit = 1,
#'     payment = 1 / 12,
#'     k = 12,
#'     premium_timing = "due"
#'   ) |>
#'   summary_mc(value_col = "L_t", by = "t")
#'
#' # Joint-life reserve using a multiple-life status
#' lt |>
#'   simulate_lifetimes(
#'     x = c(60, 58),
#'     n_sim = 25,
#'     seed = 123
#'   ) |>
#'   mc_multilife_status(status = "joint") |>
#'   mc_reserve(
#'     t = c(0, 5),
#'     i = 0.04,
#'     type = "whole",
#'     annuity_type = "whole",
#'     benefit = 1,
#'     payment = 1,
#'     k = 1,
#'     col_K = "K_status",
#'     col_T = "T_status"
#'   )
#'
#' @export
mc_reserve <- function(
    .data = NULL,
    t = 0,
    i = NULL,
    P = NULL,
    col_P = "P",
    benefit = 1,
    payment = 1,
    k = 1L,
    type = c(
      "whole",
      "term",
      "deferred",
      "deferred_term",
      "pure_endowment",
      "endowment",
      "whole_life",
      "deferred_temporary"
    ),
    annuity_type = c(
      "whole",
      "temporary",
      "deferred",
      "deferred_temporary",
      "certain",
      "guaranteed",
      "whole_life"
    ),
    n = NULL,
    h = 0,
    n_guar = NULL,
    timing = c("end_of_year", "moment_of_death"),
    premium_timing = c("due", "immediate"),
    reserve_timing = c("before_payment", "after_payment"),
    i_type = c(
      "effective",
      "nominal_interest",
      "nominal_discount",
      "force",
      "nominal"
    ),
    m = 1,
    col_K = "Kx",
    col_T = "Tx",
    in_force_basis = c("auto", "complete", "curtate"),
    not_in_force = c("na", "zero"),
    col_L = "L_t",
    ...
) {
  dots <- list(...)
  t_missing <- missing(t)
  type_missing <- missing(type)
  annuity_type_missing <- missing(annuity_type)
  i_type_missing <- missing(i_type)

  # -------------------------------------------------------------------------
  # Transitional compatibility with the previous public API
  # -------------------------------------------------------------------------

  allowed_old <- c(
    "data",
    "duration",
    "rate",
    "premium",
    "premium_col",
    "payments_per_year",
    "insurance",
    "annuity",
    "term",
    "deferral_years",
    "guarantee_years",
    "payment_timing",
    "interest_type",
    "k_col",
    "tx_col",
    "reserve_col"
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

  if (!is.null(dots$data)) {
    if (!is.null(.data)) {
      stop("Provide only one of `.data` or deprecated `data`.", call. = FALSE)
    }

    .data <- dots$data
  }

  if (!is.null(dots$duration)) {
    if (!t_missing) {
      stop("Provide only one of `t` or deprecated `duration`.", call. = FALSE)
    }

    t <- dots$duration
  }

  if (!is.null(dots$rate)) {
    if (!is.null(i)) {
      stop("Provide only one of `i` or deprecated `rate`.", call. = FALSE)
    }

    i <- dots$rate
  }

  if (!is.null(dots$premium)) {
    if (!is.null(P)) {
      stop("Provide only one of `P` or deprecated `premium`.", call. = FALSE)
    }

    P <- dots$premium
  }

  if (!is.null(dots$premium_col)) {
    if (!identical(col_P, "P")) {
      stop("Provide only one of `col_P` or deprecated `premium_col`.",
           call. = FALSE)
    }

    col_P <- dots$premium_col
  }

  if (!is.null(dots$payments_per_year)) {
    if (!identical(k, 1L) && !identical(k, 1)) {
      stop("Provide only one of `k` or deprecated `payments_per_year`.",
           call. = FALSE)
    }

    k <- dots$payments_per_year
  }

  if (!is.null(dots$insurance)) {
    if (!type_missing) {
      stop("Provide only one of `type` or deprecated `insurance`.", call. = FALSE)
    }

    type <- dots$insurance
  }

  if (!is.null(dots$annuity)) {
    if (!annuity_type_missing) {
      stop("Provide only one of `annuity_type` or deprecated `annuity`.",
           call. = FALSE)
    }

    annuity_type <- dots$annuity
  }

  if (!is.null(dots$term)) {
    if (!is.null(n)) {
      stop("Provide only one of `n` or deprecated `term`.", call. = FALSE)
    }

    n <- dots$term
  }

  if (!is.null(dots$deferral_years)) {
    if (!identical(h, 0) && !identical(h, 0L)) {
      stop("Provide only one of `h` or deprecated `deferral_years`.", call. = FALSE)
    }

    h <- dots$deferral_years
  }

  if (!is.null(dots$guarantee_years)) {
    if (!is.null(n_guar)) {
      stop("Provide only one of `n_guar` or deprecated `guarantee_years`.",
           call. = FALSE)
    }

    n_guar <- dots$guarantee_years
  }

  if (!is.null(dots$payment_timing)) {
    if (!missing(timing)) {
      stop("Provide only one of `timing` or deprecated `payment_timing`.",
           call. = FALSE)
    }

    timing <- dots$payment_timing
  }

  if (!is.null(dots$interest_type)) {
    if (!i_type_missing) {
      stop("Provide only one of `i_type` or deprecated `interest_type`.",
           call. = FALSE)
    }

    i_type <- dots$interest_type
  }

  if (!is.null(dots$k_col)) {
    if (!identical(col_K, "Kx")) {
      stop("Provide only one of `col_K` or deprecated `k_col`.", call. = FALSE)
    }

    col_K <- dots$k_col
  }

  if (!is.null(dots$tx_col)) {
    if (!identical(col_T, "Tx")) {
      stop("Provide only one of `col_T` or deprecated `tx_col`.", call. = FALSE)
    }

    col_T <- dots$tx_col
  }

  if (!is.null(dots$reserve_col)) {
    if (!identical(col_L, "L_t")) {
      stop("Provide only one of `col_L` or deprecated `reserve_col`.",
           call. = FALSE)
    }

    col_L <- dots$reserve_col
  }

  type <- match.arg(
    type,
    choices = c(
      "whole",
      "term",
      "deferred",
      "deferred_term",
      "pure_endowment",
      "endowment",
      "whole_life",
      "deferred_temporary"
    )
  )

  if (identical(type, "whole_life")) {
    type <- "whole"
  }

  if (identical(type, "deferred_temporary")) {
    type <- "deferred_term"
  }

  annuity_type <- match.arg(
    annuity_type,
    choices = c(
      "whole",
      "temporary",
      "deferred",
      "deferred_temporary",
      "certain",
      "guaranteed",
      "whole_life"
    )
  )

  if (identical(annuity_type, "whole_life")) {
    annuity_type <- "whole"
  }

  timing <- match.arg(timing)
  premium_timing <- match.arg(premium_timing)
  reserve_timing <- match.arg(reserve_timing)
  i_type <- match.arg(i_type)
  in_force_basis <- match.arg(in_force_basis)
  not_in_force <- match.arg(not_in_force)

  if (identical(i_type, "nominal")) {
    i_type <- "nominal_interest"
  }

  # -------------------------------------------------------------------------
  # Validation
  # -------------------------------------------------------------------------

  if (!is.data.frame(.data)) {
    stop("`.data` must be a data frame or tibble.", call. = FALSE)
  }

  if (!is.numeric(t) ||
      length(t) == 0L ||
      anyNA(t) ||
      any(!is.finite(t)) ||
      any(t < 0)) {
    stop("`t` must be a non-negative numeric vector.", call. = FALSE)
  }

  if (is.null(i)) {
    stop("`i` must be provided.", call. = FALSE)
  }

  .mc_assert_numeric_scalar(i, "i")
  .mc_assert_numeric_scalar(benefit, "benefit", min = 0)
  .mc_assert_numeric_scalar(payment, "payment", min = 0)
  .mc_assert_positive_integer(k, "k")
  .mc_assert_numeric_scalar(m, "m", min = 0, strict_min = TRUE)
  .mc_assert_numeric_column(.data, col_K, "col_K")
  .mc_assert_character_scalar(col_P, "col_P")
  .mc_assert_character_scalar(col_L, "col_L")
  .mc_assert_character_scalar(col_T, "col_T")

  if (type %in% c("term", "deferred_term", "pure_endowment", "endowment") ||
      annuity_type %in% c("temporary", "deferred_temporary", "certain")) {
    .mc_assert_numeric_scalar(n, "n", min = 0, strict_min = TRUE)
  }

  if (type %in% c("deferred", "deferred_term") ||
      annuity_type %in% c("deferred", "deferred_temporary")) {
    .mc_assert_numeric_scalar(h, "h", min = 0, strict_min = TRUE)
  } else {
    .mc_assert_numeric_scalar(h, "h", min = 0)
  }

  if (annuity_type == "guaranteed") {
    .mc_assert_numeric_scalar(
      n_guar,
      "n_guar",
      min = 0,
      strict_min = TRUE
    )
  } else {
    .mc_assert_numeric_scalar(
      n_guar,
      "n_guar",
      min = 0,
      allow_null = TRUE
    )
  }

  if (!is.null(P)) {
    .mc_assert_numeric_scalar(P, "P")
  }

  k <- as.integer(round(k))
  m <- as.integer(round(m))

  Kx <- .data[[col_K]]

  if (any(Kx < 0, na.rm = TRUE)) {
    stop("`col_K` must contain non-negative simulated lifetimes.", call. = FALSE)
  }

  has_T <- col_T %in% names(.data) &&
    is.numeric(.data[[col_T]]) &&
    !all(is.na(.data[[col_T]]))

  if (timing == "moment_of_death" && !has_T) {
    stop(
      "Benefits payable at the moment of death require a valid complete ",
      "future lifetime column identified by `col_T`.",
      call. = FALSE
    )
  }

  if (k > 1L && annuity_type != "certain" && !has_T) {
    stop(
      "Fractional life-contingent premium payments require a valid complete ",
      "future lifetime column identified by `col_T`.",
      call. = FALSE
    )
  }

  has_fractional_t <- any(
    abs(t - round(t)) > sqrt(.Machine$double.eps)
  )

  if (has_fractional_t && !has_T) {
    stop(
      "Fractional reserve durations require a valid complete future lifetime ",
      "column identified by `col_T`.",
      call. = FALSE
    )
  }

  if (in_force_basis == "auto") {
    in_force_basis <- if (has_T) "complete" else "curtate"
  }

  if (in_force_basis == "complete") {
    .mc_assert_numeric_column(.data, col_T, "col_T")
  }

  if (in_force_basis == "curtate" && has_fractional_t) {
    stop(
      "Fractional reserve durations require `in_force_basis = 'complete' ",
      "and a valid `col_T`.",
      call. = FALSE
    )
  }

  i_effective <- .mc_effective_rate(
    i = i,
    i_type = i_type,
    m = m
  )

  if (!is.finite(i_effective) || i_effective <= -1) {
    stop(
      "The annual effective interest rate implied by `i`, `i_type`, and `m` ",
      "must be greater than -1.",
      call. = FALSE
    )
  }

  v <- 1 / (1 + i_effective)

  Tx <- if (has_T) {
    .data[[col_T]]
  } else {
    rep(NA_real_, length(Kx))
  }

  death_time <- if (timing == "end_of_year") {
    Kx + 1
  } else {
    Tx
  }

  benefit_info <- .mc_reserve_benefit_info(
    Kx = Kx,
    Tx = Tx,
    death_time = death_time,
    type = type,
    n = n,
    h = h,
    timing = timing
  )

  premium_payment_times <- Map(
    function(K, T) {
      if (k == 1L) {
        return(
          .mc_annuity_payment_times_annual(
            K = K,
            type = annuity_type,
            n = n,
            h = h,
            n_guar = n_guar,
            timing = premium_timing
          )
        )
      }

      .mc_annuity_payment_times_fractional(
        T = T,
        type = annuity_type,
        n = n,
        h = h,
        n_guar = n_guar,
        timing = premium_timing,
        payment_interval = 1 / k
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

  if (!is.null(P)) {
    P_used <- rep(P, length(Kx))
  } else {
    if (!col_P %in% names(.data) &&
        identical(col_P, "P") &&
        "premium" %in% names(.data)) {
      col_P <- "premium"
    }

    if (col_P %in% names(.data)) {
      if (!is.numeric(.data[[col_P]])) {
        stop("`col_P` must identify a numeric column.", call. = FALSE)
      }

      P_used <- .data[[col_P]]
    } else {
      denominator <- mean(issue_pv_premiums)

      if (is.na(denominator) || denominator <= 0) {
        stop(
          "The mean simulated premium annuity present value must be positive ",
          "when estimating the net premium internally.",
          call. = FALSE
        )
      }

      P_hat <- mean(issue_pv_benefit) / denominator
      P_used <- rep(P_hat, length(Kx))
    }
  }

  row_id <- seq_len(nrow(.data))
  n_t <- length(t)

  expanded <- tibble::as_tibble(.data)[
    rep(row_id, each = n_t),
    ,
    drop = FALSE
  ]

  row_rep <- rep(row_id, each = n_t)
  t_rep <- rep(t, times = length(row_id))

  Kx_rep <- Kx[row_rep]
  Tx_rep <- Tx[row_rep]
  P_rep <- P_used[row_rep]

  in_force <- if (in_force_basis == "complete") {
    Tx_rep > t_rep
  } else {
    Kx_rep >= t_rep
  }

  benefit_time_rep <- benefit_info$benefit_time[row_rep]
  benefit_indicator_rep <- benefit_info$benefit_indicator[row_rep]

  include_benefit <- if (reserve_timing == "before_payment") {
    benefit_time_rep >= t_rep
  } else {
    benefit_time_rep > t_rep
  }

  Z_t <- ifelse(
    benefit_indicator_rep == 1 & include_benefit,
    benefit * v^(benefit_time_rep - t_rep),
    0
  )

  Y_t <- vapply(
    seq_along(row_rep),
    function(idx) {
      times <- premium_payment_times[[row_rep[idx]]]

      future_times <- if (reserve_timing == "before_payment") {
        times[times >= t_rep[idx]]
      } else {
        times[times > t_rep[idx]]
      }

      payment * sum(v^(future_times - t_rep[idx]))
    },
    numeric(1)
  )

  if (not_in_force == "na") {
    Z_t[!in_force] <- NA_real_
    Y_t[!in_force] <- NA_real_
  } else {
    Z_t[!in_force] <- 0
    Y_t[!in_force] <- 0
  }

  L_t <- Z_t - P_rep * Y_t

  out <- expanded |>
    dplyr::mutate(
      t = t_rep,
      duration = t_rep,
      in_force = in_force,
      reserve_timing = reserve_timing,
      not_in_force = not_in_force,
      i = i,
      rate = i,
      i_type = i_type,
      interest_type = i_type,
      m = m,
      i_effective = i_effective,
      effective_rate = i_effective,
      v = v,
      discount_factor = v,
      type = type,
      insurance = ifelse(type == "whole", "whole_life", type),
      annuity_type = annuity_type,
      annuity = ifelse(annuity_type == "whole", "whole_life", annuity_type),
      benefit = benefit,
      payment = payment,
      k = k,
      payments_per_year = k,
      n = if (is.null(n)) NA_real_ else n,
      term = if (is.null(n)) NA_real_ else n,
      h = h,
      deferral_years = h,
      n_guar = if (is.null(n_guar)) NA_real_ else n_guar,
      guarantee_years = if (is.null(n_guar)) NA_real_ else n_guar,
      timing = timing,
      payment_timing = timing,
      premium_timing = premium_timing,
      Z_t = Z_t,
      Y_t = Y_t,
      future_pv_benefit = Z_t,
      future_pv_premiums = Y_t,
      P = P_rep,
      premium = P_rep,
      "{col_L}" := L_t
    )

  if (!identical(col_L, "reserve_loss") && !"reserve_loss" %in% names(out)) {
    out <- out |>
      dplyr::mutate(
        reserve_loss = L_t
      )
  }

  out
}


#' Internal helper: benefit payment information for reserve calculations
#'
#' Determines whether a benefit is paid and at what time for each simulated
#' scenario.
#'
#' @param Kx Numeric vector of simulated curtate future lifetimes.
#' @param Tx Numeric vector of simulated complete future lifetimes.
#' @param death_time Numeric vector of death benefit payment times.
#' @param type Character string specifying the insurance type.
#' @param n Numeric scalar or `NULL`.
#' @param h Numeric scalar.
#' @param timing Character string. Either `"end_of_year"` or
#'   `"moment_of_death"`.
#'
#' @return A list with `benefit_indicator` and `benefit_time`.
#'
#' @keywords internal
.mc_reserve_benefit_info <- function(
    Kx,
    Tx,
    death_time,
    type,
    n,
    h,
    timing
) {
  benefit_indicator <- rep(0, length(Kx))
  benefit_time <- rep(NA_real_, length(Kx))

  if (type == "whole") {
    benefit_indicator <- rep(1, length(Kx))
    benefit_time <- death_time
  }

  if (type == "term") {
    if (timing == "end_of_year") {
      benefit_indicator <- as.numeric(Kx < n)
    } else {
      benefit_indicator <- as.numeric(Tx <= n)
    }

    benefit_time <- ifelse(benefit_indicator == 1, death_time, NA_real_)
  }

  if (type == "deferred") {
    if (timing == "end_of_year") {
      benefit_indicator <- as.numeric(Kx >= h)
    } else {
      benefit_indicator <- as.numeric(Tx > h)
    }

    benefit_time <- ifelse(benefit_indicator == 1, death_time, NA_real_)
  }

  if (type == "deferred_term") {
    if (timing == "end_of_year") {
      benefit_indicator <- as.numeric(
        Kx >= h & Kx < h + n
      )
    } else {
      benefit_indicator <- as.numeric(
        Tx > h & Tx <= h + n
      )
    }

    benefit_time <- ifelse(benefit_indicator == 1, death_time, NA_real_)
  }

  if (type == "pure_endowment") {
    if (timing == "end_of_year") {
      benefit_indicator <- as.numeric(Kx >= n)
    } else {
      benefit_indicator <- as.numeric(Tx >= n)
    }

    benefit_time <- ifelse(benefit_indicator == 1, n, NA_real_)
  }

  if (type == "endowment") {
    if (timing == "end_of_year") {
      death_before_term <- Kx < n
    } else {
      death_before_term <- Tx <= n
    }

    benefit_indicator <- rep(1, length(Kx))
    benefit_time <- ifelse(death_before_term, death_time, n)
  }

  list(
    benefit_indicator = benefit_indicator,
    benefit_time = benefit_time
  )
}
