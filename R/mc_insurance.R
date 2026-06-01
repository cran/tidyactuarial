#' Compute simulated present values for life insurance benefits
#'
#' Computes Monte Carlo simulated present values of life insurance benefits from
#' simulated future lifetimes, using compact actuarial notation.
#'
#' This function is designed to be used after \code{\link{simulate_lifetime}},
#' \code{\link{simulate_lifetimes}}, or \code{\link{mc_multilife_status}}. It
#' takes simulated values of the curtate future lifetime \eqn{K_x}, and when
#' needed the complete future lifetime \eqn{T_x}, and evaluates the present
#' value random variable associated with classical life insurance benefits.
#'
#' @param .data A data frame or tibble containing simulated future lifetimes,
#'   typically returned by \code{\link{simulate_lifetime}} or by
#'   \code{\link{mc_multilife_status}} when working with multiple-life statuses.
#' @param i Numeric scalar. Interest-rate input used for discounting.
#' @param benefit Numeric scalar. Benefit amount payable under the insurance.
#'   Default is \code{1}.
#' @param type Character string specifying the insurance type. Canonical
#'   options are \code{"whole"}, \code{"term"}, \code{"deferred"},
#'   \code{"deferred_term"}, \code{"pure_endowment"}, and \code{"endowment"}.
#'   Transitional aliases \code{"whole_life"} and \code{"deferred_temporary"}
#'   are also accepted.
#' @param n Numeric scalar. Insurance term in years. Required for
#'   \code{"term"}, \code{"deferred_term"}, \code{"pure_endowment"}, and
#'   \code{"endowment"}.
#' @param h Numeric scalar. Deferral period in years. Default is \code{0}.
#'   Required to be positive for \code{"deferred"} and
#'   \code{"deferred_term"} insurance.
#' @param timing Character string specifying when death benefits are paid.
#'   Available options are \code{"end_of_year"} and
#'   \code{"moment_of_death"}. Default is \code{"end_of_year"}.
#' @param i_type Character string specifying the interest-rate convention.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}. The transitional value
#'   \code{"nominal"} is accepted and treated as
#'   \code{"nominal_interest"}.
#' @param m Positive integer. Number of interest conversion periods per year
#'   for nominal annual rates. Default is \code{1}. This argument controls the
#'   interest-rate conversion frequency only. It does not represent benefit
#'   frequency or premium payment frequency.
#' @param col_K Character string. Name of the column containing simulated
#'   curtate future lifetimes. Default is \code{"Kx"}.
#' @param col_T Character string. Name of the column containing simulated
#'   complete future lifetimes. Required when
#'   \code{timing = "moment_of_death"}. Default is \code{"Tx"}.
#' @param col_pv Character string. Name of the output column containing
#'   simulated present values of benefits. Default is \code{"pv_benefit"}.
#' @param ... Transitional compatibility for older calls using \code{data},
#'   \code{rate}, \code{insurance}, \code{term}, \code{deferral_years},
#'   \code{payment_timing}, \code{interest_type}, \code{k_col},
#'   \code{tx_col}, and \code{benefit_col}.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{i} is the interest-rate input, \code{i_type}
#' is the interest-rate convention, \code{m} is the nominal conversion
#' frequency, \code{n} is the insurance term, and \code{h} is the deferral
#' period.
#'
#' The arguments \code{m} and \code{col_K} have deliberately different roles:
#' \itemize{
#'   \item \code{m} is used only for nominal interest-rate conversion.
#'   \item \code{col_K} identifies the simulated curtate future lifetime column.
#' }
#'
#' If \code{timing = "end_of_year"}, death benefits are discounted using
#' \eqn{K_x + 1}. If \code{timing = "moment_of_death"}, death benefits are
#' discounted using \eqn{T_x}.
#'
#' The following insurance types are supported:
#' \itemize{
#'   \item \code{"whole"}: benefit is paid whenever death occurs.
#'   \item \code{"term"}: benefit is paid if death occurs within \code{n}
#'   years.
#'   \item \code{"deferred"}: benefit is paid if death occurs after the
#'   deferral period \code{h}.
#'   \item \code{"deferred_term"}: benefit is paid if death occurs after
#'   \code{h} and within the following \code{n} years.
#'   \item \code{"pure_endowment"}: benefit is paid at time \code{n} if the
#'   life survives to that time.
#'   \item \code{"endowment"}: death benefit is paid if death occurs within
#'   \code{n} years; otherwise, a survival benefit is paid at time \code{n}.
#' }
#'
#' The function returns simulated present values, not only their expected value.
#' Therefore the resulting column can be summarized with
#' \code{\link{summary_mc}}, plotted with \code{ggplot2}, or used to construct
#' premiums, losses, and reserves.
#'
#' @return A tibble with the original simulation columns and additional columns:
#' \describe{
#'   \item{i}{Original interest-rate input.}
#'   \item{i_type}{Interest-rate convention.}
#'   \item{m}{Interest conversion frequency.}
#'   \item{i_effective}{Equivalent annual effective interest rate.}
#'   \item{v}{Annual discount factor.}
#'   \item{type}{Canonical insurance type.}
#'   \item{benefit}{Benefit amount.}
#'   \item{n}{Insurance term, if applicable.}
#'   \item{h}{Deferral period.}
#'   \item{timing}{Timing used for death benefits.}
#'   \item{benefit_time}{Simulated payment time of the benefit.}
#'   \item{benefit_indicator}{Indicator that the benefit is paid.}
#'   \item{pv_benefit}{Simulated present value of the benefit, or another name
#'   supplied through \code{col_pv}.}
#' }
#'
#' For transition, the output also includes legacy columns such as
#' \code{rate}, \code{interest_type}, \code{effective_rate},
#' \code{discount_factor}, \code{insurance}, \code{term},
#' \code{deferral_years}, and \code{payment_timing}.
#'
#' @seealso
#' \code{\link{simulate_lifetime}}, \code{\link{simulate_lifetimes}},
#' \code{\link{mc_multilife_status}}, \code{\link{mc_annuity}},
#' \code{\link{mc_premium}}, \code{\link{mc_loss}}, \code{\link{mc_reserve}},
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
#' # Whole-life insurance payable at the end of the year of death
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
#'   )
#'
#' # 20-year term insurance
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     seed = 123
#'   ) |>
#'   mc_insurance(
#'     i = 0.05,
#'     type = "term",
#'     n = 20,
#'     benefit = 100000
#'   )
#'
#' # 10-year deferred whole-life insurance
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     seed = 123
#'   ) |>
#'   mc_insurance(
#'     i = 0.05,
#'     type = "deferred",
#'     h = 10,
#'     benefit = 1
#'   )
#'
#' # Endowment insurance payable at the moment of death if death occurs
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     frac = "udd",
#'     seed = 123
#'   ) |>
#'   mc_insurance(
#'     i = 0.05,
#'     type = "endowment",
#'     n = 20,
#'     timing = "moment_of_death",
#'     benefit = 1
#'   )
#'
#' # Nominal interest rate convertible monthly
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     seed = 123
#'   ) |>
#'   mc_insurance(
#'     i = 0.06,
#'     i_type = "nominal_interest",
#'     m = 12,
#'     type = "whole",
#'     benefit = 1
#'   )
#'
#' # First-death insurance using a multiple-life status
#' lt |>
#'   simulate_lifetimes(
#'     x = c(60, 58),
#'     n_sim = 25,
#'     frac = "udd",
#'     seed = 123
#'   ) |>
#'   mc_multilife_status(status = "joint") |>
#'   mc_insurance(
#'     i = 0.04,
#'     type = "whole",
#'     benefit = 100000,
#'     col_K = "K_status",
#'     col_T = "T_status"
#'   )
#'
#' @export
mc_insurance <- function(
    .data = NULL,
    i = NULL,
    benefit = 1,
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
    n = NULL,
    h = 0,
    timing = c("end_of_year", "moment_of_death"),
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
    col_pv = "pv_benefit",
    ...
) {
  dots <- list(...)
  type_missing <- missing(type)
  i_type_missing <- missing(i_type)

  # -------------------------------------------------------------------------
  # Transitional compatibility with the previous public API
  # -------------------------------------------------------------------------

  allowed_old <- c(
    "data",
    "rate",
    "insurance",
    "term",
    "deferral_years",
    "payment_timing",
    "interest_type",
    "k_col",
    "tx_col",
    "benefit_col"
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

  if (!is.null(dots$rate)) {
    if (!is.null(i)) {
      stop("Provide only one of `i` or deprecated `rate`.", call. = FALSE)
    }

    i <- dots$rate
  }

  if (!is.null(dots$insurance)) {
    if (!type_missing) {
      stop("Provide only one of `type` or deprecated `insurance`.", call. = FALSE)
    }

    type <- dots$insurance
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

  if (!is.null(dots$benefit_col)) {
    if (!identical(col_pv, "pv_benefit")) {
      stop("Provide only one of `col_pv` or deprecated `benefit_col`.",
           call. = FALSE)
    }

    col_pv <- dots$benefit_col
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

  timing <- match.arg(timing)
  i_type <- match.arg(i_type)

  if (identical(i_type, "nominal")) {
    i_type <- "nominal_interest"
  }

  # -------------------------------------------------------------------------
  # Validation
  # -------------------------------------------------------------------------

  if (!is.data.frame(.data)) {
    stop("`.data` must be a data frame or tibble.", call. = FALSE)
  }

  if (is.null(i)) {
    stop("`i` must be provided.", call. = FALSE)
  }

  .mc_assert_numeric_scalar(i, "i")
  .mc_assert_numeric_scalar(benefit, "benefit", min = 0)
  .mc_assert_numeric_scalar(m, "m", min = 0, strict_min = TRUE)
  .mc_assert_numeric_column(.data, col_K, "col_K")
  .mc_assert_character_scalar(col_T, "col_T")
  .mc_assert_character_scalar(col_pv, "col_pv")

  if (timing == "moment_of_death") {
    .mc_assert_numeric_column(.data, col_T, "col_T")
  }

  if (type %in% c("term", "deferred_term", "pure_endowment", "endowment")) {
    .mc_assert_numeric_scalar(n, "n", min = 0, strict_min = TRUE)
  }

  if (type %in% c("deferred", "deferred_term")) {
    .mc_assert_numeric_scalar(h, "h", min = 0, strict_min = TRUE)
  } else {
    .mc_assert_numeric_scalar(h, "h", min = 0)
  }

  m <- as.integer(round(m))

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

  Kx <- .data[[col_K]]

  if (any(Kx < 0, na.rm = TRUE)) {
    stop("`col_K` must contain non-negative simulated lifetimes.", call. = FALSE)
  }

  Tx <- if (col_T %in% names(.data) && is.numeric(.data[[col_T]])) {
    .data[[col_T]]
  } else {
    rep(NA_real_, length(Kx))
  }

  if (timing == "moment_of_death") {
    if (any(Tx < 0, na.rm = TRUE)) {
      stop("`col_T` must contain non-negative simulated lifetimes.", call. = FALSE)
    }

    if (all(is.na(Tx))) {
      stop(
        "`col_T` contains only missing values. Benefits payable at the moment ",
        "of death require complete future lifetimes.",
        call. = FALSE
      )
    }
  }

  death_time <- if (timing == "end_of_year") {
    Kx + 1
  } else {
    Tx
  }

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
      benefit_indicator <- as.numeric(Kx >= h & Kx < h + n)
    } else {
      benefit_indicator <- as.numeric(Tx > h & Tx <= h + n)
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

  pv_benefit <- ifelse(
    benefit_indicator == 1,
    benefit * v^benefit_time,
    0
  )

  n_out <- if (is.null(n)) NA_real_ else n

  .data |>
    dplyr::mutate(
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
      benefit = benefit,
      n = n_out,
      term = n_out,
      h = h,
      deferral_years = h,
      timing = timing,
      payment_timing = timing,
      benefit_time = benefit_time,
      benefit_indicator = benefit_indicator,
      "{col_pv}" := pv_benefit
    )
}
