#' Compute simulated present values for life annuities
#'
#' Computes Monte Carlo simulated present values of life annuity payments from
#' simulated future lifetimes, using compact actuarial notation.
#'
#' This function is designed to be used after \code{\link{simulate_lifetime}},
#' \code{\link{simulate_lifetimes}}, or \code{\link{mc_multilife_status}}. It
#' takes simulated values of the curtate future lifetime \eqn{K_x}, and when
#' needed the complete future lifetime \eqn{T_x}, and evaluates the present
#' value random variable associated with classical annuity benefits.
#'
#' @param .data A data frame or tibble containing simulated future lifetimes,
#'   typically returned by \code{\link{simulate_lifetime}} or by
#'   \code{\link{mc_multilife_status}} when working with multiple-life statuses.
#' @param i Numeric scalar. Interest-rate input used for discounting.
#' @param payment Numeric scalar. Amount of each annuity payment. Default is
#'   \code{1}.
#' @param k Positive integer. Number of annuity payments per year. Default is
#'   \code{1}, corresponding to annual payments.
#' @param type Character string specifying the annuity type. Canonical options
#'   are \code{"whole"}, \code{"temporary"}, \code{"deferred"},
#'   \code{"deferred_temporary"}, \code{"certain"}, and
#'   \code{"guaranteed"}. The transitional alias \code{"whole_life"} is also
#'   accepted and mapped to \code{"whole"}.
#' @param n Numeric scalar. Term of the annuity in years. Required for
#'   \code{"temporary"}, \code{"deferred_temporary"}, and \code{"certain"}
#'   annuities.
#' @param h Numeric scalar. Deferral period in years. Default is \code{0}.
#'   Required to be positive for \code{"deferred"} and
#'   \code{"deferred_temporary"} annuities.
#' @param n_guar Numeric scalar. Guaranteed payment period in years. Required
#'   for \code{type = "guaranteed"}.
#' @param timing Character string specifying the annuity payment timing.
#'   Available options are \code{"immediate"} and \code{"due"}.
#' @param i_type Character string specifying the interest-rate convention.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}. The transitional value
#'   \code{"nominal"} is accepted and treated as \code{"nominal_interest"}.
#' @param m Positive integer. Number of interest conversion periods per year
#'   for nominal annual rates. Default is \code{1}. This argument controls the
#'   interest-rate conversion frequency only. It does not represent annuity
#'   payment frequency.
#' @param col_K Character string. Name of the column containing simulated
#'   curtate future lifetimes. Default is \code{"Kx"}.
#' @param col_T Character string. Name of the column containing simulated
#'   complete future lifetimes. Default is \code{"Tx"}. This column is required
#'   when \code{k > 1}, except for annuities certain.
#' @param col_pv Character string. Name of the output column containing
#'   simulated present values of annuity payments. Default is
#'   \code{"pv_annuity"}.
#' @param ... Transitional compatibility for older calls using \code{data},
#'   \code{rate}, \code{payments_per_year}, \code{annuity}, \code{term},
#'   \code{deferral_years}, \code{guarantee_years}, \code{interest_type},
#'   \code{k_col}, \code{tx_col}, and \code{annuity_col}.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{i} is the interest-rate input, \code{i_type}
#' is the interest-rate convention, \code{m} is the nominal conversion
#' frequency, \code{k} is the annuity payment frequency, \code{n} is the
#' annuity term, and \code{h} is the deferral period.
#'
#' The arguments \code{m} and \code{k} have different meanings:
#' \itemize{
#'   \item \code{m} is used only for nominal interest-rate conversion.
#'   \item \code{k} controls how frequently annuity payments are made.
#' }
#'
#' The argument \code{payment} represents the amount of each annuity payment.
#' Thus, for a monthly annuity with total annual payment equal to 1, use
#' \code{payment = 1 / 12} and \code{k = 12}.
#'
#' For annual payments, \code{k = 1}, the function works directly with
#' \eqn{K_x}. For fractional payments, such as monthly, quarterly, or
#' semiannual payments, the function uses \eqn{T_x} to determine whether the
#' life is alive at each fractional payment time.
#'
#' For annual whole-life annuity-immediate, the simulated present value is
#' \deqn{
#'   Y = \sum_{j=1}^{K_x} c v^j,
#' }
#' where \eqn{c} is the amount of each payment.
#'
#' For annual whole-life annuity-due, the simulated present value is
#' \deqn{
#'   \ddot{Y} = \sum_{j=0}^{K_x} c v^j.
#' }
#'
#' The function returns simulated present values, not only their expected value.
#' Therefore the resulting column can be summarized with \code{\link{summary_mc}},
#' plotted with \code{ggplot2}, or used to construct premiums, losses, and
#' reserves.
#'
#' @return A tibble with the original simulation columns and additional columns:
#' \describe{
#'   \item{i}{Original interest-rate input.}
#'   \item{i_type}{Interest-rate convention.}
#'   \item{m}{Interest conversion frequency.}
#'   \item{i_effective}{Equivalent annual effective interest rate.}
#'   \item{v}{Annual discount factor.}
#'   \item{type}{Canonical annuity type.}
#'   \item{payment}{Amount of each annuity payment.}
#'   \item{k}{Annuity payment frequency.}
#'   \item{n}{Annuity term, if applicable.}
#'   \item{h}{Deferral period.}
#'   \item{n_guar}{Guaranteed period, if applicable.}
#'   \item{timing}{Annuity payment timing.}
#'   \item{n_payments}{Number of payments made in the simulated scenario.}
#'   \item{first_payment_time}{First payment time in the simulated scenario.}
#'   \item{last_payment_time}{Last payment time in the simulated scenario.}
#'   \item{pv_annuity}{Simulated present value of annuity payments, or another
#'   name supplied through \code{col_pv}.}
#' }
#'
#' For transition, the output also includes legacy columns such as
#' \code{rate}, \code{interest_type}, \code{effective_rate},
#' \code{discount_factor}, \code{annuity}, \code{payments_per_year},
#' \code{term}, \code{deferral_years}, and \code{guarantee_years}.
#'
#' @seealso
#' \code{\link{simulate_lifetime}}, \code{\link{simulate_lifetimes}},
#' \code{\link{mc_multilife_status}}, \code{\link{mc_insurance}},
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
#' # Annual whole-life annuity-due
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     seed = 123
#'   ) |>
#'   mc_annuity(
#'     i = 0.05,
#'     type = "whole",
#'     payment = 1,
#'     k = 1,
#'     timing = "due"
#'   )
#'
#' # Monthly whole-life annuity-due with total annual payment equal to 1
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     frac = "udd",
#'     seed = 123
#'   ) |>
#'   mc_annuity(
#'     i = 0.05,
#'     type = "whole",
#'     payment = 1 / 12,
#'     k = 12,
#'     timing = "due"
#'   )
#'
#' # Quarterly temporary life annuity-immediate
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     frac = "udd",
#'     seed = 123
#'   ) |>
#'   mc_annuity(
#'     i = 0.05,
#'     type = "temporary",
#'     n = 20,
#'     payment = 1 / 4,
#'     k = 4,
#'     timing = "immediate"
#'   )
#'
#' # Nominal interest convertible monthly, with quarterly payments
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     frac = "udd",
#'     seed = 123
#'   ) |>
#'   mc_annuity(
#'     i = 0.06,
#'     i_type = "nominal_interest",
#'     m = 12,
#'     type = "whole",
#'     payment = 1 / 4,
#'     k = 4,
#'     timing = "due"
#'   )
#'
#' # Multiple-life status workflow
#' lt |>
#'   simulate_lifetimes(
#'     x = c(60, 58),
#'     n_sim = 25,
#'     frac = "udd",
#'     seed = 123
#'   ) |>
#'   mc_multilife_status(status = "joint") |>
#'   mc_annuity(
#'     i = 0.04,
#'     type = "whole",
#'     payment = 1,
#'     k = 1,
#'     timing = "due",
#'     col_K = "K_status",
#'     col_T = "T_status"
#'   )
#'
#' @export
mc_annuity <- function(
    .data = NULL,
    i,
    payment = 1,
    k = 1L,
    type = c(
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
    timing = c("immediate", "due"),
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
    col_pv = "pv_annuity",
    ...
) {
  dots <- list(...)
  type_missing <- missing(type)

  # -------------------------------------------------------------------------
  # Transitional compatibility with the previous public API
  # -------------------------------------------------------------------------

  allowed_old <- c(
    "data",
    "rate",
    "payments_per_year",
    "annuity",
    "term",
    "deferral_years",
    "guarantee_years",
    "interest_type",
    "k_col",
    "tx_col",
    "annuity_col"
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
    if (!missing(i)) {
      stop("Provide only one of `i` or deprecated `rate`.", call. = FALSE)
    }

    i <- dots$rate
  }

  if (!is.null(dots$payments_per_year)) {
    if (!identical(k, 1L) && !identical(k, 1)) {
      stop(
        "Provide only one of `k` or deprecated `payments_per_year`.",
        call. = FALSE
      )
    }

    k <- dots$payments_per_year
  }

  if (!is.null(dots$annuity)) {
    if (!type_missing) {
      stop("Provide only one of `type` or deprecated `annuity`.", call. = FALSE)
    }

    type <- dots$annuity
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

  if (!is.null(dots$interest_type)) {
    if (!missing(i_type)) {
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

  if (!is.null(dots$annuity_col)) {
    if (!identical(col_pv, "pv_annuity")) {
      stop("Provide only one of `col_pv` or deprecated `annuity_col`.",
           call. = FALSE)
    }

    col_pv <- dots$annuity_col
  }

  type <- match.arg(
    type,
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

  if (identical(type, "whole_life")) {
    type <- "whole"
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

  if (missing(i)) {
    stop("`i` must be provided.", call. = FALSE)
  }

  .mc_assert_numeric_scalar(i, "i")
  .mc_assert_numeric_scalar(payment, "payment", min = 0)
  .mc_assert_positive_integer(k, "k")
  .mc_assert_numeric_scalar(m, "m", min = 0, strict_min = TRUE)
  .mc_assert_numeric_column(.data, col_K, "col_K")
  .mc_assert_character_scalar(col_T, "col_T")
  .mc_assert_character_scalar(col_pv, "col_pv")

  if (type %in% c("temporary", "deferred_temporary", "certain")) {
    .mc_assert_numeric_scalar(n, "n", min = 0, strict_min = TRUE)
  }

  if (type %in% c("deferred", "deferred_temporary")) {
    .mc_assert_numeric_scalar(h, "h", min = 0, strict_min = TRUE)
  } else {
    .mc_assert_numeric_scalar(h, "h", min = 0)
  }

  if (type == "guaranteed") {
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

  k <- as.integer(round(k))
  m <- as.integer(round(m))

  needs_T <- k > 1L && type != "certain"

  if (needs_T) {
    .mc_assert_numeric_column(.data, col_T, "col_T")

    if (all(is.na(.data[[col_T]]))) {
      stop(
        "`col_T` contains only missing values. Fractional life-contingent ",
        "payments require complete future lifetimes.",
        call. = FALSE
      )
    }
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

  Kx <- .data[[col_K]]

  if (any(Kx < 0, na.rm = TRUE)) {
    stop("`col_K` must contain non-negative simulated lifetimes.", call. = FALSE)
  }

  Tx <- if (col_T %in% names(.data) && is.numeric(.data[[col_T]])) {
    .data[[col_T]]
  } else {
    rep(NA_real_, length(Kx))
  }

  if (needs_T && any(Tx < 0, na.rm = TRUE)) {
    stop("`col_T` must contain non-negative simulated lifetimes.", call. = FALSE)
  }

  payment_interval <- 1 / k

  payment_times <- Map(
    function(K, T) {
      if (k == 1L) {
        return(
          .mc_annuity_payment_times_annual(
            K = K,
            type = type,
            n = n,
            h = h,
            n_guar = n_guar,
            timing = timing
          )
        )
      }

      .mc_annuity_payment_times_fractional(
        T = T,
        type = type,
        n = n,
        h = h,
        n_guar = n_guar,
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
      if (length(times) == 0L) {
        return(NA_real_)
      }

      min(times)
    },
    numeric(1)
  )

  last_payment_time <- vapply(
    payment_times,
    function(times) {
      if (length(times) == 0L) {
        return(NA_real_)
      }

      max(times)
    },
    numeric(1)
  )

  n_out <- if (is.null(n)) NA_real_ else n
  n_guar_out <- if (is.null(n_guar)) NA_real_ else n_guar

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
      annuity = ifelse(type == "whole", "whole_life", type),
      payment = payment,
      k = k,
      payments_per_year = k,
      n = n_out,
      term = n_out,
      h = h,
      deferral_years = h,
      n_guar = n_guar_out,
      guarantee_years = n_guar_out,
      timing = timing,
      n_payments = n_payments,
      first_payment_time = first_payment_time,
      last_payment_time = last_payment_time,
      "{col_pv}" := pv_annuity
    )
}


#' Internal helper: annual annuity payment times
#'
#' @param K Numeric scalar. Simulated curtate future lifetime.
#' @param type Character string. Canonical annuity type.
#' @param n Numeric scalar or `NULL`.
#' @param h Numeric scalar. Deferral period.
#' @param n_guar Numeric scalar or `NULL`.
#' @param timing Character string. Either `"immediate"` or `"due"`.
#'
#' @return Numeric vector with annual payment times.
#'
#' @keywords internal
.mc_annuity_payment_times_annual <- function(
    K,
    type,
    n,
    h,
    n_guar,
    timing
) {
  if (type != "certain" &&
      (is.na(K) || !is.finite(K) || K < 0)) {
    return(numeric(0))
  }

  if (type == "whole") {
    if (timing == "immediate") {
      return(.mc_payment_times(1, K))
    }

    if (timing == "due") {
      return(.mc_payment_times(0, K))
    }
  }

  if (type == "temporary") {
    if (timing == "immediate") {
      return(.mc_payment_times(1, min(K, n)))
    }

    if (timing == "due") {
      return(.mc_payment_times(0, min(K, n - 1)))
    }
  }

  if (type == "deferred") {
    if (timing == "immediate") {
      return(.mc_payment_times(h + 1, K))
    }

    if (timing == "due") {
      return(.mc_payment_times(h, K))
    }
  }

  if (type == "deferred_temporary") {
    if (timing == "immediate") {
      return(
        .mc_payment_times(
          h + 1,
          min(K, h + n)
        )
      )
    }

    if (timing == "due") {
      return(
        .mc_payment_times(
          h,
          min(K, h + n - 1)
        )
      )
    }
  }

  if (type == "certain") {
    if (timing == "immediate") {
      return(.mc_payment_times(1, n))
    }

    if (timing == "due") {
      return(.mc_payment_times(0, n - 1))
    }
  }

  if (type == "guaranteed") {
    if (timing == "immediate") {
      return(.mc_payment_times(1, max(K, n_guar)))
    }

    if (timing == "due") {
      return(.mc_payment_times(0, max(K, n_guar - 1)))
    }
  }

  numeric(0)
}


#' Internal helper: fractional annuity payment times
#'
#' @param T Numeric scalar. Simulated complete future lifetime.
#' @param type Character string. Canonical annuity type.
#' @param n Numeric scalar or `NULL`.
#' @param h Numeric scalar. Deferral period.
#' @param n_guar Numeric scalar or `NULL`.
#' @param timing Character string. Either `"immediate"` or `"due"`.
#' @param payment_interval Numeric scalar. Time between payments.
#'
#' @return Numeric vector with fractional payment times.
#'
#' @keywords internal
.mc_annuity_payment_times_fractional <- function(
    T,
    type,
    n,
    h,
    n_guar,
    timing,
    payment_interval
) {
  eps <- sqrt(.Machine$double.eps)

  if (type == "certain") {
    if (timing == "immediate") {
      return(.mc_payment_times(payment_interval, n, by = payment_interval))
    }

    if (timing == "due") {
      return(
        .mc_payment_times(
          0,
          n - payment_interval,
          by = payment_interval
        )
      )
    }
  }

  if (is.na(T) || !is.finite(T) || T < 0) {
    return(numeric(0))
  }

  last_alive_time <- if (timing == "immediate") {
    floor(T / payment_interval + eps) * payment_interval
  } else {
    floor((T - eps) / payment_interval) * payment_interval
  }

  if (last_alive_time < 0) {
    return(numeric(0))
  }

  if (type == "whole") {
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

  if (type == "temporary") {
    if (timing == "immediate") {
      last_time <- min(last_alive_time, n)

      return(
        .mc_payment_times(
          payment_interval,
          last_time,
          by = payment_interval
        )
      )
    }

    if (timing == "due") {
      last_time <- min(last_alive_time, n - payment_interval)

      return(.mc_payment_times(0, last_time, by = payment_interval))
    }
  }

  if (type == "deferred") {
    if (timing == "immediate") {
      first_time <- h + payment_interval

      return(.mc_payment_times(first_time, last_alive_time, by = payment_interval))
    }

    if (timing == "due") {
      first_time <- h

      return(.mc_payment_times(first_time, last_alive_time, by = payment_interval))
    }
  }

  if (type == "deferred_temporary") {
    if (timing == "immediate") {
      first_time <- h + payment_interval
      last_time <- min(last_alive_time, h + n)

      return(.mc_payment_times(first_time, last_time, by = payment_interval))
    }

    if (timing == "due") {
      first_time <- h
      last_time <- min(last_alive_time, h + n - payment_interval)

      return(.mc_payment_times(first_time, last_time, by = payment_interval))
    }
  }

  if (type == "guaranteed") {
    if (timing == "immediate") {
      guaranteed_last_time <- n_guar
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
      guaranteed_last_time <- n_guar - payment_interval
      last_time <- max(last_alive_time, guaranteed_last_time)

      return(.mc_payment_times(0, last_time, by = payment_interval))
    }
  }

  numeric(0)
}
