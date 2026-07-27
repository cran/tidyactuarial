#' Compute Monte Carlo prospective reserves for life contingencies
#'
#' Computes simulated prospective loss random variables at one or more policy
#' durations. Future benefits and future premiums are revalued from each
#' requested duration.
#'
#' @param .data A data frame or tibble containing simulated curtate future
#'   lifetimes and, when required, complete future lifetimes.
#' @param t Nonnegative numeric vector of valuation durations.
#' @param i Numeric scalar. Interest-rate input used for discounting.
#' @param P Optional nonnegative numeric scalar supplied directly as the
#'   premium input.
#' @param col_P Character scalar naming a premium column in \code{.data}.
#'   The default is \code{"P"}. When it is not supplied explicitly, modern
#'   standardized premium columns are preferred.
#' @param premium_unit Unit of \code{P} or the selected premium column:
#'   \code{"annualized"}, \code{"per_payment"}, or
#'   \code{"annuity_scale"}. With \code{"auto"},
#'   \code{premium_annualized}, \code{premium_per_payment}, and
#'   \code{premium_annuity_scale} determine their own units, while legacy
#'   \code{P} and \code{premium} retain their historical annuity-scale
#'   interpretation.
#' @param benefit Nonnegative numeric scalar. Insurance benefit amount.
#' @param payment Positive numeric scalar. Amount used at each premium-annuity
#'   payment date. For the normalized \eqn{k}-thly convention, use
#'   \code{payment = 1 / k}.
#' @param k Positive integer number of premium payments per year.
#' @param type Insurance type. Supported values are \code{"whole"},
#'   \code{"term"}, \code{"deferred"}, \code{"deferred_term"},
#'   \code{"pure_endowment"}, and \code{"endowment"}.
#' @param annuity_type Premium-annuity type. Supported values are
#'   \code{"whole"}, \code{"temporary"}, \code{"deferred"},
#'   \code{"deferred_temporary"}, \code{"certain"}, and
#'   \code{"guaranteed"}.
#' @param n Optional contract or premium term in years.
#' @param h Nonnegative deferral period in years.
#' @param n_guar Optional guaranteed premium-payment period in years.
#' @param timing Death-benefit timing: \code{"end_of_year"} or
#'   \code{"moment_of_death"}.
#' @param premium_timing Premium-annuity timing: \code{"due"} or
#'   \code{"immediate"}.
#' @param reserve_timing Whether cash flows exactly at the valuation duration
#'   are included: \code{"before_payment"} or \code{"after_payment"}.
#' @param i_type Interest-rate convention.
#' @param m Positive integer nominal conversion frequency.
#' @param col_K Curtate future lifetime column.
#' @param col_T Complete future lifetime column.
#' @param in_force_basis Basis used to determine whether a simulated policy is
#'   in force: \code{"auto"}, \code{"complete"}, or \code{"curtate"}.
#' @param not_in_force Output for known scenarios that are not in force:
#'   \code{"na"} or \code{"zero"}. Missing lifetime information remains
#'   missing under either option.
#' @param col_L Output column containing the prospective loss.
#' @param tol Nonnegative numeric tolerance for consistency checks.
#' @param premium Deprecated explicit alias for \code{P}. It remains a
#'   formal argument because otherwise R partially matches it against
#'   \code{premium_unit} and \code{premium_timing}.
#' @param ... Transitional compatibility for older argument names.
#'
#' @details
#' Let \eqn{Z_t} be the present value at duration \eqn{t} of future benefits.
#' Suppose the future premium annuity is constructed with amount \eqn{c} at
#' each of the \eqn{k} payment dates per year:
#' \deqn{
#' Y_{t,c}=c\sum_{j:t_j\ge t}v^{t_j-t}.
#' }
#'
#' If \eqn{P^{(k)}} is the annualized premium, then the installment is
#' \eqn{P^{(k)}/k} and the coefficient multiplying \eqn{Y_{t,c}} is
#' \deqn{
#' q=\frac{P^{(k)}}{kc}.
#' }
#' Therefore the simulated prospective loss is
#' \deqn{
#' L_t=Z_t-qY_{t,c}.
#' }
#'
#' Under the recommended normalization \code{payment = 1 / k},
#' \eqn{q=P^{(k)}} and
#' \deqn{
#' L_t=Z_t-P^{(k)}Y_t^{(k)}.
#' }
#'
#' If no premium is supplied, the function estimates the issue scale
#' coefficient from complete benefit-annuity simulation pairs:
#' \deqn{
#' \widehat q=\frac{\overline Z_0}{\overline Y_{0,c}}.
#' }
#' This makes the sample mean issue loss equal to zero, up to numerical error,
#' when \code{reserve_timing = "before_payment"}.
#'
#' The column \code{Y_t} is the present value of the premium-annuity basis,
#' including the amount supplied through \code{payment}. The column
#' \code{future_pv_premiums} is the actual simulated present value of premium
#' income, equal to \code{premium_annuity_scale * Y_t}.
#'
#' @return A tibble with one row per simulation and requested duration. Important
#'   standardized columns include:
#' \describe{
#'   \item{Z_t}{Future benefit present value at duration \code{t}.}
#'   \item{Y_t}{Future premium-annuity basis present value.}
#'   \item{future_pv_premiums}{Actual future premium present value.}
#'   \item{premium_annualized}{Annualized premium \eqn{P^{(k)}}.}
#'   \item{premium_per_payment}{Premium paid at each date.}
#'   \item{premium_annuity_scale}{Coefficient multiplying \code{Y_t}.}
#'   \item{L_t}{Prospective loss, or the name supplied through \code{col_L}.}
#' }
#'
#' @seealso \code{\link{mc_annuity}}, \code{\link{mc_insurance}},
#'   \code{\link{mc_premium}}, \code{\link{mc_loss}},
#'   \code{\link{reserve_x}}, \code{\link{reserve_xy}}
#'
#' @family monte-carlo
#'
#' @examples
#' simulated <- tibble::tibble(
#'   Kx = c(2, 4),
#'   Tx = c(2.4, 4.7)
#' )
#'
#' simulated |>
#'   mc_reserve(
#'     t = c(0, 1),
#'     i = 0.05,
#'     type = "whole",
#'     annuity_type = "whole",
#'     payment = 1 / 12,
#'     k = 12
#'   )
#'
#' @export
mc_reserve <- function(
    .data = NULL,
    t = 0,
    i = NULL,
    P = NULL,
    col_P = "P",
    premium_unit = c(
      "auto",
      "annualized",
      "per_payment",
      "annuity_scale"
    ),
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
    tol = 1e-10,
    premium = NULL,
    ...
) {
  data_missing <- missing(.data)
  t_missing <- missing(t)
  P_missing <- missing(P)
  col_P_missing <- missing(col_P)
  payment_missing <- missing(payment)
  k_missing <- missing(k)
  type_missing <- missing(type)
  annuity_type_missing <- missing(annuity_type)
  h_missing <- missing(h)
  timing_missing <- missing(timing)
  i_type_missing <- missing(i_type)
  col_K_missing <- missing(col_K)
  col_T_missing <- missing(col_T)
  col_L_missing <- missing(col_L)

  dots <- list(...)

  dot_has <- function(name) {
    name %in% names(dots)
  }

  dot_get <- function(name) {
    dots[[name]]
  }

  allowed_old <- c(
    "data",
    "duration",
    "rate",
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

  if (dot_has("data")) {
    if (!data_missing) {
      stop("Provide only one of `.data` or deprecated `data`.", call. = FALSE)
    }

    .data <- dot_get("data")
  }

  if (dot_has("duration")) {
    if (!t_missing) {
      stop("Provide only one of `t` or deprecated `duration`.", call. = FALSE)
    }

    t <- dot_get("duration")
  }

  if (dot_has("rate")) {
    if (!is.null(i)) {
      stop("Provide only one of `i` or deprecated `rate`.", call. = FALSE)
    }

    i <- dot_get("rate")
  }

  if (!is.null(premium)) {
    if (!P_missing) {
      stop(
        "Provide only one of `P` or deprecated `premium`.",
        call. = FALSE
      )
    }

    P <- premium
    P_missing <- FALSE
  }

  if (dot_has("premium_col")) {
    if (!col_P_missing) {
      stop(
        "Provide only one of `col_P` or deprecated `premium_col`.",
        call. = FALSE
      )
    }

    col_P <- dot_get("premium_col")
    col_P_missing <- FALSE
  }

  if (dot_has("payments_per_year")) {
    if (!k_missing) {
      stop(
        "Provide only one of `k` or deprecated `payments_per_year`.",
        call. = FALSE
      )
    }

    k <- dot_get("payments_per_year")
    k_missing <- FALSE
  }

  if (dot_has("insurance")) {
    if (!type_missing) {
      stop(
        "Provide only one of `type` or deprecated `insurance`.",
        call. = FALSE
      )
    }

    type <- dot_get("insurance")
  }

  if (dot_has("annuity")) {
    if (!annuity_type_missing) {
      stop(
        "Provide only one of `annuity_type` or deprecated `annuity`.",
        call. = FALSE
      )
    }

    annuity_type <- dot_get("annuity")
  }

  if (dot_has("term")) {
    if (!is.null(n)) {
      stop("Provide only one of `n` or deprecated `term`.", call. = FALSE)
    }

    n <- dot_get("term")
  }

  if (dot_has("deferral_years")) {
    if (!h_missing) {
      stop(
        "Provide only one of `h` or deprecated `deferral_years`.",
        call. = FALSE
      )
    }

    h <- dot_get("deferral_years")
  }

  if (dot_has("guarantee_years")) {
    if (!is.null(n_guar)) {
      stop(
        "Provide only one of `n_guar` or deprecated `guarantee_years`.",
        call. = FALSE
      )
    }

    n_guar <- dot_get("guarantee_years")
  }

  if (dot_has("payment_timing")) {
    if (!timing_missing) {
      stop(
        "Provide only one of `timing` or deprecated `payment_timing`.",
        call. = FALSE
      )
    }

    timing <- dot_get("payment_timing")
  }

  if (dot_has("interest_type")) {
    if (!i_type_missing) {
      stop(
        "Provide only one of `i_type` or deprecated `interest_type`.",
        call. = FALSE
      )
    }

    i_type <- dot_get("interest_type")
  }

  if (dot_has("k_col")) {
    if (!col_K_missing) {
      stop(
        "Provide only one of `col_K` or deprecated `k_col`.",
        call. = FALSE
      )
    }

    col_K <- dot_get("k_col")
  }

  if (dot_has("tx_col")) {
    if (!col_T_missing) {
      stop(
        "Provide only one of `col_T` or deprecated `tx_col`.",
        call. = FALSE
      )
    }

    col_T <- dot_get("tx_col")
  }

  if (dot_has("reserve_col")) {
    if (!col_L_missing) {
      stop(
        "Provide only one of `col_L` or deprecated `reserve_col`.",
        call. = FALSE
      )
    }

    col_L <- dot_get("reserve_col")
  }

  premium_unit <- match.arg(premium_unit)

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

  if (!is.data.frame(.data)) {
    stop("`.data` must be a data frame or tibble.", call. = FALSE)
  }

  if (!is.numeric(t) ||
      length(t) == 0L ||
      anyNA(t) ||
      any(!is.finite(t)) ||
      any(t < 0)) {
    stop("`t` must be a nonnegative numeric vector.", call. = FALSE)
  }

  if (is.null(i)) {
    stop("`i` must be provided.", call. = FALSE)
  }

  .mc_assert_numeric_scalar(i, "i")
  .mc_assert_numeric_scalar(benefit, "benefit", min = 0)
  .mc_assert_character_scalar(col_P, "col_P")
  .mc_assert_character_scalar(col_L, "col_L")
  .mc_assert_character_scalar(col_T, "col_T")
  .mc_assert_numeric_column(.data, col_K, "col_K")
  .mc_assert_positive_integer(m, "m")

  if (!is.numeric(tol) ||
      length(tol) != 1L ||
      is.na(tol) ||
      !is.finite(tol) ||
      tol < 0) {
    stop("`tol` must be a single nonnegative finite number.", call. = FALSE)
  }

  constant_column <- function(column) {
    values <- .data[[column]]

    if (!is.numeric(values)) {
      stop("`", column, "` must be numeric.", call. = FALSE)
    }

    observed <- unique(values[!is.na(values)])

    if (length(observed) != 1L || !is.finite(observed[[1L]])) {
      stop(
        "`", column, "` must contain one constant finite value for the ",
        "simulated contract.",
        call. = FALSE
      )
    }

    observed[[1L]]
  }

  k_from_data <- NULL

  if ("payments_per_year" %in% names(.data)) {
    k_from_data <- constant_column("payments_per_year")
  }

  if ("k" %in% names(.data)) {
    k_old <- constant_column("k")

    if (!is.null(k_from_data) && abs(k_from_data - k_old) > tol) {
      stop(
        "`payments_per_year` and `k` are inconsistent in `.data`.",
        call. = FALSE
      )
    }

    if (is.null(k_from_data)) {
      k_from_data <- k_old
    }
  }

  if (k_missing && !is.null(k_from_data)) {
    k <- k_from_data
  } else if (!is.null(k_from_data) &&
             is.numeric(k) &&
             length(k) == 1L &&
             is.finite(k) &&
             abs(k - k_from_data) > tol) {
    stop(
      "`k` does not agree with the payment frequency stored in `.data`.",
      call. = FALSE
    )
  }

  .mc_assert_positive_integer(k, "k")
  k <- as.integer(round(k))
  m <- as.integer(round(m))

  payment_from_data <- NULL

  if ("payment_per_payment" %in% names(.data)) {
    payment_from_data <- constant_column("payment_per_payment")
  } else if ("payment" %in% names(.data)) {
    payment_from_data <- constant_column("payment")
  } else if ("payment_annualized" %in% names(.data)) {
    payment_from_data <- constant_column("payment_annualized") / k
  }

  if (payment_missing && !is.null(payment_from_data)) {
    payment <- payment_from_data
  } else if (!is.null(payment_from_data) &&
             is.numeric(payment) &&
             length(payment) == 1L &&
             is.finite(payment) &&
             abs(payment - payment_from_data) >
               tol * max(1, abs(payment_from_data))) {
    stop(
      "`payment` does not agree with the annuity payment amount stored ",
      "in `.data`.",
      call. = FALSE
    )
  }

  .mc_assert_numeric_scalar(
    payment,
    "payment",
    min = 0,
    strict_min = TRUE
  )

  if ("payment_annualized" %in% names(.data)) {
    payment_annualized_data <- constant_column("payment_annualized")

    if (abs(payment_annualized_data - k * payment) >
        tol * max(1, abs(payment_annualized_data))) {
      stop(
        "`payment_annualized` is inconsistent with `k * payment`.",
        call. = FALSE
      )
    }
  }

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

  if (timing == "end_of_year") {
    assert_whole_years <- function(value, name) {
      if (!is.null(value) &&
          abs(value - round(value)) > sqrt(.Machine$double.eps)) {
        stop(
          "`", name, "` must be a whole number of years when ",
          "`timing = \"end_of_year\"`.",
          call. = FALSE
        )
      }
    }

    if (type %in% c("term", "deferred_term", "pure_endowment", "endowment")) {
      assert_whole_years(n, "n")
      n <- round(n)
    }

    if (type %in% c("deferred", "deferred_term")) {
      assert_whole_years(h, "h")
      h <- round(h)
    }
  }

  assert_payment_grid <- function(value, name) {
    if (!is.null(value) &&
        abs(value * k - round(value * k)) >
          sqrt(.Machine$double.eps)) {
      stop(
        "`", name, " * k` must be an integer so that the premium ",
        "schedule contains a whole number of payments.",
        call. = FALSE
      )
    }
  }

  if (annuity_type %in% c("temporary", "deferred_temporary", "certain")) {
    assert_payment_grid(n, "n")
    n <- round(n * k) / k
  }

  if (annuity_type == "guaranteed") {
    assert_payment_grid(n_guar, "n_guar")
    n_guar <- round(n_guar * k) / k
  }

  Kx <- .data[[col_K]]

  if (any(Kx < 0, na.rm = TRUE)) {
    stop("`col_K` must contain nonnegative simulated lifetimes.", call. = FALSE)
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
      "Fractional reserve durations require `in_force_basis = 'complete'` ",
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

  if (has_T && any(Tx < 0, na.rm = TRUE)) {
    stop("`col_T` must contain nonnegative simulated lifetimes.", call. = FALSE)
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
    ifelse(
      benefit_info$benefit_indicator == 0,
      0,
      NA_real_
    )
  )

  issue_pv_annuity <- vapply(
    premium_payment_times,
    function(times) {
      payment * sum(v^times)
    },
    numeric(1)
  )

  premium_column <- NULL
  premium_input <- NULL
  premium_was_estimated <- FALSE

  if (!is.null(P)) {
    .mc_assert_numeric_scalar(P, "P", min = 0)
    premium_input <- rep(P, length(Kx))
    premium_column <- "P"
  } else {
    if (!col_P_missing) {
      if (!col_P %in% names(.data)) {
        stop(
          "The premium column selected by `col_P` was not found in `.data`.",
          call. = FALSE
        )
      }

      premium_column <- col_P
    } else {
      candidates <- c(
        "premium_annualized",
        "premium_per_payment",
        "premium_annuity_scale",
        "P",
        "premium"
      )

      available <- candidates[candidates %in% names(.data)]

      if (length(available) > 0L) {
        premium_column <- available[[1L]]
      }
    }

    if (!is.null(premium_column)) {
      if (!is.numeric(.data[[premium_column]])) {
        stop("The selected premium column must be numeric.", call. = FALSE)
      }

      premium_input <- .data[[premium_column]]

      if (any(
        is.na(premium_input) |
          !is.finite(premium_input) |
          premium_input < 0
      )) {
        stop(
          "The selected premium column must contain nonnegative finite values.",
          call. = FALSE
        )
      }
    } else {
      keep_issue <- !is.na(issue_pv_benefit) &
        !is.na(issue_pv_annuity)

      if (!any(keep_issue)) {
        stop(
          "No complete benefit-annuity simulation pairs are available to ",
          "estimate the issue premium.",
          call. = FALSE
        )
      }

      mean_Z0 <- mean(issue_pv_benefit[keep_issue])
      mean_Y0 <- mean(issue_pv_annuity[keep_issue])

      if (!is.finite(mean_Y0) || mean_Y0 <= 0) {
        stop(
          "The mean simulated premium annuity present value must be positive ",
          "when estimating the net premium internally.",
          call. = FALSE
        )
      }

      premium_input <- rep(mean_Z0 / mean_Y0, length(Kx))
      premium_column <- "premium_annuity_scale"
      premium_was_estimated <- TRUE
    }
  }

  resolved_unit <- premium_unit

  if (identical(resolved_unit, "auto")) {
    if (identical(premium_column, "premium_annualized")) {
      resolved_unit <- "annualized"
    } else if (identical(premium_column, "premium_per_payment")) {
      resolved_unit <- "per_payment"
    } else {
      resolved_unit <- "annuity_scale"
    }
  }

  if (identical(resolved_unit, "annualized")) {
    premium_annualized <- premium_input
    premium_per_payment <- premium_annualized / k
    premium_annuity_scale <- premium_per_payment / payment
  } else if (identical(resolved_unit, "per_payment")) {
    premium_per_payment <- premium_input
    premium_annualized <- k * premium_per_payment
    premium_annuity_scale <- premium_per_payment / payment
  } else {
    premium_annuity_scale <- premium_input
    premium_per_payment <- premium_annuity_scale * payment
    premium_annualized <- k * premium_per_payment
  }

  check_premium_column <- function(column, expected) {
    if (column %in% names(.data)) {
      observed <- .data[[column]]

      if (!is.numeric(observed) ||
          anyNA(observed) ||
          any(!is.finite(observed)) ||
          any(
            abs(observed - expected) >
              tol * pmax(1, abs(expected))
          )) {
        stop(
          "`", column, "` is inconsistent with the resolved premium units.",
          call. = FALSE
        )
      }
    }
  }

  check_premium_column(
    "premium_annualized",
    premium_annualized
  )
  check_premium_column(
    "premium_per_payment",
    premium_per_payment
  )
  check_premium_column(
    "premium_annuity_scale",
    premium_annuity_scale
  )

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

  premium_annualized_rep <- premium_annualized[row_rep]
  premium_per_payment_rep <- premium_per_payment[row_rep]
  premium_annuity_scale_rep <- premium_annuity_scale[row_rep]

  in_force <- if (in_force_basis == "complete") {
    ifelse(
      is.na(Tx_rep),
      NA,
      Tx_rep > t_rep
    )
  } else {
    ifelse(
      is.na(Kx_rep),
      NA,
      Kx_rep >= t_rep
    )
  }

  benefit_time_rep <- benefit_info$benefit_time[row_rep]
  benefit_indicator_rep <- benefit_info$benefit_indicator[row_rep]

  include_benefit <- if (reserve_timing == "before_payment") {
    benefit_time_rep >= t_rep
  } else {
    benefit_time_rep > t_rep
  }

  Z_t <- ifelse(
    is.na(benefit_indicator_rep) | is.na(include_benefit),
    NA_real_,
    ifelse(
      benefit_indicator_rep == 1 & include_benefit,
      benefit * v^(benefit_time_rep - t_rep),
      0
    )
  )

  Y_t <- vapply(
    seq_along(row_rep),
    function(idx) {
      times <- premium_payment_times[[row_rep[[idx]]]]

      if (length(times) == 0L &&
          (is.na(Kx_rep[[idx]]) ||
           (k > 1L && annuity_type != "certain" &&
            is.na(Tx_rep[[idx]])))) {
        return(NA_real_)
      }

      future_times <- if (reserve_timing == "before_payment") {
        times[times >= t_rep[[idx]]]
      } else {
        times[times > t_rep[[idx]]]
      }

      payment * sum(v^(future_times - t_rep[[idx]]))
    },
    numeric(1)
  )

  known_not_in_force <- !is.na(in_force) & !in_force
  unknown_in_force <- is.na(in_force)

  if (not_in_force == "na") {
    Z_t[known_not_in_force] <- NA_real_
    Y_t[known_not_in_force] <- NA_real_
  } else {
    Z_t[known_not_in_force] <- 0
    Y_t[known_not_in_force] <- 0
  }

  Z_t[unknown_in_force] <- NA_real_
  Y_t[unknown_in_force] <- NA_real_

  future_pv_premiums <- premium_annuity_scale_rep * Y_t
  L_t <- Z_t - future_pv_premiums

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
      annuity = ifelse(
        annuity_type == "whole",
        "whole_life",
        annuity_type
      ),
      benefit = benefit,
      payment = payment,
      payment_per_payment = payment,
      payment_annualized = k * payment,
      k = k,
      payments_per_year = k,
      n = if (is.null(n)) NA_real_ else n,
      term = if (is.null(n)) NA_real_ else n,
      h = h,
      deferral_years = h,
      n_guar = if (is.null(n_guar)) NA_real_ else n_guar,
      guarantee_years = if (is.null(n_guar)) NA_real_ else n_guar,
      timing = timing,
      payment_timing = premium_timing,
      premium_timing = premium_timing,
      Z_t = Z_t,
      Y_t = Y_t,
      future_pv_benefit = Z_t,
      future_pv_premium_annuity = Y_t,
      future_pv_premiums = future_pv_premiums,
      premium_annualized = premium_annualized_rep,
      premium_per_payment = premium_per_payment_rep,
      premium_annuity_scale = premium_annuity_scale_rep,
      premium_was_estimated = premium_was_estimated,
      P = premium_annuity_scale_rep,
      premium = premium_annuity_scale_rep,
      "{col_L}" := L_t
    )

  if (!identical(col_L, "reserve_loss")) {
    out <- out |>
      dplyr::mutate(
        reserve_loss = L_t
      )
  }

  out
}


#' Internal helper: benefit payment information for reserve calculations
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
  n_sim <- length(Kx)
  benefit_indicator <- rep(0, n_sim)
  benefit_time <- rep(NA_real_, n_sim)

  available <- if (timing == "end_of_year") {
    !is.na(Kx) & is.finite(Kx)
  } else {
    !is.na(Tx) & is.finite(Tx)
  }

  benefit_indicator[!available] <- NA_real_

  if (type == "whole") {
    benefit_indicator[available] <- 1
    benefit_time[available] <- death_time[available]
  }

  if (type == "term") {
    eligible <- if (timing == "end_of_year") {
      Kx < n
    } else {
      Tx <= n
    }

    benefit_indicator[available] <- as.numeric(eligible[available])
    benefit_time <- ifelse(
      benefit_indicator == 1,
      death_time,
      NA_real_
    )
  }

  if (type == "deferred") {
    eligible <- if (timing == "end_of_year") {
      Kx >= h
    } else {
      Tx > h
    }

    benefit_indicator[available] <- as.numeric(eligible[available])
    benefit_time <- ifelse(
      benefit_indicator == 1,
      death_time,
      NA_real_
    )
  }

  if (type == "deferred_term") {
    eligible <- if (timing == "end_of_year") {
      Kx >= h & Kx < h + n
    } else {
      Tx > h & Tx <= h + n
    }

    benefit_indicator[available] <- as.numeric(eligible[available])
    benefit_time <- ifelse(
      benefit_indicator == 1,
      death_time,
      NA_real_
    )
  }

  if (type == "pure_endowment") {
    eligible <- if (timing == "end_of_year") {
      Kx >= n
    } else {
      Tx >= n
    }

    benefit_indicator[available] <- as.numeric(eligible[available])
    benefit_time <- ifelse(
      benefit_indicator == 1,
      n,
      NA_real_
    )
  }

  if (type == "endowment") {
    death_before_term <- if (timing == "end_of_year") {
      Kx < n
    } else {
      Tx <= n
    }

    benefit_indicator[available] <- 1
    benefit_time[available] <- ifelse(
      death_before_term[available],
      death_time[available],
      n
    )
  }

  list(
    benefit_indicator = benefit_indicator,
    benefit_time = benefit_time
  )
}
