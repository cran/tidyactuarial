#' Geometric annuity factor
#'
#' Computes the present value or accumulated value factor for a geometric
#' annuity.
#'
#' This function replaces the more specific geometric present value and
#' accumulated value functions.
#'
#' Supported timing conventions:
#' \itemize{
#'   \item \code{"immediate"}: payments at the end of each period.
#'   \item \code{"due"}: payments at the beginning of each period.
#' }
#'
#' @param n_years Numeric vector of payment durations in years.
#'   Ignored only when \code{perpetuity = TRUE} and
#'   \code{valuation = "present"}. Otherwise, each value must be positive
#'   and finite.
#' @param payments_per_year Positive integer vector giving the number of
#'   payments per year.
#' @param rate Numeric vector of annual interest-rate values.
#' @param rate_type Character vector indicating the interest-rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the compounding frequency for
#'   nominal interest-rate inputs.
#' @param growth_rate Numeric vector of annual growth-rate values for the
#'   payments.
#' @param growth_rate_type Character vector indicating the growth-rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param growth_m Positive integer vector giving the compounding frequency for
#'   nominal growth-rate inputs.
#' @param deferral_years Numeric vector of deferral times in years.
#'   Must be greater than or equal to 0. For present values, the deferral
#'   discounts the payment block. For accumulated values, it is recorded but
#'   does not change the factor under the adopted terminal-horizon convention.
#' @param timing Character vector. One of \code{"immediate"} or \code{"due"}.
#' @param first_payment Numeric vector giving the first payment of the
#'   geometric sequence.
#' @param perpetuity Logical vector. If \code{TRUE}, computes the geometric
#'   perpetuity present value. Perpetuities are not supported for
#'   \code{valuation = "accumulated"}.
#' @param valuation Character string. Use \code{"present"} for present value
#'   or \code{"accumulated"} for accumulated value at the end of the term.
#' @param output Character string. Use \code{"value"} to return a numeric
#'   value, or \code{"table"} to return a tibble with intermediate quantities.
#'
#' @return
#' If \code{output = "value"}, a numeric vector.
#'
#' If \code{output = "table"}, a tibble with input values, equivalent rates,
#' period quantities, and both present and accumulated value factors.
#'
#' @details
#' Let \eqn{i_p} be the effective interest rate per payment period,
#' \eqn{g_p} the effective growth rate per payment period, and
#' \eqn{v_p = (1+i_p)^{-1}}.
#'
#' For a finite geometric annuity-immediate with \eqn{N} payment periods:
#' \deqn{
#' ga_N = \sum_{r=1}^{N}\frac{(1+g_p)^{r-1}}{(1+i_p)^r}
#' }
#'
#' If \eqn{i_p \neq g_p}, then
#' \deqn{
#' ga_N = \frac{1-\left(\frac{1+g_p}{1+i_p}\right)^N}{i_p-g_p}
#' }
#'
#' The accumulated value factor is computed at the standard terminal horizon:
#' \deqn{
#' gs_N = \sum_{r=1}^{N}(1+g_p)^{r-1}(1+i_p)^{N-r}.
#' }
#'
#' For annuities-due, the corresponding immediate factor is multiplied by
#' \eqn{1+i_p}.
#'
#' @seealso \code{\link{a_angle}}, \code{\link{s_angle}},
#'   \code{\link{annuity_arith}}, \code{\link{standardize_interest}}
#'
#' @family annuities
#'
#' @examples
#' # Present value of a geometric annuity
#' annuity_geom(
#'   n_years = 10,
#'   rate = 0.05,
#'   growth_rate = 0.02,
#'   valuation = "present"
#' )
#'
#' # Accumulated value of a geometric annuity
#' annuity_geom(
#'   n_years = 10,
#'   rate = 0.05,
#'   growth_rate = 0.02,
#'   valuation = "accumulated"
#' )
#'
#' # Tibble output
#' annuity_geom(
#'   n_years = 10,
#'   rate = 0.05,
#'   growth_rate = 0.02,
#'   output = "table"
#' )
#'
#' @export
annuity_geom <- function(
    n_years = NULL,
    payments_per_year = 1L,
    rate,
    rate_type = "effective",
    m = 1L,
    growth_rate = 0,
    growth_rate_type = "effective",
    growth_m = 1L,
    deferral_years = 0,
    timing = "immediate",
    first_payment = 1,
    perpetuity = FALSE,
    valuation = c("present", "accumulated"),
    output = c("value", "table")
) {
  valuation <- match.arg(valuation)
  output <- match.arg(output)

  if (missing(rate)) {
    stop("`rate` must be provided.", call. = FALSE)
  }

  # --- Early type validation ---
  if (!is.null(n_years) && !is.numeric(n_years)) {
    stop("`n_years` must be numeric or NULL.", call. = FALSE)
  }
  if (!is.numeric(payments_per_year)) {
    stop("`payments_per_year` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(rate)) {
    stop("`rate` must be a numeric vector.", call. = FALSE)
  }
  if (!is.character(rate_type)) {
    stop("`rate_type` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(m)) {
    stop("`m` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(growth_rate)) {
    stop("`growth_rate` must be a numeric vector.", call. = FALSE)
  }
  if (!is.character(growth_rate_type)) {
    stop("`growth_rate_type` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(growth_m)) {
    stop("`growth_m` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(deferral_years)) {
    stop("`deferral_years` must be numeric.", call. = FALSE)
  }
  if (!is.character(timing)) {
    stop("`timing` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(first_payment)) {
    stop("`first_payment` must be numeric.", call. = FALSE)
  }
  if (!is.logical(perpetuity)) {
    stop("`perpetuity` must be a logical vector.", call. = FALSE)
  }

  # --- Determine common size ---
  size <- max(
    if (is.null(n_years)) 0L else length(n_years),
    length(payments_per_year),
    length(rate),
    length(rate_type),
    length(m),
    length(growth_rate),
    length(growth_rate_type),
    length(growth_m),
    length(deferral_years),
    length(timing),
    length(first_payment),
    length(perpetuity),
    1L
  )

  valid_size <- function(x) is.null(x) || length(x) %in% c(1L, size)

  if (!valid_size(n_years) ||
      !valid_size(payments_per_year) ||
      !valid_size(rate) ||
      !valid_size(rate_type) ||
      !valid_size(m) ||
      !valid_size(growth_rate) ||
      !valid_size(growth_rate_type) ||
      !valid_size(growth_m) ||
      !valid_size(deferral_years) ||
      !valid_size(timing) ||
      !valid_size(first_payment) ||
      !valid_size(perpetuity)) {
    stop(
      "`n_years`, `payments_per_year`, `rate`, `rate_type`, `m`, ",
      "`growth_rate`, `growth_rate_type`, `growth_m`, `deferral_years`, ",
      "`timing`, `first_payment`, and `perpetuity` must have length 1 ",
      "or a common length.",
      call. = FALSE
    )
  }

  # --- Recycle ---
  if (is.null(n_years)) {
    n_years <- rep(NA_real_, size)
  } else {
    n_years <- rep_len(n_years, size)
  }

  payments_per_year <- rep_len(payments_per_year, size)
  rate <- rep_len(rate, size)
  rate_type <- rep_len(rate_type, size)
  m <- rep_len(m, size)
  growth_rate <- rep_len(growth_rate, size)
  growth_rate_type <- rep_len(growth_rate_type, size)
  growth_m <- rep_len(growth_m, size)
  deferral_years <- rep_len(deferral_years, size)
  timing <- rep_len(timing, size)
  first_payment <- rep_len(first_payment, size)
  perpetuity <- rep_len(perpetuity, size)

  # --- Value-level validation ---
  if (any(is.na(perpetuity))) {
    stop("`perpetuity` must not contain NA.", call. = FALSE)
  }

  if (valuation == "accumulated" && any(perpetuity)) {
    stop(
      "`perpetuity = TRUE` is not supported when `valuation = 'accumulated'`.",
      call. = FALSE
    )
  }

  bad_n <- !is.na(n_years) & (!is.finite(n_years) | n_years <= 0)
  if (any(bad_n & !perpetuity)) {
    stop(
      "When `perpetuity = FALSE`, `n_years` must contain only finite values greater than 0.",
      call. = FALSE
    )
  }

  bad_deferral <- !is.na(deferral_years) &
    (!is.finite(deferral_years) | deferral_years < 0)

  if (any(bad_deferral)) {
    stop("`deferral_years` must contain only finite values >= 0 or NA.", call. = FALSE)
  }

  bad_first_payment <- !is.na(first_payment) & !is.finite(first_payment)
  if (any(bad_first_payment)) {
    stop("`first_payment` must contain only finite numeric values or NA.", call. = FALSE)
  }

  timing <- tolower(timing)
  valid_timing <- c("immediate", "due")

  bad_timing <- !is.na(timing) & !(timing %in% valid_timing)
  if (any(bad_timing)) {
    stop("`timing` must be 'immediate' or 'due'.", call. = FALSE)
  }

  # --- Convert rates ---
  i_effective <- standardize_interest(type = rate_type, rate = rate, m = m)
  g_effective <- standardize_interest(
    type = growth_rate_type,
    rate = growth_rate,
    m = growth_m
  )

  present_value_factor <- rep(NA_real_, size)
  accumulated_value_factor <- rep(NA_real_, size)
  selected_factor <- rep(NA_real_, size)

  i_period <- rep(NA_real_, size)
  g_period <- rep(NA_real_, size)
  v_period <- rep(NA_real_, size)
  n_periods <- rep(NA_real_, size)
  deferral_periods <- rep(NA_real_, size)

  eps <- 1e-12

  is_na_elem <- is.na(rate) |
    is.na(growth_rate) |
    is.na(i_effective) |
    is.na(g_effective) |
    is.na(deferral_years) |
    is.na(timing) |
    is.na(first_payment)

  computable_idx <- which(!is_na_elem)

  for (j in computable_idx) {
    k <- payments_per_year[j]

    if (is.na(k) || !is.finite(k) || k < 1 || k != floor(k)) {
      stop("`payments_per_year` must be a positive integer.", call. = FALSE)
    }

    k <- as.integer(k)

    deferral_raw <- deferral_years[j] * k
    h <- round(deferral_raw)

    if (abs(deferral_raw - h) > 1e-10) {
      stop(
        "For discrete annuities, `deferral_years * payments_per_year` must be an integer.",
        call. = FALSE
      )
    }

    i_period[j] <- (1 + i_effective[j])^(1 / k) - 1
    g_period[j] <- (1 + g_effective[j])^(1 / k) - 1
    v_period[j] <- 1 / (1 + i_period[j])
    deferral_periods[j] <- h

    if (perpetuity[j]) {
      if (i_period[j] <= g_period[j] + eps) {
        stop(
          "For a geometric perpetuity, the per-period interest rate must be strictly greater than the per-period growth rate.",
          call. = FALSE
        )
      }

      base_pv <- 1 / (i_period[j] - g_period[j])

      if (timing[j] == "due") {
        base_pv <- (1 + i_period[j]) * base_pv
      }

      present_value_factor[j] <- first_payment[j] * base_pv * v_period[j]^h
      selected_factor[j] <- present_value_factor[j]
      next
    }

    n_raw <- n_years[j] * k
    n <- round(n_raw)

    if (abs(n_raw - n) > 1e-10) {
      stop(
        "For discrete annuities, `n_years * payments_per_year` must be an integer.",
        call. = FALSE
      )
    }

    n_periods[j] <- n

    if (abs(i_period[j] - g_period[j]) < eps) {
      base_pv <- n / (1 + i_period[j])
      base_av <- n * (1 + i_period[j])^(n - 1)
    } else {
      q <- (1 + g_period[j]) / (1 + i_period[j])

      base_pv <- (1 - q^n) / (i_period[j] - g_period[j])

      base_av <- ((1 + i_period[j])^n - (1 + g_period[j])^n) /
        (i_period[j] - g_period[j])
    }

    if (timing[j] == "due") {
      base_pv <- (1 + i_period[j]) * base_pv
      base_av <- (1 + i_period[j]) * base_av
    }

    present_value_factor[j] <- first_payment[j] * base_pv * v_period[j]^h
    accumulated_value_factor[j] <- first_payment[j] * base_av

    selected_factor[j] <- if (valuation == "present") {
      present_value_factor[j]
    } else {
      accumulated_value_factor[j]
    }
  }

  if (output == "value") {
    return(selected_factor)
  }

  tibble::tibble(
    n_years = n_years,
    payments_per_year = as.integer(payments_per_year),
    deferral_years = deferral_years,
    timing = timing,
    perpetuity = perpetuity,
    valuation = valuation,
    rate_input = rate,
    rate_type = rate_type,
    m = m,
    growth_rate_input = growth_rate,
    growth_rate_type = growth_rate_type,
    growth_m = growth_m,
    i_effective = i_effective,
    g_effective = g_effective,
    i_period = i_period,
    g_period = g_period,
    v_period = v_period,
    n_periods = n_periods,
    deferral_periods = deferral_periods,
    first_payment = first_payment,
    present_value_factor = present_value_factor,
    accumulated_value_factor = accumulated_value_factor,
    annuity_value = selected_factor
  )
}
