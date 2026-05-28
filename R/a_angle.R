#' Level annuity factor a-angle-n
#'
#' Computes the actuarial present value factor for a level annuity.
#'
#' Supported timing conventions:
#' \itemize{
#'   \item \code{"immediate"}: annuity-immediate with discrete payments.
#'   \item \code{"due"}: annuity-due with discrete payments.
#'   \item \code{"continuous"}: continuous annuity.
#' }
#'
#' For discrete annuities, \code{payments_per_year = k} means payments are made
#' every \eqn{1/k} year. The function returns the annuity factor, assuming a
#' unit payment at each payment time.
#'
#' Deferral is supported through \code{deferral_years = h}. For discrete
#' annuities, the deferral must align with the payment grid, that is,
#' \eqn{h k} must be an integer.
#'
#' If \code{perpetuity = TRUE}, the infinite-term annuity factor is returned.
#'
#' @param n_years Numeric vector of payment durations in years.
#'   Ignored when \code{perpetuity = TRUE}. If \code{perpetuity = FALSE},
#'   each value must be positive and finite.
#' @param payments_per_year Positive integer vector giving the number of
#'   discrete payments per year. Ignored for continuous annuities.
#' @param rate Numeric vector of rate values.
#' @param rate_type Character vector indicating the rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the compounding frequency
#'   for nominal rates. Ignored for \code{"effective"} and \code{"force"}.
#' @param deferral_years Numeric vector of deferral times in years.
#'   Must be greater than or equal to 0.
#' @param timing Character vector. One of \code{"immediate"},
#'   \code{"due"}, or \code{"continuous"}.
#' @param perpetuity Logical vector. If \code{TRUE}, computes the perpetuity
#'   factor.
#' @param payment Numeric vector of level payment amounts. Used only when
#'   \code{output = "table"} to report the corresponding present value.
#'   The annuity factor itself is always computed for unit payments.
#' @param output Character string. Use \code{"value"} to return a numeric
#'   annuity factor, or \code{"table"} to return a tibble with intermediate
#'   calculations.
#'
#' @return
#' If \code{output = "value"}, a numeric vector of annuity factors.
#'
#' If \code{output = "table"}, a tibble with input values, equivalent rates,
#' annuity factors, payment amounts, and present values.
#'
#' @details
#' The function first converts the supplied rate to the equivalent annual
#' effective interest rate using \code{\link{standardize_interest}}.
#'
#' For finite discrete annuities:
#' \deqn{a_{\overline{n|}} = \frac{1 - v^n}{i}}
#'
#' For due annuities:
#' \deqn{\ddot{a}_{\overline{n|}} = (1+i)a_{\overline{n|}}}
#'
#' For continuous annuities:
#' \deqn{\bar{a}_{\overline{n|}} = \frac{1 - e^{-\delta n}}{\delta}}
#'
#' Input vectors must have length 1 or a common length. Missing values are
#' propagated.
#'
#' @seealso \code{\link{s_angle}}, \code{\link{standardize_interest}},
#'   \code{\link{present_value}}
#'
#' @family annuities
#'
#' @examples
#' # Numeric annuity factor
#' a_angle(n_years = 10, rate = 0.05)
#'
#' # Nominal interest converted monthly, with monthly payments
#' a_angle(
#'   n_years = 10,
#'   rate = 0.06,
#'   rate_type = "nominal_interest",
#'   m = 12,
#'   payments_per_year = 12
#' )
#'
#' # Continuous annuity
#' a_angle(
#'   n_years = 15,
#'   rate = 0.04,
#'   rate_type = "force",
#'   timing = "continuous"
#' )
#'
#' # Tibble output for teaching or auditing
#' a_angle(
#'   n_years = 10,
#'   rate = 0.05,
#'   payment = 1000,
#'   output = "table"
#' )
#'
#' # Vectorized example
#' a_angle(
#'   n_years = c(5, 10, 20),
#'   payments_per_year = c(1, 12, 1),
#'   rate = c(0.05, 0.06, 0.04),
#'   rate_type = c("effective", "nominal_interest", "force"),
#'   m = c(1, 12, 1),
#'   deferral_years = c(0, 0, 2),
#'   timing = c("immediate", "immediate", "continuous"),
#'   perpetuity = c(FALSE, FALSE, FALSE)
#' )
#'
#' @export
a_angle <- function(
    n_years = NULL,
    payments_per_year = 1L,
    rate,
    rate_type = "effective",
    m = 1L,
    deferral_years = 0,
    timing = "immediate",
    perpetuity = FALSE,
    payment = 1,
    output = c("value", "table")
) {
  output <- match.arg(output)

  if (missing(rate)) {
    stop("`rate` must be provided.", call. = FALSE)
  }

  # --- Early type validation ---
  if (!is.numeric(rate)) {
    stop("`rate` must be a numeric vector.", call. = FALSE)
  }
  if (!is.null(n_years) && !is.numeric(n_years)) {
    stop("`n_years` must be numeric or NULL.", call. = FALSE)
  }
  if (!is.numeric(payments_per_year)) {
    stop("`payments_per_year` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(m)) {
    stop("`m` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(deferral_years)) {
    stop("`deferral_years` must be numeric.", call. = FALSE)
  }
  if (!is.character(rate_type)) {
    stop("`rate_type` must be a character vector.", call. = FALSE)
  }
  if (!is.character(timing)) {
    stop("`timing` must be a character vector.", call. = FALSE)
  }
  if (!is.logical(perpetuity)) {
    stop("`perpetuity` must be a logical vector.", call. = FALSE)
  }
  if (!is.numeric(payment)) {
    stop("`payment` must be a numeric vector.", call. = FALSE)
  }

  # --- Determine common size and validate lengths ---
  size <- max(
    if (is.null(n_years)) 0L else length(n_years),
    length(payments_per_year),
    length(rate),
    length(rate_type),
    length(m),
    length(deferral_years),
    length(timing),
    length(perpetuity),
    length(payment),
    1L
  )

  valid_size <- function(x) is.null(x) || length(x) %in% c(1L, size)

  if (!valid_size(n_years) ||
      !valid_size(payments_per_year) ||
      !valid_size(rate) ||
      !valid_size(rate_type) ||
      !valid_size(m) ||
      !valid_size(deferral_years) ||
      !valid_size(timing) ||
      !valid_size(perpetuity) ||
      !valid_size(payment)) {
    stop(
      "`n_years`, `payments_per_year`, `rate`, `rate_type`, `m`, ",
      "`deferral_years`, `timing`, `perpetuity`, and `payment` must have ",
      "length 1 or a common length.",
      call. = FALSE
    )
  }

  # --- Recycle to common size ---
  if (is.null(n_years)) {
    n_years <- rep(NA_real_, size)
  } else {
    n_years <- rep_len(n_years, size)
  }

  payments_per_year <- rep_len(payments_per_year, size)
  rate              <- rep_len(rate, size)
  rate_type         <- rep_len(rate_type, size)
  m                 <- rep_len(m, size)
  deferral_years    <- rep_len(deferral_years, size)
  timing            <- rep_len(timing, size)
  perpetuity        <- rep_len(perpetuity, size)
  payment           <- rep_len(payment, size)

  # --- Value-level validation ---
  bad_deferral <- !is.na(deferral_years) &
    (!is.finite(deferral_years) | deferral_years < 0)

  if (any(bad_deferral)) {
    stop("`deferral_years` must contain only finite values >= 0 or NA.", call. = FALSE)
  }

  if (any(is.na(perpetuity))) {
    stop("`perpetuity` must not contain NA.", call. = FALSE)
  }

  bad_payment <- !is.na(payment) & !is.finite(payment)
  if (any(bad_payment)) {
    stop("`payment` must contain only finite numeric values or NA.", call. = FALSE)
  }

  timing <- tolower(timing)
  valid_timing <- c("immediate", "due", "continuous")

  bad_timing <- !is.na(timing) & !(timing %in% valid_timing)
  if (any(bad_timing)) {
    stop("`timing` must be 'immediate', 'due', or 'continuous'.", call. = FALSE)
  }

  # --- Convert to effective annual rate ---
  i_effective <- standardize_interest(type = rate_type, rate = rate, m = m)
  delta <- log1p(i_effective)

  annuity_factor <- rep(NA_real_, size)
  eps <- 1e-12

  # --- Identify element-level NA ---
  is_na_elem <- is.na(rate) |
    is.na(i_effective) |
    is.na(deferral_years) |
    is.na(timing)

  # --- Validate n_years for finite annuities ---
  needs_n <- !is_na_elem & !perpetuity
  bad_n <- needs_n & (is.na(n_years) | !is.finite(n_years) | n_years <= 0)

  if (any(bad_n)) {
    stop(
      "When `perpetuity = FALSE`, `n_years` must be positive and finite.",
      call. = FALSE
    )
  }

  # =========================================================
  # Continuous timing
  # =========================================================
  is_cont <- !is_na_elem & timing == "continuous"

  cont_perp <- is_cont & perpetuity
  if (any(cont_perp)) {
    if (any(delta[cont_perp] <= eps)) {
      stop(
        "Continuous perpetuity requires a strictly positive force of interest.",
        call. = FALSE
      )
    }

    annuity_factor[cont_perp] <-
      exp(-delta[cont_perp] * deferral_years[cont_perp]) / delta[cont_perp]
  }

  cont_fin <- is_cont & !perpetuity
  if (any(cont_fin)) {
    d_cf <- delta[cont_fin]
    n_cf <- n_years[cont_fin]
    h_cf <- deferral_years[cont_fin]

    near_zero <- abs(d_cf) < eps

    base_cf <- ifelse(
      near_zero,
      n_cf,
      (1 - exp(-d_cf * n_cf)) / d_cf
    )

    annuity_factor[cont_fin] <- base_cf * exp(-d_cf * h_cf)
  }

  # =========================================================
  # Discrete timing: immediate / due
  # =========================================================
  is_disc <- !is_na_elem & (timing == "immediate" | timing == "due")
  disc_idx <- which(is_disc)

  i_period <- rep(NA_real_, size)
  v_period <- rep(NA_real_, size)
  n_periods <- rep(NA_real_, size)
  deferral_periods <- rep(NA_real_, size)

  if (length(disc_idx) > 0L) {
    for (j in disc_idx) {
      k <- payments_per_year[j]

      if (is.na(k) || !is.finite(k) || k < 1 || k != floor(k)) {
        stop(
          "`payments_per_year` must be a positive integer for discrete annuities.",
          call. = FALSE
        )
      }

      k <- as.integer(k)

      deferral_periods_raw <- deferral_years[j] * k
      deferral_periods[j] <- round(deferral_periods_raw)

      if (abs(deferral_periods_raw - deferral_periods[j]) > 1e-10) {
        stop(
          "For discrete annuities, `deferral_years * payments_per_year` ",
          "must be an integer.",
          call. = FALSE
        )
      }

      i_period[j] <- (1 + i_effective[j])^(1 / k) - 1
      v_period[j] <- 1 / (1 + i_period[j])

      if (perpetuity[j]) {
        if (i_period[j] <= eps) {
          stop(
            "Discrete perpetuity requires a strictly positive per-period rate.",
            call. = FALSE
          )
        }

        base <- 1 / i_period[j]

        if (timing[j] == "due") {
          base <- base * (1 + i_period[j])
        }

        annuity_factor[j] <- base * v_period[j]^deferral_periods[j]
        next
      }

      n_periods_raw <- n_years[j] * k
      n_periods[j] <- round(n_periods_raw)

      if (abs(n_periods_raw - n_periods[j]) > 1e-10) {
        stop(
          "For discrete annuities, `n_years * payments_per_year` must be an integer.",
          call. = FALSE
        )
      }

      base <- if (abs(i_period[j]) < eps) {
        n_periods[j]
      } else {
        (1 - v_period[j]^n_periods[j]) / i_period[j]
      }

      if (timing[j] == "due") {
        base <- base * (1 + i_period[j])
      }

      annuity_factor[j] <- base * v_period[j]^deferral_periods[j]
    }
  }

  if (output == "value") {
    return(annuity_factor)
  }

  tibble::tibble(
    n_years = n_years,
    payments_per_year = ifelse(
      timing == "continuous",
      NA_integer_,
      as.integer(payments_per_year)
    ),
    deferral_years = deferral_years,
    timing = timing,
    perpetuity = perpetuity,
    rate_input = rate,
    rate_type = rate_type,
    m = m,
    i_effective = i_effective,
    delta = delta,
    i_period = i_period,
    v_period = v_period,
    n_periods = n_periods,
    deferral_periods = deferral_periods,
    annuity_factor = annuity_factor,
    payment = payment,
    present_value = payment * annuity_factor
  )
}
