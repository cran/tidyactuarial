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
#' every \eqn{1/k} year. The function returns the annuity factor only, assuming
#' a unit payment at each payment time.
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
#' @param type Character vector indicating the rate type.
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
#'
#' @return Numeric vector of annuity factors.
#'
#' @details
#' The function first converts the supplied rate to the equivalent annual
#' effective interest rate using \code{\link{standardize_interest}}.
#'
#' For finite discrete annuities:
#' \deqn{a_{\overline{n|}} = \frac{1 - v^n}{i}}{a_n| = (1 - v^n)/(i)}
#'
#' For due annuities:
#' \deqn{\ddot{a}_{\overline{n|}} = (1+i)a_{\overline{n|}}}{a_n| = (1+i)a_n|}
#'
#' For continuous annuities:
#' \deqn{\bar{a}_{\overline{n|}} = \frac{1 - e^{-\delta n}}{\delta}}{a_n| = 1 - e^-delta ndelta}
#'
#' Input vectors must have length 1 or a common length.
#' Missing values are propagated.
#'
#' @seealso \code{\link{a_angle_tbl}}, \code{\link{s_angle}}, \code{\link{Da_angle}}, \code{\link{Ia_angle}},
#'   \code{\link{standardize_interest}}, \code{\link{present_value}}
#'
#' @family annuities
#'
#' @examples
#' # Simple scalar examples
#' a_angle(n_years = 10, rate = 0.05, type = "effective")
#' a_angle(n_years = 10, rate = 0.06, type = "nominal_interest", m = 12,
#'         payments_per_year = 12)
#' a_angle(n_years = 15, rate = 0.04, type = "force", timing = "continuous")
#'
#' # Medium vectorized example
#' a_angle(
#'   n_years = c(5, 10, 20),
#'   payments_per_year = c(1, 12, 1),
#'   rate = c(0.05, 0.06, 0.04),
#'   type = c("effective", "nominal_interest", "force"),
#'   m = c(1, 12, 1),
#'   deferral_years = c(0, 0, 2),
#'   timing = c("immediate", "immediate", "continuous"),
#'   perpetuity = c(FALSE, FALSE, FALSE)
#' )
#'
#' # Use inside a data pipeline
#' if (requireNamespace("dplyr", quietly = TRUE) &&
#'     requireNamespace("tibble", quietly = TRUE)) {
#'   contracts <- tibble::tibble(
#'     n_years = c(10, 15, 20),
#'     rate = c(0.05, 0.06, 0.04),
#'     type = c("effective", "nominal_interest", "force"),
#'     m = c(1, 12, 1),
#'     payments_per_year = c(1, 12, NA),
#'     deferral_years = c(0, 0, 3),
#'     timing = c("immediate", "immediate", "continuous"),
#'     perpetuity = c(FALSE, FALSE, FALSE)
#'   )
#'
#'   dplyr::mutate(
#'     contracts,
#'     factor = a_angle(
#'       n_years = n_years,
#'       payments_per_year = dplyr::coalesce(payments_per_year, 1L),
#'       rate = rate,
#'       type = type,
#'       m = m,
#'       deferral_years = deferral_years,
#'       timing = timing,
#'       perpetuity = perpetuity
#'     )
#'   )
#' }
#'
#' @export
a_angle <- function(
    n_years = NULL,
    payments_per_year = 1L,
    rate,
    type = "effective",
    m = 1L,
    deferral_years = 0,
    timing = "immediate",
    perpetuity = FALSE
) {
  if (missing(rate)) {
    stop("`rate` must be provided.", call. = FALSE)
  }

  # --- Early type validation (before recycling) ---
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
  if (!is.character(type)) {
    stop("`type` must be a character vector.", call. = FALSE)
  }
  if (!is.character(timing)) {
    stop("`timing` must be a character vector.", call. = FALSE)
  }
  if (!is.logical(perpetuity)) {
    stop("`perpetuity` must be a logical vector.", call. = FALSE)
  }

  # --- Determine common size and validate lengths ---
  size <- max(
    length(n_years),
    length(payments_per_year),
    length(rate),
    length(type),
    length(m),
    length(deferral_years),
    length(timing),
    length(perpetuity),
    1L
  )

  valid_size <- function(x) is.null(x) || length(x) %in% c(1L, size)

  if (!valid_size(n_years) || !valid_size(payments_per_year) || !valid_size(rate) ||
      !valid_size(type) || !valid_size(m) || !valid_size(deferral_years) ||
      !valid_size(timing) || !valid_size(perpetuity)) {
    stop(
      "`n_years`, `payments_per_year`, `rate`, `type`, `m`, ",
      "`deferral_years`, `timing`, and `perpetuity` must have length 1 ",
      "or a common length.",
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
  type              <- rep_len(type, size)
  m                 <- rep_len(m, size)
  deferral_years    <- rep_len(deferral_years, size)
  timing            <- rep_len(timing, size)
  perpetuity        <- rep_len(perpetuity, size)

  # --- Value-level validation ---
  bad_deferral <- !is.na(deferral_years) & (!is.finite(deferral_years) | deferral_years < 0)
  if (any(bad_deferral)) {
    stop("`deferral_years` must contain only finite values >= 0 or NA.", call. = FALSE)
  }

  if (any(is.na(perpetuity))) {
    stop("`perpetuity` must not contain NA.", call. = FALSE)
  }

  timing <- tolower(timing)
  valid_timing <- c("immediate", "due", "continuous")
  bad_timing <- !is.na(timing) & !(timing %in% valid_timing)
  if (any(bad_timing)) {
    stop("`timing` must be 'immediate', 'due', or 'continuous'.", call. = FALSE)
  }

  # --- Convert to effective annual rate ---
  i_effective <- standardize_interest(type = type, rate = rate, m = m)
  delta <- log1p(i_effective)

  out <- rep(NA_real_, size)
  eps <- 1e-12

  # --- Identify element-level NA (skip entirely) ---
  is_na_elem <- is.na(rate) | is.na(i_effective) | is.na(deferral_years) | is.na(timing)

  # --- Validate n_years for non-perpetuity, non-NA elements ---
  needs_n <- !is_na_elem & !perpetuity
  bad_n <- needs_n & (is.na(n_years) | !is.finite(n_years) | n_years <= 0)
  if (any(bad_n)) {
    stop(
      "When `perpetuity = FALSE`, `n_years` must be positive and finite.",
      call. = FALSE
    )
  }

  # =========================================================
  # CONTINUOUS timing (vectorized)
  # =========================================================
  is_cont <- !is_na_elem & timing == "continuous"

  # Continuous perpetuity
  cont_perp <- is_cont & perpetuity
  if (any(cont_perp)) {
    if (any(delta[cont_perp] <= eps)) {
      stop(
        "Continuous perpetuity requires a strictly positive force of interest.",
        call. = FALSE
      )
    }
    out[cont_perp] <- exp(-delta[cont_perp] * deferral_years[cont_perp]) / delta[cont_perp]
  }

  # Continuous finite
  cont_fin <- is_cont & !perpetuity
  if (any(cont_fin)) {
    d_cf <- delta[cont_fin]
    n_cf <- n_years[cont_fin]
    h_cf <- deferral_years[cont_fin]
    near_zero <- abs(d_cf) < eps
    base_cf <- ifelse(near_zero, n_cf, (1 - exp(-d_cf * n_cf)) / d_cf)
    out[cont_fin] <- base_cf * exp(-d_cf * h_cf)
  }

  # =========================================================
  # DISCRETE timing (immediate / due) - element-wise loop
  #   (alignment checks require per-element validation)
  # =========================================================
  is_disc <- !is_na_elem & (timing == "immediate" | timing == "due")
  disc_idx <- which(is_disc)

  if (length(disc_idx) > 0L) {
    for (j in disc_idx) {
      k <- payments_per_year[j]
      if (is.na(k) || !is.finite(k) || k < 1 || k != floor(k)) {
        stop("`payments_per_year` must be a positive integer for discrete annuities.",
             call. = FALSE)
      }
      k <- as.integer(k)

      deferral_periods_raw <- deferral_years[j] * k
      deferral_periods <- round(deferral_periods_raw)
      if (abs(deferral_periods_raw - deferral_periods) > 1e-10) {
        stop(
          "For discrete annuities, `deferral_years * payments_per_year` ",
          "must be an integer.",
          call. = FALSE
        )
      }

      i_period <- (1 + i_effective[j])^(1 / k) - 1
      v_period <- 1 / (1 + i_period)

      if (perpetuity[j]) {
        if (i_period <= eps) {
          stop(
            "Discrete perpetuity requires a strictly positive per-period rate.",
            call. = FALSE
          )
        }
        base <- 1 / i_period
        if (timing[j] == "due") {
          base <- base * (1 + i_period)
        }
        out[j] <- base * v_period^deferral_periods
        next
      }

      n_periods_raw <- n_years[j] * k
      n_periods <- round(n_periods_raw)
      if (abs(n_periods_raw - n_periods) > 1e-10) {
        stop(
          "For discrete annuities, `n_years * payments_per_year` must be an integer.",
          call. = FALSE
        )
      }

      base <- if (abs(i_period) < eps) {
        n_periods
      } else {
        (1 - v_period^n_periods) / i_period
      }

      if (timing[j] == "due") {
        base <- base * (1 + i_period)
      }

      out[j] <- base * v_period^deferral_periods
    }
  }

  out
}


#' Level annuity details in tibble form
#'
#' Computes the actuarial present value factor for a level annuity and returns
#' a tibble with the main input values, implied rates, annuity factor,
#' payment amount, and present value.
#'
#' This is a reporting wrapper around \code{\link{a_angle}}. The annuity factor assumes
#' unit payments. The present value is then computed as
#' \deqn{PV = R \times a}{PV = R * a}
#' where \eqn{R} is the payment amount and \eqn{a} is the annuity factor.
#'
#' @param n_years Numeric vector of payment durations in years.
#' @param payments_per_year Positive integer vector giving the number of
#'   discrete payments per year. Ignored for continuous annuities.
#' @param rate Numeric vector of rate values.
#' @param type Character vector indicating the rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the compounding frequency
#'   for nominal rates.
#' @param deferral_years Numeric vector of deferral times in years.
#' @param timing Character vector. One of \code{"immediate"},
#'   \code{"due"}, or \code{"continuous"}.
#' @param payment Numeric vector of level payment amounts.
#' @param perpetuity Logical vector. If \code{TRUE}, computes the perpetuity
#'   factor.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{n_years}{Payment duration in years.}
#'   \item{payments_per_year}{Number of payments per year.}
#'   \item{deferral_years}{Deferral period in years.}
#'   \item{timing}{Payment timing convention.}
#'   \item{perpetuity}{Whether the annuity is perpetual.}
#'   \item{rate_input}{Original supplied rate.}
#'   \item{rate_type}{Type of supplied rate.}
#'   \item{m}{Compounding frequency for nominal rates.}
#'   \item{i_effective}{Equivalent annual effective rate.}
#'   \item{delta}{Equivalent force of interest.}
#'   \item{i_period}{Equivalent per-payment effective rate for discrete annuities.}
#'   \item{v_period}{Equivalent per-payment discount factor for discrete annuities.}
#'   \item{n_periods}{Number of payment periods for finite discrete annuities.}
#'   \item{deferral_periods}{Number of deferred periods for discrete annuities.}
#'   \item{annuity_factor}{Computed annuity factor.}
#'   \item{payment}{Level payment amount.}
#'   \item{present_value}{Present value of the annuity.}
#' }
#'
#' @seealso \code{\link{a_angle}}, \code{\link{s_angle}}, \code{\link{standardize_interest}}
#'
#' @family annuities
#'
#' @examples
#' # Simple scalar example
#' a_angle_tbl(n_years = 10, rate = 0.05, payment = 1000)
#'
#' # Medium vectorized example
#' a_angle_tbl(
#'   n_years = c(10, 15, 20),
#'   payments_per_year = c(1, 12, 1),
#'   rate = c(0.05, 0.06, 0.04),
#'   type = c("effective", "nominal_interest", "force"),
#'   m = c(1, 12, 1),
#'   deferral_years = c(0, 0, 2),
#'   timing = c("immediate", "immediate", "continuous"),
#'   payment = c(1000, 200, 5000),
#'   perpetuity = c(FALSE, FALSE, FALSE)
#' )
#'
#' @export
a_angle_tbl <- function(
    n_years = NULL,
    payments_per_year = 1L,
    rate,
    type = "effective",
    m = 1L,
    deferral_years = 0,
    timing = "immediate",
    payment = 1,
    perpetuity = FALSE
) {
  if (missing(rate)) {
    stop("`rate` must be provided.", call. = FALSE)
  }
  if (!is.numeric(payment)) {
    stop("`payment` must be a numeric vector.", call. = FALSE)
  }

  # --- Determine common size (fixed: use length, not value, for n_years) ---
  size <- max(
    if (is.null(n_years)) 0L else length(n_years),
    length(payments_per_year),
    length(rate),
    length(type),
    length(m),
    length(deferral_years),
    length(timing),
    length(payment),
    length(perpetuity),
    1L
  )

  valid_size <- function(x) is.null(x) || length(x) %in% c(1L, size)

  if (!valid_size(n_years) || !valid_size(payments_per_year) || !valid_size(rate) ||
      !valid_size(type) || !valid_size(m) || !valid_size(deferral_years) ||
      !valid_size(timing) || !valid_size(payment) || !valid_size(perpetuity)) {
    stop(
      "`n_years`, `payments_per_year`, `rate`, `type`, `m`, ",
      "`deferral_years`, `timing`, `payment`, and `perpetuity` must have ",
      "length 1 or a common length.",
      call. = FALSE
    )
  }

  if (is.null(n_years)) {
    n_years <- rep(NA_real_, size)
  } else {
    n_years <- rep_len(n_years, size)
  }

  payments_per_year <- rep_len(payments_per_year, size)
  rate              <- rep_len(rate, size)
  type              <- rep_len(type, size)
  m                 <- rep_len(m, size)
  deferral_years    <- rep_len(deferral_years, size)
  timing            <- rep_len(timing, size)
  payment           <- rep_len(payment, size)
  perpetuity        <- rep_len(perpetuity, size)

  bad_payment <- !is.na(payment) & !is.finite(payment)
  if (any(bad_payment)) {
    stop("`payment` must contain only finite numeric values or NA.", call. = FALSE)
  }

  annuity_factor <- a_angle(
    n_years = n_years,
    payments_per_year = payments_per_year,
    rate = rate,
    type = type,
    m = m,
    deferral_years = deferral_years,
    timing = timing,
    perpetuity = perpetuity
  )

  i_effective <- standardize_interest(type = type, rate = rate, m = m)
  delta <- log1p(i_effective)

  discrete <- tolower(timing) %in% c("immediate", "due")

  i_period <- rep(NA_real_, size)
  v_period <- rep(NA_real_, size)
  n_periods <- rep(NA_real_, size)
  deferral_periods <- rep(NA_real_, size)

  if (any(discrete, na.rm = TRUE)) {
    k <- payments_per_year[discrete]
    i_period[discrete] <- (1 + i_effective[discrete])^(1 / k) - 1
    v_period[discrete] <- 1 / (1 + i_period[discrete])
    deferral_periods[discrete] <- deferral_years[discrete] * k

    finite_disc <- discrete & !perpetuity
    n_periods[finite_disc] <- n_years[finite_disc] * payments_per_year[finite_disc]
  }

  tibble::tibble(
    n_years = n_years,
    payments_per_year = ifelse(tolower(timing) == "continuous", NA_integer_, as.integer(payments_per_year)),
    deferral_years = deferral_years,
    timing = tolower(timing),
    perpetuity = perpetuity,
    rate_input = rate,
    rate_type = type,
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
