#' Level annuity factor a-angle-n
#'
#' Computes the actuarial present value factor for a level annuity using compact
#' actuarial notation.
#'
#' Supported timing conventions:
#' \itemize{
#'   \item \code{"immediate"}: annuity-immediate with discrete payments.
#'   \item \code{"due"}: annuity-due with discrete payments.
#'   \item \code{"continuous"}: continuous annuity.
#' }
#'
#' For discrete annuities, \code{k} is the number of payments per year, so
#' payments are made every \eqn{1/k} year. The function returns the annuity
#' factor, assuming a unit payment at each payment time.
#'
#' Deferment is supported through \code{h}. For discrete annuities, the
#' deferment must align with the payment grid, that is, \eqn{hk} must be an
#' integer.
#'
#' If \code{perpetuity = TRUE}, the infinite-term annuity factor is returned.
#'
#' @param n Numeric vector of payment durations in years. Ignored when
#'   \code{perpetuity = TRUE}. If \code{perpetuity = FALSE}, each value must be
#'   positive and finite.
#' @param k Positive integer vector giving the number of discrete payments per
#'   year. Ignored for continuous annuities.
#' @param i Numeric vector of interest-rate values.
#' @param i_type Character vector indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the conversion frequency for nominal
#'   rates. Ignored for \code{"effective"} and \code{"force"}.
#' @param h Numeric vector of deferment times in years. Must be greater than or
#'   equal to 0.
#' @param timing Character vector. One of \code{"immediate"}, \code{"due"}, or
#'   \code{"continuous"}.
#' @param perpetuity Logical vector. If \code{TRUE}, computes the perpetuity
#'   factor.
#' @param payment Numeric vector of level payment amounts. Used only when
#'   \code{tidy = TRUE} to report the corresponding present value. The annuity
#'   factor itself is always computed for unit payments.
#' @param tidy Logical scalar. If \code{FALSE}, returns a numeric annuity factor.
#'   If \code{TRUE}, returns a tibble with intermediate calculations.
#'
#' @return
#' If \code{tidy = FALSE}, a numeric vector of annuity factors.
#'
#' If \code{tidy = TRUE}, a tibble with input values, equivalent rates, annuity
#' factors, payment amounts, and present values.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}:
#'
#' \itemize{
#'   \item \code{n}: annuity term;
#'   \item \code{k}: payment frequency;
#'   \item \code{i}: interest rate;
#'   \item \code{i_type}: interest-rate type;
#'   \item \code{m}: conversion frequency for nominal rates;
#'   \item \code{h}: deferment period.
#' }
#'
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
#' a_angle(n = 10, i = 0.05)
#'
#' # Nominal interest converted monthly, with monthly payments
#' a_angle(
#'   n = 10,
#'   i = 0.06,
#'   i_type = "nominal_interest",
#'   m = 12,
#'   k = 12
#' )
#'
#' # Continuous annuity
#' a_angle(
#'   n = 15,
#'   i = 0.04,
#'   i_type = "force",
#'   timing = "continuous"
#' )
#'
#' # Tibble output for teaching or auditing
#' a_angle(
#'   n = 10,
#'   i = 0.05,
#'   payment = 1000,
#'   tidy = TRUE
#' )
#'
#' # Vectorized example
#' a_angle(
#'   n = c(5, 10, 20),
#'   k = c(1, 12, 1),
#'   i = c(0.05, 0.06, 0.04),
#'   i_type = c("effective", "nominal_interest", "force"),
#'   m = c(1, 12, 1),
#'   h = c(0, 0, 2),
#'   timing = c("immediate", "immediate", "continuous"),
#'   perpetuity = c(FALSE, FALSE, FALSE)
#' )
#'
#' @export
a_angle <- function(
    n = NULL,
    k = 1L,
    i,
    i_type = "effective",
    m = 1L,
    h = 0,
    timing = "immediate",
    perpetuity = FALSE,
    payment = 1,
    tidy = FALSE
) {
  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  if (missing(i)) {
    stop("`i` must be provided.", call. = FALSE)
  }

  # --- Early type validation ---
  if (!is.numeric(i)) {
    stop("`i` must be a numeric vector.", call. = FALSE)
  }
  if (!is.null(n) && !is.numeric(n)) {
    stop("`n` must be numeric or NULL.", call. = FALSE)
  }
  if (!is.numeric(k)) {
    stop("`k` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(m)) {
    stop("`m` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(h)) {
    stop("`h` must be numeric.", call. = FALSE)
  }
  if (!is.character(i_type)) {
    stop("`i_type` must be a character vector.", call. = FALSE)
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
    if (is.null(n)) 0L else length(n),
    length(k),
    length(i),
    length(i_type),
    length(m),
    length(h),
    length(timing),
    length(perpetuity),
    length(payment),
    1L
  )

  valid_size <- function(x) is.null(x) || length(x) %in% c(1L, size)

  if (!valid_size(n) ||
      !valid_size(k) ||
      !valid_size(i) ||
      !valid_size(i_type) ||
      !valid_size(m) ||
      !valid_size(h) ||
      !valid_size(timing) ||
      !valid_size(perpetuity) ||
      !valid_size(payment)) {
    stop(
      "`n`, `k`, `i`, `i_type`, `m`, `h`, `timing`, `perpetuity`, ",
      "and `payment` must have length 1 or a common length.",
      call. = FALSE
    )
  }

  # --- Recycle to common size ---
  if (is.null(n)) {
    n <- rep(NA_real_, size)
  } else {
    n <- rep_len(n, size)
  }

  k          <- rep_len(k, size)
  i          <- rep_len(i, size)
  i_type     <- rep_len(i_type, size)
  m          <- rep_len(m, size)
  h          <- rep_len(h, size)
  timing     <- rep_len(timing, size)
  perpetuity <- rep_len(perpetuity, size)
  payment    <- rep_len(payment, size)

  # --- Value-level validation ---
  bad_h <- !is.na(h) & (!is.finite(h) | h < 0)

  if (any(bad_h)) {
    stop("`h` must contain only finite values >= 0 or NA.", call. = FALSE)
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
  i_effective <- standardize_interest(type = i_type, rate = i, m = m)
  delta <- log1p(i_effective)

  annuity_factor <- rep(NA_real_, size)
  eps <- 1e-12

  # --- Identify element-level NA ---
  is_na_elem <- is.na(i) |
    is.na(i_effective) |
    is.na(h) |
    is.na(timing)

  # --- Validate n for finite annuities ---
  needs_n <- !is_na_elem & !perpetuity
  bad_n <- needs_n & (is.na(n) | !is.finite(n) | n <= 0)

  if (any(bad_n)) {
    stop(
      "When `perpetuity = FALSE`, `n` must be positive and finite.",
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
      exp(-delta[cont_perp] * h[cont_perp]) / delta[cont_perp]
  }

  cont_fin <- is_cont & !perpetuity
  if (any(cont_fin)) {
    d_cf <- delta[cont_fin]
    n_cf <- n[cont_fin]
    h_cf <- h[cont_fin]

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
  h_periods <- rep(NA_real_, size)

  if (length(disc_idx) > 0L) {
    for (j in disc_idx) {
      kk <- k[j]

      if (is.na(kk) || !is.finite(kk) || kk < 1 || kk != floor(kk)) {
        stop(
          "`k` must be a positive integer for discrete annuities.",
          call. = FALSE
        )
      }

      kk <- as.integer(kk)

      h_periods_raw <- h[j] * kk
      h_periods[j] <- round(h_periods_raw)

      if (abs(h_periods_raw - h_periods[j]) > 1e-10) {
        stop(
          "For discrete annuities, `h * k` must be an integer.",
          call. = FALSE
        )
      }

      i_period[j] <- (1 + i_effective[j])^(1 / kk) - 1
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

        annuity_factor[j] <- base * v_period[j]^h_periods[j]
        next
      }

      n_periods_raw <- n[j] * kk
      n_periods[j] <- round(n_periods_raw)

      if (abs(n_periods_raw - n_periods[j]) > 1e-10) {
        stop(
          "For discrete annuities, `n * k` must be an integer.",
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

      annuity_factor[j] <- base * v_period[j]^h_periods[j]
    }
  }

  if (!tidy) {
    return(annuity_factor)
  }

  tibble::tibble(
    n = n,
    k = ifelse(
      timing == "continuous",
      NA_integer_,
      as.integer(k)
    ),
    h = h,
    timing = timing,
    perpetuity = perpetuity,
    i_input = i,
    i_type = i_type,
    m = m,
    i_effective = i_effective,
    delta = delta,
    i_period = i_period,
    v_period = v_period,
    n_periods = n_periods,
    h_periods = h_periods,
    annuity_factor = annuity_factor,
    payment = payment,
    present_value = payment * annuity_factor
  )
}
