#' Geometric annuity factor
#'
#' Computes the actuarial present value or accumulated value factor for a
#' geometric annuity, using compact actuarial notation.
#'
#' This function covers finite geometric annuities and geometric perpetuities.
#'
#' Supported timing conventions:
#' \itemize{
#'   \item \code{"immediate"}: payments at the end of each period.
#'   \item \code{"due"}: payments at the beginning of each period.
#' }
#'
#' @param n Numeric vector of payment durations in years. Ignored only when
#'   \code{perpetuity = TRUE} and \code{valuation = "present"}. Otherwise,
#'   each value must be positive and finite.
#' @param k Positive integer vector giving the number of payments per year.
#' @param i Numeric vector of interest-rate values.
#' @param i_type Character vector indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the conversion frequency for nominal
#'   interest-rate inputs. Observed values are validated for every rate family,
#'   although only nominal conversions use this frequency.
#' @param g Numeric vector of annual growth-rate values for the payments.
#' @param g_type Character vector indicating the growth-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param g_m Positive integer vector giving the conversion frequency for
#'   nominal growth-rate inputs. Observed values are validated for every rate
#'   family, although only nominal conversions use this frequency.
#' @param h Numeric vector of deferment times in years. Must be greater than or
#'   equal to 0. For present values, the deferment discounts the payment block.
#'   For accumulated values, it is recorded but does not change the factor under
#'   the adopted terminal-horizon convention.
#' @param timing Character vector. One of \code{"immediate"} or \code{"due"}.
#' @param P1 Numeric vector giving the first payment of the geometric sequence.
#' @param perpetuity Logical vector. If \code{TRUE}, computes the geometric
#'   perpetuity present value. Perpetuities are not supported for
#'   \code{valuation = "accumulated"}.
#' @param valuation Character string. Use \code{"present"} for present value
#'   or \code{"accumulated"} for accumulated value at the end of the term.
#' @param tidy Logical scalar. If \code{FALSE}, returns a numeric factor. If
#'   \code{TRUE}, returns a tibble with intermediate quantities.
#'
#' @return
#' If \code{tidy = FALSE}, a numeric vector.
#'
#' If \code{tidy = TRUE}, a tibble with input values, equivalent rates, period
#' quantities, and both present and accumulated value factors.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{n} is the term, \code{k} is the payment
#' frequency, \code{i} is the interest-rate input, \code{i_type} is the
#' interest-rate type, \code{m} is the conversion frequency for nominal
#' interest rates, \code{g} is the growth-rate input, \code{g_type} is the
#' growth-rate type, \code{g_m} is the conversion frequency for nominal growth
#' rates, \code{h} is the deferment period, and \code{P1} is the first payment.
#'
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
#' ga_N = \frac{1-\left(\frac{1+g_p}{1+i_p}\right)^N}{i_p-g_p}.
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
#' Valuation retains the original annual log rates \eqn{\delta_i} and
#' \eqn{\delta_g}. Let \eqn{r=\delta_i/k}, \eqn{s=\delta_g/k},
#' \eqn{H=hk}, and \eqn{b=0} for immediate or \eqn{b=1} for due payments.
#' With \eqn{G_N(z)=\sum_{j=0}^{N-1}\exp(zj)}, the finite present value is
#' \deqn{P_1\exp[-r(H+1-b)]G_N(s-r).}
#' The accumulated value at \eqn{h+n} years is
#' \deqn{P_1\exp[br+(N-1)\max(r,s)]G_N(-|s-r|).}
#' Stable logarithmic sums avoid subtracting nearly equal effective rates or
#' accumulation levels. The spread is computed before period conversion as
#' \eqn{(\delta_g-\delta_i)/k}; no nonzero spread is replaced by zero.
#' For matching effective rates or nominal rates with equal conversion
#' frequencies, log-ratio identities preserve close input differences before
#' taking logarithms. Widely separated ratios that lose range fall back to
#' the original log rates. Different rate families retain the accuracy limits
#' of their separate conversions.
#'
#' A perpetuity requires \eqn{\delta_i>\delta_g}, without an arbitrary
#' numerical gap. Its value is
#' \eqn{P_1\exp[-r(H+1-b)]/[1-\exp(s-r)]}.
#' Missing finite terms propagate missing values. Signed and zero first
#' payments retain their meaning. Standalone audit rates may saturate, while
#' monetary values combine the first payment with log factors before declaring
#' range loss. No arbitrary precision or universal extreme-range guarantee is
#' implied. Unlike \code{g} in \code{annuity_arith}, this \code{g} is an annual
#' growth-rate specification and not a monetary increment per payment.
#'
#' @seealso \code{\link{a_angle}}, \code{\link{s_angle}},
#'   \code{\link{annuity_arith}}, \code{\link{standardize_interest}}
#'
#' @family annuities
#'
#' @examples
#' # Present value of a geometric annuity
#' annuity_geom(
#'   n = 10,
#'   i = 0.05,
#'   g = 0.02,
#'   valuation = "present"
#' )
#'
#' # Accumulated value of a geometric annuity
#' annuity_geom(
#'   n = 10,
#'   i = 0.05,
#'   g = 0.02,
#'   valuation = "accumulated"
#' )
#'
#' # Nominal interest and nominal growth
#' annuity_geom(
#'   n = 10,
#'   k = 12,
#'   i = 0.06,
#'   i_type = "nominal_interest",
#'   m = 12,
#'   g = 0.024,
#'   g_type = "nominal_interest",
#'   g_m = 12
#' )
#'
#' # Tibble output
#' annuity_geom(
#'   n = 10,
#'   i = 0.05,
#'   g = 0.02,
#'   tidy = TRUE
#' )
#'
#' @export
annuity_geom <- function(
    n = NULL,
    k = 1L,
    i,
    i_type = "effective",
    m = 1L,
    g = 0,
    g_type = "effective",
    g_m = 1L,
    h = 0,
    timing = "immediate",
    P1 = 1,
    perpetuity = FALSE,
    valuation = c("present", "accumulated"),
    tidy = FALSE
) {
  valuation <- match.arg(valuation)

  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  if (missing(i)) {
    stop("`i` must be provided.", call. = FALSE)
  }

  # --- Early type validation ---
  if (!is.null(n) && !is.numeric(n)) {
    stop("`n` must be numeric or NULL.", call. = FALSE)
  }
  if (!is.numeric(k)) {
    stop("`k` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(i)) {
    stop("`i` must be a numeric vector.", call. = FALSE)
  }
  if (!is.character(i_type)) {
    stop("`i_type` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(m)) {
    stop("`m` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(g)) {
    stop("`g` must be a numeric vector.", call. = FALSE)
  }
  if (!is.character(g_type)) {
    stop("`g_type` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(g_m)) {
    stop("`g_m` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(h)) {
    stop("`h` must be numeric.", call. = FALSE)
  }
  if (!is.character(timing)) {
    stop("`timing` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(P1)) {
    stop("`P1` must be numeric.", call. = FALSE)
  }
  if (!is.logical(perpetuity)) {
    stop("`perpetuity` must be a logical vector.", call. = FALSE)
  }

  # --- Determine common size ---
  size <- max(
    if (is.null(n)) 0L else length(n),
    length(k),
    length(i),
    length(i_type),
    length(m),
    length(g),
    length(g_type),
    length(g_m),
    length(h),
    length(timing),
    length(P1),
    length(perpetuity),
    1L
  )

  valid_size <- function(x) is.null(x) || length(x) %in% c(1L, size)

  if (!valid_size(n) ||
      !valid_size(k) ||
      !valid_size(i) ||
      !valid_size(i_type) ||
      !valid_size(m) ||
      !valid_size(g) ||
      !valid_size(g_type) ||
      !valid_size(g_m) ||
      !valid_size(h) ||
      !valid_size(timing) ||
      !valid_size(P1) ||
      !valid_size(perpetuity)) {
    stop(
      "`n`, `k`, `i`, `i_type`, `m`, `g`, `g_type`, `g_m`, `h`, ",
      "`timing`, `P1`, and `perpetuity` must have length 1 or a common length.",
      call. = FALSE
    )
  }

  # --- Recycle ---
  if (is.null(n)) {
    n <- rep(NA_real_, size)
  } else {
    n <- rep_len(n, size)
  }

  k <- rep_len(k, size)
  i <- rep_len(i, size)
  i_type <- rep_len(i_type, size)
  m <- rep_len(m, size)
  g <- rep_len(g, size)
  g_type <- rep_len(g_type, size)
  g_m <- rep_len(g_m, size)
  h <- rep_len(h, size)
  timing <- rep_len(timing, size)
  P1 <- rep_len(P1, size)
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

  bad_n <- !is.na(n) & (!is.finite(n) | n <= 0)
  if (any(bad_n & !perpetuity)) {
    stop(
      "When `perpetuity = FALSE`, `n` must contain only finite values greater than 0.",
      call. = FALSE
    )
  }

  bad_h <- !is.na(h) & (!is.finite(h) | h < 0)
  if (any(bad_h)) {
    stop("`h` must contain only finite values >= 0 or NA.", call. = FALSE)
  }

  bad_P1 <- !is.na(P1) & !is.finite(P1)
  if (any(bad_P1)) {
    stop("`P1` must contain only finite numeric values or NA.", call. = FALSE)
  }

  bad_m <- !is.na(m) & (!is.finite(m) | m < 1 | abs(m - round(m)) > 1e-10)
  if (any(bad_m)) {
    stop("`m` must contain positive integer values or NA.", call. = FALSE)
  }

  bad_g_m <- !is.na(g_m) & (!is.finite(g_m) | g_m < 1 | abs(g_m - round(g_m)) > 1e-10)
  if (any(bad_g_m)) {
    stop("`g_m` must contain positive integer values or NA.", call. = FALSE)
  }

  timing <- tolower(timing)
  valid_timing <- c("immediate", "due")

  bad_timing <- !is.na(timing) & !(timing %in% valid_timing)
  if (any(bad_timing)) {
    stop("`timing` must be 'immediate' or 'due'.", call. = FALSE)
  }

  # --- Convert rates ---
  interest_basis <- .interest_basis(i_type = i_type, i = i, m = m)
  growth_basis <- .interest_basis(i_type = g_type, i = g, m = g_m)
  i_effective <- interest_basis$effective
  g_effective <- growth_basis$effective

  present_value_factor <- rep(NA_real_, size)
  accumulated_value_factor <- rep(NA_real_, size)
  selected_factor <- rep(NA_real_, size)

  i_period <- rep(NA_real_, size)
  g_period <- rep(NA_real_, size)
  v_period <- rep(NA_real_, size)
  n_periods <- rep(NA_real_, size)
  h_periods <- rep(NA_real_, size)

  is_na_elem <- (!perpetuity & is.na(n)) |
    is.na(i) |
    is.na(g) |
    is.na(i_effective) |
    is.na(g_effective) |
    is.na(h) |
    is.na(timing) |
    is.na(P1)

  computable_idx <- which(!is_na_elem)

  for (idx in computable_idx) {
    kk <- k[idx]

    if (is.na(kk) || !is.finite(kk) || kk < 1 || kk != floor(kk)) {
      stop("`k` must be a positive integer.", call. = FALSE)
    }

    kk <- as.integer(kk)

    h_raw <- h[idx] * kk
    H <- round(h_raw)

    if (abs(h_raw - H) > 1e-10) {
      stop(
        "For discrete annuities, `h * k` must be an integer.",
        call. = FALSE
      )
    }

    r <- interest_basis$growth[idx] / kk
    s <- growth_basis$growth[idx] / kk
    # Preserve the annual spread before period conversion. For matching rate
    # families, log-ratio identities avoid cancellation of rounded logarithms.
    annual_spread <- growth_basis$growth[idx] - interest_basis$growth[idx]
    candidate_spread <- annual_spread
    if (tolower(i_type[idx]) == tolower(g_type[idx])) {
      if (tolower(i_type[idx]) == "effective") {
        candidate_spread <- log1p((g[idx] - i[idx]) / (1 + i[idx]))
      } else if (tolower(i_type[idx]) == "nominal_interest" && m[idx] == g_m[idx] &&
                 is.finite(m[idx] + i[idx])) {
        candidate_spread <- m[idx] * log1p((g[idx] - i[idx]) / (m[idx] + i[idx]))
      } else if (tolower(i_type[idx]) == "nominal_discount" && m[idx] == g_m[idx] &&
                 is.finite(m[idx] - i[idx])) {
        candidate_spread <- -m[idx] * log1p((i[idx] - g[idx]) / (m[idx] - i[idx]))
      }
    }
    if (is.finite(candidate_spread)) annual_spread <- candidate_spread
    spread <- annual_spread / kk
    i_period[idx] <- expm1(r)
    g_period[idx] <- expm1(s)
    v_period[idx] <- exp(-r)
    h_periods[idx] <- H
    payment_shift <- H + as.integer(timing[idx] == "immediate")

    if (perpetuity[idx]) {
      if (annual_spread >= 0) {
        stop(
          "For a geometric perpetuity, the per-period interest rate must be strictly greater than the per-period growth rate.",
          call. = FALSE
        )
      }

      log_pv <- -r * payment_shift - log(-expm1(spread))
      present_value_factor[idx] <- .cash_flow_value(P1[idx], log_pv)
      selected_factor[idx] <- present_value_factor[idx]
      next
    }

    n_raw <- n[idx] * kk
    N <- round(n_raw)

    if (abs(n_raw - N) > 1e-10) {
      stop(
        "For discrete annuities, `n * k` must be an integer.",
        call. = FALSE
      )
    }

    n_periods[idx] <- N

    # One payment has no growth exponent. Anchor accumulated terms at the
    # larger rate to avoid subtracting two overflowing accumulation levels.
    log_pv_sum <- if (N == 1) 0 else .annuity_log_sum(spread, N)
    log_fv_sum <- if (N == 1) 0 else .annuity_log_sum(-abs(spread), N)
    log_pv <- log_pv_sum - r * payment_shift
    log_fv <- log_fv_sum + max(r, s) * (N - 1) +
      if (timing[idx] == "due") r else 0

    present_value_factor[idx] <- .cash_flow_value(P1[idx], log_pv)
    accumulated_value_factor[idx] <- .cash_flow_value(P1[idx], log_fv)

    selected_factor[idx] <- if (valuation == "present") {
      present_value_factor[idx]
    } else {
      accumulated_value_factor[idx]
    }
  }

  if (!tidy) {
    return(selected_factor)
  }

  tibble::tibble(
    n = n,
    k = as.integer(k),
    h = h,
    timing = timing,
    perpetuity = perpetuity,
    valuation = valuation,
    i_input = i,
    i_type = i_type,
    m = m,
    g_input = g,
    g_type = g_type,
    g_m = g_m,
    i_effective = i_effective,
    g_effective = g_effective,
    i_period = i_period,
    g_period = g_period,
    v_period = v_period,
    n_periods = n_periods,
    h_periods = h_periods,
    P1 = P1,
    present_value_factor = present_value_factor,
    accumulated_value_factor = accumulated_value_factor,
    annuity_factor = selected_factor
  )
}
