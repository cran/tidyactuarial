#' Arithmetic annuity factor
#'
#' Computes the actuarial present value or accumulated value factor for an
#' arithmetic annuity, using compact actuarial notation.
#'
#' This function covers increasing, decreasing, and custom arithmetic payment
#' patterns. It replaces the more specific increasing and decreasing annuity
#' functions.
#'
#' Supported payment patterns:
#' \itemize{
#'   \item \code{"increasing"}: payments \eqn{1, 2, \ldots, N}.
#'   \item \code{"decreasing"}: payments \eqn{N, N-1, \ldots, 1}.
#'   \item \code{"custom"}: payments following
#'     \eqn{P_r = P_1 + (r - 1)g}.
#' }
#'
#' Supported timing conventions:
#' \itemize{
#'   \item \code{"immediate"}: payments at the end of each period.
#'   \item \code{"due"}: payments at the beginning of each period.
#' }
#'
#' @param n Numeric vector of payment durations in years. Each value must be
#'   positive and finite.
#' @param k Positive integer vector giving the number of payments per year.
#' @param i Numeric vector of interest-rate values.
#' @param i_type Character vector indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the conversion frequency for nominal
#'   rates. Ignored for \code{"effective"} and \code{"force"}.
#' @param h Numeric vector of deferment times in years. Must be greater than or
#'   equal to 0.
#' @param timing Character vector. One of \code{"immediate"} or \code{"due"}.
#' @param pattern Character vector. One of \code{"increasing"},
#'   \code{"decreasing"}, or \code{"custom"}.
#' @param P1 Numeric vector. First payment for \code{pattern = "custom"}.
#'   Ignored for \code{"increasing"} and \code{"decreasing"}.
#' @param g Numeric vector. Arithmetic increment for
#'   \code{pattern = "custom"}. Ignored for \code{"increasing"} and
#'   \code{"decreasing"}.
#' @param valuation Character string. Use \code{"present"} for present value
#'   or \code{"accumulated"} for accumulated value at the end of the term.
#' @param tidy Logical scalar. If \code{FALSE}, returns a numeric factor. If
#'   \code{TRUE}, returns a tibble with intermediate quantities.
#'
#' @return
#' If \code{tidy = FALSE}, a numeric vector of arithmetic annuity factors.
#'
#' If \code{tidy = TRUE}, a tibble with input values, equivalent rates, period
#' quantities, and both present and accumulated value factors.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{n} is the term, \code{k} is the payment
#' frequency, \code{i} is the interest-rate input, \code{i_type} is the
#' interest-rate type, \code{m} is the conversion frequency for nominal rates,
#' and \code{h} is the deferment period.
#'
#' The function first converts the supplied rate to the equivalent annual
#' effective interest rate using \code{\link{standardize_interest}}.
#'
#' For each scenario, the total number of payments is
#' \deqn{N = n k.}
#'
#' The present value is computed by summing the discounted payment stream. The
#' accumulated value is computed at the end of the annuity term. Under this
#' convention, \code{h} affects the present value factor but not the
#' accumulated value factor.
#'
#' @seealso \code{\link{a_angle}}, \code{\link{s_angle}},
#'   \code{\link{standardize_interest}}
#'
#' @family annuities
#'
#' @examples
#' # Increasing arithmetic annuity
#' annuity_arith(n = 10, i = 0.05, pattern = "increasing")
#'
#' # Decreasing arithmetic annuity
#' annuity_arith(n = 10, i = 0.05, pattern = "decreasing")
#'
#' # Custom arithmetic annuity
#' annuity_arith(
#'   n = 10,
#'   i = 0.05,
#'   pattern = "custom",
#'   P1 = 100,
#'   g = 25
#' )
#'
#' # Accumulated value factor
#' annuity_arith(
#'   n = 10,
#'   i = 0.05,
#'   pattern = "increasing",
#'   valuation = "accumulated"
#' )
#'
#' # Tibble output
#' annuity_arith(
#'   n = 10,
#'   i = 0.05,
#'   pattern = "decreasing",
#'   tidy = TRUE
#' )
#'
#' @export
annuity_arith <- function(
    n,
    k = 1L,
    i,
    i_type = "effective",
    m = 1L,
    h = 0,
    timing = "immediate",
    pattern = "increasing",
    P1 = 1,
    g = 1,
    valuation = c("present", "accumulated"),
    tidy = FALSE
) {
  valuation <- match.arg(valuation)

  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  if (missing(n)) {
    stop("`n` must be provided.", call. = FALSE)
  }
  if (missing(i)) {
    stop("`i` must be provided.", call. = FALSE)
  }

  # --- Early type validation ---
  if (!is.numeric(n)) {
    stop("`n` must be a numeric vector.", call. = FALSE)
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
  if (!is.numeric(h)) {
    stop("`h` must be numeric.", call. = FALSE)
  }
  if (!is.character(timing)) {
    stop("`timing` must be a character vector.", call. = FALSE)
  }
  if (!is.character(pattern)) {
    stop("`pattern` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(P1)) {
    stop("`P1` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(g)) {
    stop("`g` must be numeric.", call. = FALSE)
  }

  # --- Determine common size ---
  size <- max(
    length(n),
    length(k),
    length(i),
    length(i_type),
    length(m),
    length(h),
    length(timing),
    length(pattern),
    length(P1),
    length(g),
    1L
  )

  valid_size <- function(x) length(x) %in% c(1L, size)

  if (!valid_size(n) ||
      !valid_size(k) ||
      !valid_size(i) ||
      !valid_size(i_type) ||
      !valid_size(m) ||
      !valid_size(h) ||
      !valid_size(timing) ||
      !valid_size(pattern) ||
      !valid_size(P1) ||
      !valid_size(g)) {
    stop(
      "`n`, `k`, `i`, `i_type`, `m`, `h`, `timing`, `pattern`, ",
      "`P1`, and `g` must have length 1 or a common length.",
      call. = FALSE
    )
  }

  # --- Recycle ---
  n       <- rep_len(n, size)
  k       <- rep_len(k, size)
  i       <- rep_len(i, size)
  i_type  <- rep_len(i_type, size)
  m       <- rep_len(m, size)
  h       <- rep_len(h, size)
  timing  <- rep_len(timing, size)
  pattern <- rep_len(pattern, size)
  P1      <- rep_len(P1, size)
  g       <- rep_len(g, size)

  # --- Value-level validation ---
  bad_n <- !is.na(n) & (!is.finite(n) | n <= 0)
  if (any(bad_n)) {
    stop("`n` must contain only finite values greater than 0 or NA.", call. = FALSE)
  }

  bad_h <- !is.na(h) & (!is.finite(h) | h < 0)
  if (any(bad_h)) {
    stop("`h` must contain only finite values >= 0 or NA.", call. = FALSE)
  }

  bad_P1 <- !is.na(P1) & !is.finite(P1)
  if (any(bad_P1)) {
    stop("`P1` must contain only finite numeric values or NA.", call. = FALSE)
  }

  bad_g <- !is.na(g) & !is.finite(g)
  if (any(bad_g)) {
    stop("`g` must contain only finite numeric values or NA.", call. = FALSE)
  }

  bad_m <- !is.na(m) & (!is.finite(m) | m < 1 | abs(m - round(m)) > 1e-10)
  if (any(bad_m)) {
    stop("`m` must contain positive integer values or NA.", call. = FALSE)
  }

  timing <- tolower(timing)
  valid_timing <- c("immediate", "due")

  bad_timing <- !is.na(timing) & !(timing %in% valid_timing)
  if (any(bad_timing)) {
    stop("`timing` must be 'immediate' or 'due'.", call. = FALSE)
  }

  pattern <- tolower(pattern)
  valid_pattern <- c("increasing", "decreasing", "custom")

  bad_pattern <- !is.na(pattern) & !(pattern %in% valid_pattern)
  if (any(bad_pattern)) {
    stop("`pattern` must be 'increasing', 'decreasing', or 'custom'.", call. = FALSE)
  }

  # --- Convert to effective annual rate ---
  i_effective <- standardize_interest(i_type = i_type, i = i, m = m)

  if (any(!is.na(i_effective) & (!is.finite(i_effective) | i_effective <= -1))) {
    stop(
      "The standardized annual effective interest rates must be greater than -1.",
      call. = FALSE
    )
  }

  present_value_factor <- rep(NA_real_, size)
  accumulated_value_factor <- rep(NA_real_, size)
  selected_factor <- rep(NA_real_, size)

  i_period <- rep(NA_real_, size)
  v_period <- rep(NA_real_, size)
  n_periods <- rep(NA_real_, size)
  h_periods <- rep(NA_real_, size)
  P1_used <- rep(NA_real_, size)
  g_used <- rep(NA_real_, size)

  is_na_elem <- is.na(n) |
    is.na(k) |
    is.na(i) |
    is.na(i_effective) |
    is.na(h) |
    is.na(timing) |
    is.na(pattern)

  computable_idx <- which(!is_na_elem)

  for (idx in computable_idx) {
    kk <- k[idx]

    if (is.na(kk) || !is.finite(kk) || kk < 1 || kk != floor(kk)) {
      stop("`k` must be a positive integer.", call. = FALSE)
    }

    kk <- as.integer(kk)

    n_raw <- n[idx] * kk
    N <- round(n_raw)

    if (abs(n_raw - N) > 1e-10) {
      stop(
        "For discrete annuities, `n * k` must be an integer.",
        call. = FALSE
      )
    }

    h_raw <- h[idx] * kk
    H <- round(h_raw)

    if (abs(h_raw - H) > 1e-10) {
      stop(
        "For discrete annuities, `h * k` must be an integer.",
        call. = FALSE
      )
    }

    i_period[idx] <- (1 + i_effective[idx])^(1 / kk) - 1
    v_period[idx] <- 1 / (1 + i_period[idx])
    n_periods[idx] <- N
    h_periods[idx] <- H

    if (pattern[idx] == "increasing") {
      p1 <- 1
      inc <- 1
    } else if (pattern[idx] == "decreasing") {
      p1 <- N
      inc <- -1
    } else {
      if (is.na(P1[idx]) || is.na(g[idx])) {
        stop(
          "`P1` and `g` must not be NA when `pattern = 'custom'`.",
          call. = FALSE
        )
      }

      p1 <- P1[idx]
      inc <- g[idx]
    }

    P1_used[idx] <- p1
    g_used[idx] <- inc

    payment_index <- seq_len(N)
    payments <- p1 + (payment_index - 1) * inc

    local_times <- if (timing[idx] == "immediate") {
      payment_index
    } else {
      payment_index - 1
    }

    present_value_factor[idx] <- sum(payments * v_period[idx]^(H + local_times))

    accumulated_value_factor[idx] <- sum(
      payments * (1 + i_period[idx])^(N - local_times)
    )

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
    pattern = pattern,
    valuation = valuation,
    i_input = i,
    i_type = i_type,
    m = m,
    i_effective = i_effective,
    i_period = i_period,
    v_period = v_period,
    n_periods = n_periods,
    h_periods = h_periods,
    P1 = P1_used,
    g = g_used,
    present_value_factor = present_value_factor,
    accumulated_value_factor = accumulated_value_factor,
    annuity_factor = selected_factor
  )
}

