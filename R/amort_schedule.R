#' Amortization schedule with optional prepayment adjustment
#'
#' Builds an amortization schedule under a fixed annual interest-rate
#' specification, allowing extra principal payments and optional adjustment of
#' either the remaining term or the remaining payment amount.
#'
#' The annual rate is converted internally to an effective rate per schedule
#' period using \code{k}.
#'
#' Adjustment policies after extra principal payments:
#' \itemize{
#'   \item \code{"none"}: keep the original payment and contractual term,
#'         unless the loan is fully repaid early.
#'   \item \code{"term"}: keep the regular payment and shorten the term.
#'         No special logic is needed: the loop exits naturally when the
#'         outstanding balance reaches zero.
#'   \item \code{"payment"}: keep the remaining contractual term and
#'         recalculate the regular payment after each period.
#' }
#'
#' @param principal Numeric scalar. Initial outstanding balance.
#' @param n Positive integer. Number of contractual periods.
#' @param i Numeric scalar. Annual interest-rate input.
#' @param i_type Character string indicating the annual interest-rate type:
#'   \code{"effective"}, \code{"nominal_interest"}, \code{"nominal_discount"},
#'   or \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal annual rates.
#' @param k Positive integer. Number of amortization periods per year.
#' @param timing Character string. One of \code{"immediate"} or \code{"due"}.
#' @param payment Optional numeric scalar. Initial regular payment per period.
#'   If \code{NULL}, it is computed from the loan data.
#' @param extra_principal Optional extra principal payments. Can be:
#'   \itemize{
#'     \item \code{NULL},
#'     \item a scalar,
#'     \item an unnamed numeric vector of length \code{n},
#'     \item a named numeric vector with names interpreted as period numbers.
#'   }
#' @param adjust Character string. One of \code{"none"}, \code{"term"}, or
#'   \code{"payment"}.
#' @param tol Numeric tolerance for zero-balance detection.
#'
#' @return A tibble with one row per realized period and columns:
#' \describe{
#'   \item{period}{Period index.}
#'   \item{ob_start}{Outstanding balance at the start of the period.}
#'   \item{interest}{Interest charged during the period.}
#'   \item{payment}{Regular payment in the period.}
#'   \item{extra_principal}{Extra principal paid in the period.}
#'   \item{principal}{Principal repaid through the regular payment.}
#'   \item{total_principal}{Total principal repaid in the period.}
#'   \item{cashflow}{Total payment made in the period.}
#'   \item{ob_end}{Outstanding balance at the end of the period.}
#'   \item{i_effective_annual}{Equivalent annual effective rate.}
#'   \item{i_effective_period}{Equivalent effective rate per schedule period.}
#'   \item{k}{Schedule frequency.}
#'   \item{timing}{Payment timing convention.}
#'   \item{adjust}{Adjustment rule used.}
#' }
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{i} is the interest rate, \code{i_type} is the
#' interest-rate type, \code{m} is the conversion frequency for nominal annual
#' rates, and \code{k} is the number of amortization periods per year.
#'
#' For \code{timing = "immediate"} (annuity-immediate), interest accrues on the
#' outstanding balance during the period, and the payment is made at the end.
#' For \code{timing = "due"} (annuity-due), the payment is made at the start of
#' the period and interest accrues on the balance after the payment.
#'
#' If the user supplies a custom \code{payment} that is smaller than the
#' periodic interest, principal repayment will be negative (negative
#' amortization). This is permitted but the user should be aware.
#'
#' @seealso \code{\link{a_angle}}, \code{\link{present_value}},
#'   \code{\link{pv_flow}}, \code{\link{standardize_interest}}
#'
#' @family amortization
#'
#' @examples
#' amort_schedule(
#'   principal = 100000,
#'   n = 12,
#'   i = 0.12,
#'   i_type = "nominal_interest",
#'   m = 12,
#'   k = 12
#' )
#'
#' amort_schedule(
#'   principal = 100000,
#'   n = 24,
#'   i = 0.12,
#'   i_type = "nominal_interest",
#'   m = 12,
#'   k = 12,
#'   extra_principal = c("6" = 5000, "12" = 3000),
#'   adjust = "term"
#' )
#'
#' amort_schedule(
#'   principal = 100000,
#'   n = 24,
#'   i = 0.12,
#'   i_type = "nominal_interest",
#'   m = 12,
#'   k = 12,
#'   extra_principal = c("6" = 5000, "12" = 3000),
#'   adjust = "payment"
#' )
#'
#' @export
amort_schedule <- function(
    principal,
    n,
    i,
    i_type = "effective",
    m = 1L,
    k = 1L,
    timing = c("immediate", "due"),
    payment = NULL,
    extra_principal = NULL,
    adjust = c("none", "term", "payment"),
    tol = 1e-8
) {
  timing <- match.arg(timing)
  adjust <- match.arg(adjust)

  # --- Scalar input validation ---
  if (!is.numeric(principal) || length(principal) != 1L || is.na(principal) ||
      !is.finite(principal) || principal <= 0) {
    stop("`principal` must be a single finite positive number.", call. = FALSE)
  }

  if (!is.numeric(n) || length(n) != 1L || is.na(n) || !is.finite(n) ||
      n <= 0 || n != floor(n)) {
    stop("`n` must be a single positive integer.", call. = FALSE)
  }
  n <- as.integer(n)

  if (missing(i) || !is.numeric(i) || length(i) != 1L || is.na(i) || !is.finite(i)) {
    stop("`i` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.character(i_type) || length(i_type) != 1L || is.na(i_type)) {
    stop("`i_type` must be a single character string.", call. = FALSE)
  }

  if (!is.numeric(m) || length(m) != 1L || is.na(m) || !is.finite(m) ||
      m <= 0 || m != floor(m)) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }
  m <- as.integer(m)

  if (!is.numeric(k) || length(k) != 1L ||
      is.na(k) || !is.finite(k) ||
      k <= 0 || k != floor(k)) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }
  k <- as.integer(k)

  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("`tol` must be a single positive numeric value.", call. = FALSE)
  }

  # --- Rate conversion ---
  i_effective_annual <- standardize_interest(
    type = i_type,
    rate = i,
    m = m
  )

  i_effective_period <- (1 + i_effective_annual)^(1 / k) - 1

  if (!is.finite(i_effective_period) || i_effective_period <= -1) {
    stop(
      "The derived effective rate per period must be finite and greater than -1.",
      call. = FALSE
    )
  }

  # --- Helper: compute level payment for a given balance / remaining term ---
  level_payment <- function(balance, n_remaining, i_p, timing, tol) {
    if (n_remaining <= 0) return(0)

    if (abs(i_p) < tol) {
      return(balance / n_remaining)
    }

    a_n <- (1 - (1 + i_p)^(-n_remaining)) / i_p

    if (timing == "immediate") {
      balance / a_n
    } else {
      balance / ((1 + i_p) * a_n)
    }
  }

  # --- Determine initial regular payment ---
  if (is.null(payment)) {
    current_payment <- level_payment(
      balance = principal,
      n_remaining = n,
      i_p = i_effective_period,
      timing = timing,
      tol = tol
    )
  } else {
    if (!is.numeric(payment) || length(payment) != 1L || is.na(payment) ||
        !is.finite(payment) || payment <= 0) {
      stop("`payment` must be a single finite positive number.", call. = FALSE)
    }
    current_payment <- payment
  }

  # --- Build extra-principal vector ---
  extra_vec <- numeric(n)

  if (!is.null(extra_principal)) {
    if (!is.numeric(extra_principal) || any(is.na(extra_principal)) ||
        any(!is.finite(extra_principal))) {
      stop("`extra_principal` must be numeric and finite.", call. = FALSE)
    }
    if (any(extra_principal < 0)) {
      stop("`extra_principal` must be nonnegative.", call. = FALSE)
    }

    if (length(extra_principal) == 1L && is.null(names(extra_principal))) {
      extra_vec[] <- extra_principal

    } else if (is.null(names(extra_principal))) {
      if (length(extra_principal) != n) {
        stop("If `extra_principal` is unnamed, its length must equal `n`.", call. = FALSE)
      }
      extra_vec <- as.numeric(extra_principal)

    } else {
      idx <- suppressWarnings(as.integer(names(extra_principal)))
      if (any(is.na(idx)) || any(idx < 1L) || any(idx > n)) {
        stop(
          "Names of `extra_principal` must be integers between 1 and `n`.",
          call. = FALSE
        )
      }

      tmp <- tapply(as.numeric(extra_principal), idx, sum)
      extra_vec[as.integer(names(tmp))] <- as.numeric(tmp)
    }
  }

  # --- Pre-allocate output vectors (avoid tibble creation inside loop) ---
  out_period     <- integer(n)
  out_ob_start   <- numeric(n)
  out_interest   <- numeric(n)
  out_payment    <- numeric(n)
  out_extra      <- numeric(n)
  out_principal  <- numeric(n)
  out_total_prin <- numeric(n)
  out_cashflow   <- numeric(n)
  out_ob_end     <- numeric(n)

  ob <- principal
  n_realized <- 0L

  for (period_idx in seq_len(n)) {
    if (ob <= tol) break

    n_realized <- n_realized + 1L
    extra_req <- extra_vec[period_idx]
    payment_period <- current_payment

    if (timing == "immediate") {
      # Interest accrues on opening balance; payment at end of period
      interest_period <- ob * i_effective_period
      principal_reg_period <- payment_period - interest_period
      total_principal_period <- principal_reg_period + extra_req
      cashflow_period <- payment_period + extra_req
      ob_new <- ob + interest_period - cashflow_period

      # Cap at zero if overpaid
      if (ob_new <= tol) {
        total_needed <- ob + interest_period
        extra_used <- min(extra_req, total_needed)
        payment_period <- total_needed - extra_used
        principal_reg_period <- payment_period - interest_period
        total_principal_period <- principal_reg_period + extra_used
        cashflow_period <- payment_period + extra_used
        ob_new <- 0
        extra_req <- extra_used
      }

    } else {
      # Annuity-due: payment at start of period; interest on balance after payment
      cashflow_req <- payment_period + extra_req
      ob_after_pay <- ob - cashflow_req

      if (ob_after_pay <= tol) {
        # Full liquidation at start of period - no interest accrues
        total_needed_now <- ob
        extra_used <- min(extra_req, total_needed_now)
        payment_period <- total_needed_now - extra_used
        interest_period <- 0
        principal_reg_period <- payment_period
        total_principal_period <- principal_reg_period + extra_used
        cashflow_period <- payment_period + extra_used
        ob_new <- 0
        extra_req <- extra_used
      } else {
        interest_period <- ob_after_pay * i_effective_period
        principal_reg_period <- payment_period
        total_principal_period <- principal_reg_period + extra_req
        cashflow_period <- payment_period + extra_req
        ob_new <- ob_after_pay + interest_period
      }
    }

    out_period[period_idx]     <- period_idx
    out_ob_start[period_idx]   <- ob
    out_interest[period_idx]   <- interest_period
    out_payment[period_idx]    <- payment_period
    out_extra[period_idx]      <- extra_req
    out_principal[period_idx]  <- principal_reg_period
    out_total_prin[period_idx] <- total_principal_period
    out_cashflow[period_idx]   <- cashflow_period
    out_ob_end[period_idx]     <- ob_new

    ob <- ob_new

    if (ob <= tol) break

    # For "payment" adjustment: recalculate level payment with new balance.
    # For "term" adjustment: no recalculation needed - the loop exits naturally
    # once the balance reaches zero, effectively shortening the term.
    if (adjust == "payment") {
      remaining_periods <- n - period_idx
      current_payment <- level_payment(
        balance = ob,
        n_remaining = remaining_periods,
        i_p = i_effective_period,
        timing = timing,
        tol = tol
      )
    }
  }

  # --- Build output tibble from pre-allocated vectors ---
  if (n_realized == 0L) {
    out <- tibble::tibble(
      period = integer(0),
      ob_start = numeric(0),
      interest = numeric(0),
      payment = numeric(0),
      extra_principal = numeric(0),
      principal = numeric(0),
      total_principal = numeric(0),
      cashflow = numeric(0),
      ob_end = numeric(0),
      i_effective_annual = numeric(0),
      i_effective_period = numeric(0),
      k = integer(0),
      timing = character(0),
      adjust = character(0)
    )
  } else {
    idx <- seq_len(n_realized)
    out <- tibble::tibble(
      period          = out_period[idx],
      ob_start        = out_ob_start[idx],
      interest        = out_interest[idx],
      payment         = out_payment[idx],
      extra_principal = out_extra[idx],
      principal       = out_principal[idx],
      total_principal = out_total_prin[idx],
      cashflow        = out_cashflow[idx],
      ob_end          = out_ob_end[idx],
      i_effective_annual = i_effective_annual,
      i_effective_period = i_effective_period,
      k = k,
      timing = timing,
      adjust = adjust
    )
  }

  out
}
