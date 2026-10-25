#' Sinking fund amortization schedule for a loan
#'
#' Builds a sinking fund schedule under a fixed loan rate and a fixed
#' accumulation rate for the sinking fund, using compact actuarial notation.
#'
#' The borrower pays:
#' \itemize{
#'   \item interest on the loan each period, and
#'   \item a level deposit into the sinking fund.
#' }
#'
#' At maturity, the sinking fund is used to redeem the principal.
#'
#' @param principal Numeric scalar. Initial loan amount.
#' @param n Positive integer. Number of schedule periods.
#' @param i Numeric scalar. Annual interest-rate input for the loan.
#' @param j Numeric scalar. Annual accumulation-rate input for the sinking fund.
#' @param k Positive integer. Number of schedule periods per year.
#' @param i_type Character string indicating the loan interest-rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param j_type Character string indicating the sinking-fund accumulation-rate
#'   type. Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer. Conversion frequency for nominal loan rates.
#'   Validated for every rate family, although it affects conversion only for
#'   nominal rates.
#' @param j_m Positive integer. Conversion frequency for nominal sinking-fund
#'   rates. Validated for every rate family, although it affects conversion
#'   only for nominal rates.
#' @param deposit Optional numeric scalar. Level sinking-fund deposit per
#'   period. If \code{NULL}, it is computed so that the fund accumulates to
#'   \code{principal} at the end of period \code{n} (\code{n / k} years).
#' @param tol Finite positive numeric tolerance in monetary units used for
#'   final-balance checks.
#'
#' @return A tibble with one row per period and columns:
#' \describe{
#'   \item{period}{Period index.}
#'   \item{loan_balance_start}{Outstanding loan balance at the start of the period.}
#'   \item{interest_loan}{Interest paid on the loan during the period.}
#'   \item{sinking_deposit}{Deposit made into the sinking fund.}
#'   \item{fund_balance_start}{Fund balance at the start of the period.}
#'   \item{interest_fund}{Interest earned by the fund during the period.}
#'   \item{fund_balance_end_before_redemption}{Fund balance before final redemption.}
#'   \item{redemption_from_fund}{Amount withdrawn from the fund to redeem the loan at maturity.}
#'   \item{fund_balance_end}{Fund balance after redemption.}
#'   \item{loan_balance_end}{Outstanding loan balance after redemption.}
#'   \item{total_cashflow_borrower}{Borrower's external cash outflow during the period.}
#'   \item{i_effective_loan_annual}{Equivalent annual effective loan rate.}
#'   \item{i_effective_fund_annual}{Equivalent annual effective fund rate.}
#'   \item{i_loan_period}{Effective loan rate per schedule period.}
#'   \item{i_fund_period}{Effective fund rate per schedule period.}
#'   \item{k}{Schedule frequency.}
#' }
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{i} is the loan interest rate, \code{j} is the
#' sinking-fund accumulation rate, \code{k} is the number of schedule periods
#' per year, \code{m} is the conversion frequency for the loan rate, and
#' \code{j_m} is the conversion frequency for the sinking-fund rate.
#'
#' The supplied annual rate specifications are converted to effective annual
#' rates and then to effective rates per schedule period:
#' \deqn{i_p = (1+i_e)^{1/k} - 1,\qquad
#'       j_p = (1+j_e)^{1/k} - 1.}
#'
#' If \code{deposit} is \code{NULL}, the unit-payment accumulation sum is
#' evaluated stably from the original log fund rate. Only an exactly zero
#' rate uses the zero-rate limit \code{principal / n}. The monetary balance
#' tolerance does not determine whether a rate is zero. The existing
#' requirement of finite annual and period equivalents is retained.
#' Non-finite monetary ledger entries cause an explicit error.
#' A user-supplied deposit is not adjusted to force exact redemption: any
#' uncovered loan balance or fund surplus remains in the corresponding final
#' column and triggers the existing mismatch warning when it exceeds
#' \code{tol}. Redemption from the fund is an internal transfer and is not
#' counted again in the borrower's external cash flow.
#' Storage grows with \code{n}; feasible terms are limited by integer indexing
#' and available memory.
#'
#' Otherwise, the standard sinking-fund formula is used:
#' \deqn{deposit = \frac{principal}{s_{\overline{n}|j_p}}}
#' where
#' \deqn{s_{\overline{n}|j_p} =
#' \frac{(1+j_p)^n - 1}{j_p}.}
#'
#' @seealso \code{\link{amort_schedule}}, \code{\link{s_angle}}
#'
#' @family amortization
#'
#' @examples
#' sinking_fund_schedule(
#'   principal = 100000,
#'   n = 12,
#'   i = 0.12,
#'   j = 0.096,
#'   i_type = "nominal_interest",
#'   j_type = "nominal_interest",
#'   m = 12,
#'   j_m = 12,
#'   k = 12
#' )
#'
#' sinking_fund_schedule(
#'   principal = 50000,
#'   n = 10,
#'   i = 0.02,
#'   j = 0,
#'   deposit = NULL
#' )
#'
#' @export
sinking_fund_schedule <- function(
    principal,
    n,
    i,
    j,
    k = 1L,
    i_type = "effective",
    j_type = "effective",
    m = 1L,
    j_m = 1L,
    deposit = NULL,
    tol = 1e-8
) {
  if (!is.numeric(principal) || length(principal) != 1L || is.na(principal) ||
      !is.finite(principal) || principal <= 0) {
    stop("`principal` must be a single finite positive number.", call. = FALSE)
  }

  if (!is.numeric(n) || length(n) != 1L || is.na(n) || !is.finite(n) ||
      n <= 0 || n != floor(n)) {
    stop("`n` must be a single positive integer.", call. = FALSE)
  }
  n <- as.integer(n)

  if (missing(i) || !is.numeric(i) || length(i) != 1L || is.na(i) ||
      !is.finite(i)) {
    stop("`i` must be a single finite numeric value.", call. = FALSE)
  }

  if (missing(j) || !is.numeric(j) || length(j) != 1L || is.na(j) ||
      !is.finite(j)) {
    stop("`j` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.character(i_type) || length(i_type) != 1L || is.na(i_type)) {
    stop("`i_type` must be a single character string.", call. = FALSE)
  }

  if (!is.character(j_type) || length(j_type) != 1L || is.na(j_type)) {
    stop("`j_type` must be a single character string.", call. = FALSE)
  }

  valid_type <- c(
    "effective",
    "nominal_interest",
    "nominal_discount",
    "force"
  )

  if (!i_type %in% valid_type) {
    stop(
      "`i_type` must be one of: ",
      paste(sprintf("'%s'", valid_type), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!j_type %in% valid_type) {
    stop(
      "`j_type` must be one of: ",
      paste(sprintf("'%s'", valid_type), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!is.numeric(k) || length(k) != 1L || is.na(k) || !is.finite(k) ||
      k <= 0 || k != floor(k)) {
    stop("`k` must be a single positive integer.", call. = FALSE)
  }
  k <- as.integer(k)

  if (!is.numeric(m) || length(m) != 1L || is.na(m) || !is.finite(m) ||
      m <= 0 || m != floor(m)) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }
  m <- as.integer(m)

  if (!is.numeric(j_m) || length(j_m) != 1L || is.na(j_m) ||
      !is.finite(j_m) || j_m <= 0 || j_m != floor(j_m)) {
    stop("`j_m` must be a single positive integer.", call. = FALSE)
  }
  j_m <- as.integer(j_m)

  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a single finite positive numeric value.", call. = FALSE)
  }

  # --- Rate conversion ---
  loan_basis <- .interest_basis(
    i_type = i_type,
    i = i,
    m = m
  )

  fund_basis <- .interest_basis(
    i_type = j_type,
    i = j,
    m = j_m
  )
  i_effective_loan_annual <- loan_basis$effective
  i_effective_fund_annual <- fund_basis$effective

  if (!is.finite(i_effective_loan_annual) || i_effective_loan_annual <= -1) {
    stop(
      "The equivalent annual effective loan rate must be finite and greater than -1.",
      call. = FALSE
    )
  }

  if (!is.finite(i_effective_fund_annual) || i_effective_fund_annual <= -1) {
    stop(
      "The equivalent annual effective sinking-fund rate must be finite and greater than -1.",
      call. = FALSE
    )
  }

  i_loan_period <- expm1(loan_basis$growth / k)
  fund_growth <- fund_basis$growth / k
  i_fund_period <- expm1(fund_growth)

  if (!is.finite(i_loan_period) || i_loan_period <= -1) {
    stop("The effective loan rate per period must be finite and greater than -1.",
         call. = FALSE)
  }

  if (!is.finite(i_fund_period) || i_fund_period <= -1) {
    stop("The effective fund rate per period must be finite and greater than -1.",
         call. = FALSE)
  }

  # --- Compute deposit if not provided ---
  if (is.null(deposit)) {
    deposit <- .cash_flow_value(principal, -.annuity_log_sum(fund_growth, n))
  }

  if (!is.numeric(deposit) || length(deposit) != 1L || is.na(deposit) ||
      !is.finite(deposit) || deposit <= 0) {
    stop("`deposit` must be a single finite positive number.", call. = FALSE)
  }

  # --- Pre-allocate vectors ---
  v_period              <- seq_len(n)
  v_loan_bal_start      <- rep(principal, n)
  v_interest_loan       <- rep(principal * i_loan_period, n)
  v_sinking_deposit     <- rep(deposit, n)
  v_fund_bal_start      <- numeric(n)
  v_interest_fund       <- numeric(n)
  v_fund_bal_end_before <- numeric(n)
  v_redemption          <- numeric(n)
  v_fund_bal_end        <- numeric(n)
  v_loan_bal_end        <- rep(principal, n)
  v_total_cf            <- rep(principal * i_loan_period + deposit, n)

  if (any(!is.finite(c(v_interest_loan[1L], v_total_cf[1L])))) {
    stop("The sinking-fund schedule exceeds the finite monetary range.", call. = FALSE)
  }

  fund_balance <- 0

  for (period_idx in seq_len(n)) {
    v_fund_bal_start[period_idx] <- fund_balance
    v_interest_fund[period_idx] <- fund_balance * i_fund_period
    v_fund_bal_end_before[period_idx] <-
      fund_balance + deposit + v_interest_fund[period_idx]

    if (period_idx == n) {
      v_redemption[period_idx] <- min(v_fund_bal_end_before[period_idx], principal)
      v_loan_bal_end[period_idx] <- principal - v_redemption[period_idx]
    }

    v_fund_bal_end[period_idx] <-
      v_fund_bal_end_before[period_idx] - v_redemption[period_idx]

    if (any(!is.finite(c(v_interest_fund[period_idx], v_fund_bal_end_before[period_idx],
                        v_fund_bal_end[period_idx], v_redemption[period_idx],
                        v_loan_bal_end[period_idx])))) {
      stop("The sinking-fund schedule exceeds the finite monetary range.", call. = FALSE)
    }
    fund_balance <- v_fund_bal_end[period_idx]
  }

  out <- tibble::tibble(
    period = v_period,
    loan_balance_start = v_loan_bal_start,
    interest_loan = v_interest_loan,
    sinking_deposit = v_sinking_deposit,
    fund_balance_start = v_fund_bal_start,
    interest_fund = v_interest_fund,
    fund_balance_end_before_redemption = v_fund_bal_end_before,
    redemption_from_fund = v_redemption,
    fund_balance_end = v_fund_bal_end,
    loan_balance_end = v_loan_bal_end,
    total_cashflow_borrower = v_total_cf,
    i_effective_loan_annual = i_effective_loan_annual,
    i_effective_fund_annual = i_effective_fund_annual,
    i_loan_period = i_loan_period,
    i_fund_period = i_fund_period,
    k = k
  )

  if (abs(v_loan_bal_end[n]) > tol || abs(v_fund_bal_end[n]) > tol) {
    warning(
      sprintf(
        paste0(
          "Final mismatch detected: loan_balance_end = %.10f, ",
          "fund_balance_end = %.10f (tol = %.10f)."
        ),
        v_loan_bal_end[n],
        v_fund_bal_end[n],
        tol
      ),
      call. = FALSE
    )
  }

  out
}
