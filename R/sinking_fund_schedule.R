#' Sinking fund amortization schedule for a loan
#'
#' Builds a sinking fund schedule under a fixed loan rate and a fixed
#' accumulation rate for the sinking fund.
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
#' @param n Positive integer. Number of periods.
#' @param i_loan Numeric scalar. Effective interest rate per period on the loan.
#' @param i_fund Numeric scalar. Effective interest rate per period on the
#'   sinking fund.
#' @param deposit Optional numeric scalar. Level sinking-fund deposit per
#'   period. If \code{NULL}, it is computed so that the fund accumulates to
#'   \code{principal} at time \code{n}.
#' @param tol Numeric tolerance used for zero checks and final-balance checks.
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
#' }
#'
#' @details
#' If \code{deposit} is \code{NULL} and \code{i_fund} is approximately zero,
#' the deposit is computed as \code{principal / n}.
#'
#' Otherwise, the standard sinking-fund formula is used:
#' \deqn{\text{deposit} = \frac{\text{principal}}{s_n}}{deposit = principal / s_n}
#' where
#' \deqn{s_n = \frac{(1+i_{\text{fund}})^n - 1}{i_{\text{fund}}}}{s_n = ((1+i_fund)^n - 1) / i_fund}
#'
#' @seealso \code{\link{amort_schedule}}, \code{\link{s_angle}}
#'
#' @family amortization
#'
#' @examples
#' sinking_fund_schedule(
#'   principal = 100000,
#'   n = 12,
#'   i_loan = 0.01,
#'   i_fund = 0.008
#' )
#'
#' sinking_fund_schedule(
#'   principal = 50000,
#'   n = 10,
#'   i_loan = 0.02,
#'   i_fund = 0,
#'   deposit = NULL
#' )
#'
#' @export
sinking_fund_schedule <- function(
    principal,
    n,
    i_loan,
    i_fund,
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

  if (!is.numeric(i_loan) || length(i_loan) != 1L || is.na(i_loan) ||
      !is.finite(i_loan) || i_loan <= -1) {
    stop("`i_loan` must be a single finite numeric value greater than -1.", call. = FALSE)
  }

  if (!is.numeric(i_fund) || length(i_fund) != 1L || is.na(i_fund) ||
      !is.finite(i_fund) || i_fund <= -1) {
    stop("`i_fund` must be a single finite numeric value greater than -1.", call. = FALSE)
  }

  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("`tol` must be a single positive numeric value.", call. = FALSE)
  }

  # --- Compute deposit if not provided ---
  if (is.null(deposit)) {
    if (abs(i_fund) < tol) {
      deposit <- principal / n
    } else {
      s_n <- ((1 + i_fund)^n - 1) / i_fund
      deposit <- principal / s_n
    }
  }

  if (!is.numeric(deposit) || length(deposit) != 1L || is.na(deposit) ||
      !is.finite(deposit) || deposit <= 0) {
    stop("`deposit` must be a single finite positive number.", call. = FALSE)
  }

  # --- Pre-allocate vectors ---
  v_period              <- seq_len(n)
  v_loan_bal_start      <- rep(principal, n)
  v_interest_loan       <- rep(principal * i_loan, n)
  v_sinking_deposit     <- rep(deposit, n)
  v_fund_bal_start      <- numeric(n)
  v_interest_fund       <- numeric(n)
  v_fund_bal_end_before <- numeric(n)
  v_redemption          <- numeric(n)
  v_fund_bal_end        <- numeric(n)
  v_loan_bal_end        <- rep(principal, n)
  v_total_cf            <- rep(principal * i_loan + deposit, n)

  fund_balance <- 0

  for (k in seq_len(n)) {
    v_fund_bal_start[k] <- fund_balance
    v_interest_fund[k] <- fund_balance * i_fund
    v_fund_bal_end_before[k] <- fund_balance + deposit + v_interest_fund[k]

    if (k == n) {
      v_redemption[k] <- min(v_fund_bal_end_before[k], principal)
      v_loan_bal_end[k] <- principal - v_redemption[k]
    }

    v_fund_bal_end[k] <- v_fund_bal_end_before[k] - v_redemption[k]
    fund_balance <- v_fund_bal_end[k]
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
    total_cashflow_borrower = v_total_cf
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
