#' General amortization schedule with variable rates and payments
#'
#' Builds a loan amortization schedule allowing the effective interest-rate
#' pattern, regular-payment pattern, and extra-principal pattern to vary by
#' period.
#'
#' The function extends the fixed-pattern workflow of
#' \code{\link{amort_schedule}} without modifying that function. It supports:
#' \itemize{
#'   \item a scalar interest rate recycled over all periods;
#'   \item a vector of interest rates, one for each period;
#'   \item a level payment computed automatically;
#'   \item a scalar or period-specific payment vector;
#'   \item a payment function depending on the current loan state.
#' }
#'
#' @param principal Positive numeric scalar. Initial outstanding balance.
#' @param n Positive integer scalar. Maximum number of amortization periods.
#' @param i Numeric vector of interest-rate inputs. It must have length 1 or
#'   \code{n}. A scalar is recycled over all periods.
#' @param i_type Character vector indicating the interest-rate type. Allowed
#'   values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}. It must have length 1 or
#'   \code{n}.
#' @param m Positive integer vector giving the conversion frequency for nominal
#'   rates. It must have length 1 or \code{n}.
#' @param k Positive integer scalar. Number of amortization periods per year.
#' @param timing Character scalar. Either \code{"immediate"} for payments at
#'   the end of each period or \code{"due"} for payments at the beginning.
#' @param payment Payment specification. It may be:
#'   \itemize{
#'     \item \code{NULL}, in which case a level payment is calculated;
#'     \item a nonnegative numeric scalar;
#'     \item a nonnegative numeric vector of length \code{n};
#'     \item a function with arguments \code{period}, \code{ob_start}, and
#'       \code{i_effective_period}, returning one nonnegative payment.
#'   }
#' @param extra_principal Optional nonnegative extra-principal payments. It may
#'   be a scalar or a numeric vector of length \code{n}. Defaults to zero.
#' @param output Character scalar. Either \code{"schedule"} for the complete
#'   amortization table or \code{"summary"} for a compact actuarial summary.
#' @param tol Positive numeric scalar used for zero-balance detection.
#'
#' @return
#' With \code{output = "schedule"}, a tibble with one row per realized period.
#'
#' With \code{output = "summary"}, a one-row tibble containing the initial
#' principal, number of realized periods, total interest, total paid, ending
#' balance, and an indicator of negative amortization.
#'
#' @details
#' Each value of \code{i} is interpreted as an annual rate under its associated
#' \code{i_type} and \code{m}. It is first converted to an annual effective
#' rate and then to the effective rate for one amortization period:
#' \deqn{i_j^{(p)} = (1 + i_j)^{1/k} - 1.}
#'
#' When \code{payment = NULL}, the level payment is determined from the
#' period-specific discount factors. For payments in arrears,
#' \deqn{P =
#' \frac{L}{\sum_{j=1}^{n}
#' \prod_{h=1}^{j}(1+i_h^{(p)})^{-1}}.}
#'
#' For payments in advance, the discount factors correspond to payment times
#' \eqn{0,1,\ldots,n-1}.
#'
#' A payment smaller than the interest charged may produce negative
#' amortization. The schedule identifies such periods explicitly.
#'
#' @seealso \code{\link{amort_schedule}}, \code{\link{a_angle}},
#'   \code{\link{standardize_interest}}
#'
#' @family amortization
#'
#' @examples
#' # Level payments and a constant rate
#' amort_schedule_general(
#'   principal = 100000,
#'   n = 12,
#'   i = 0.12,
#'   i_type = "nominal_interest",
#'   m = 12,
#'   k = 12
#' )
#'
#' # Increasing payments
#' amort_schedule_general(
#'   principal = 10000,
#'   n = 5,
#'   i = 0.06,
#'   payment = c(1800, 1900, 2000, 2200, 3000)
#' )
#'
#' # Period-specific annual effective rates
#' amort_schedule_general(
#'   principal = 10000,
#'   n = 4,
#'   i = c(0.04, 0.05, 0.06, 0.07),
#'   payment = NULL
#' )
#'
#' # Payment determined from the current balance
#' payment_rule <- function(
#'   period,
#'   ob_start,
#'   i_effective_period
#' ) {
#'   1500 + 100 * (period - 1)
#' }
#'
#' amort_schedule_general(
#'   principal = 8000,
#'   n = 6,
#'   i = 0.05,
#'   payment = payment_rule
#' )
#'
#' # Compact summary
#' amort_schedule_general(
#'   principal = 100000,
#'   n = 60,
#'   i = 0.08,
#'   k = 12,
#'   output = "summary"
#' )
#'
#' @export
amort_schedule_general <- function(
    principal,
    n,
    i,
    i_type = "effective",
    m = 1L,
    k = 1L,
    timing = c("immediate", "due"),
    payment = NULL,
    extra_principal = 0,
    output = c("schedule", "summary"),
    tol = 1e-8
) {
  timing <- match.arg(timing)
  output <- match.arg(output)

  if (!is.numeric(principal) ||
      length(principal) != 1L ||
      is.na(principal) ||
      !is.finite(principal) ||
      principal <= 0) {
    stop(
      "`principal` must be a single finite positive number.",
      call. = FALSE
    )
  }

  if (!is.numeric(n) ||
      length(n) != 1L ||
      is.na(n) ||
      !is.finite(n) ||
      n <= 0 ||
      n != floor(n)) {
    stop(
      "`n` must be a single positive integer.",
      call. = FALSE
    )
  }
  n <- as.integer(n)

  if (missing(i) || !is.numeric(i)) {
    stop(
      "`i` must be a numeric vector of length 1 or `n`.",
      call. = FALSE
    )
  }

  if (!is.character(i_type)) {
    stop(
      "`i_type` must be a character vector of length 1 or `n`.",
      call. = FALSE
    )
  }

  if (!is.numeric(m)) {
    stop(
      "`m` must be a numeric vector of length 1 or `n`.",
      call. = FALSE
    )
  }

  valid_period_length <- function(x) {
    length(x) %in% c(1L, n)
  }

  if (!valid_period_length(i) ||
      !valid_period_length(i_type) ||
      !valid_period_length(m)) {
    stop(
      "`i`, `i_type`, and `m` must have length 1 or `n`.",
      call. = FALSE
    )
  }

  if (!is.numeric(k) ||
      length(k) != 1L ||
      is.na(k) ||
      !is.finite(k) ||
      k <= 0 ||
      k != floor(k)) {
    stop(
      "`k` must be a single positive integer.",
      call. = FALSE
    )
  }
  k <- as.integer(k)

  if (!is.numeric(tol) ||
      length(tol) != 1L ||
      is.na(tol) ||
      !is.finite(tol) ||
      tol <= 0) {
    stop(
      "`tol` must be a single positive finite number.",
      call. = FALSE
    )
  }

  i_input <- rep_len(i, n)
  i_type_input <- rep_len(i_type, n)
  m_input <- rep_len(m, n)

  i_effective_annual <- standardize_interest(
    i_type = i_type_input,
    i = i_input,
    m = m_input
  )

  if (any(is.na(i_effective_annual))) {
    stop(
      "Interest-rate inputs cannot contain missing values.",
      call. = FALSE
    )
  }

  i_effective_period <-
    (1 + i_effective_annual)^(1 / k) - 1

  if (any(!is.finite(i_effective_period)) ||
      any(i_effective_period <= -1)) {
    stop(
      "Each effective rate per period must be finite and greater than -1.",
      call. = FALSE
    )
  }

  extra_vec <- .amort_general_recycle_nonnegative(
    x = extra_principal,
    n = n,
    argument = "extra_principal"
  )

  payment_source <- if (is.null(payment)) {
    "calculated_level"
  } else if (is.function(payment)) {
    "function"
  } else {
    "supplied"
  }

  if (is.null(payment)) {
    period_growth <- cumprod(1 + i_effective_period)

    discount_payments <- if (timing == "immediate") {
      1 / period_growth
    } else {
      c(
        1,
        if (n > 1L) {
          1 / period_growth[seq_len(n - 1L)]
        } else {
          numeric(0)
        }
      )
    }

    payment_vec <- rep(
      principal / sum(discount_payments),
      n
    )
    payment_function <- NULL
  } else if (is.function(payment)) {
    payment_vec <- rep(NA_real_, n)
    payment_function <- payment
  } else {
    payment_vec <- .amort_general_recycle_nonnegative(
      x = payment,
      n = n,
      argument = "payment"
    )
    payment_function <- NULL
  }

  out_period <- integer(n)
  out_ob_start <- numeric(n)
  out_i_annual <- numeric(n)
  out_i_period <- numeric(n)
  out_interest <- numeric(n)
  out_payment <- numeric(n)
  out_extra <- numeric(n)
  out_principal <- numeric(n)
  out_total_principal <- numeric(n)
  out_cashflow <- numeric(n)
  out_ob_end <- numeric(n)
  out_negative <- logical(n)

  ob <- principal
  n_realized <- 0L

  for (period_idx in seq_len(n)) {
    if (ob <= tol) {
      break
    }

    rate_period <- i_effective_period[[period_idx]]

    payment_requested <- if (is.function(payment_function)) {
      value <- payment_function(
        period = period_idx,
        ob_start = ob,
        i_effective_period = rate_period
      )

      if (!is.numeric(value) ||
          length(value) != 1L ||
          is.na(value) ||
          !is.finite(value) ||
          value < 0) {
        stop(
          "`payment` function must return one finite nonnegative number.",
          call. = FALSE
        )
      }

      as.numeric(value)
    } else {
      payment_vec[[period_idx]]
    }

    extra_requested <- extra_vec[[period_idx]]

    if (timing == "immediate") {
      interest_period <- ob * rate_period
      total_due <- ob + interest_period

      payment_used <- min(payment_requested, total_due)
      remaining_due <- total_due - payment_used
      extra_used <- min(extra_requested, remaining_due)

      cashflow_period <- payment_used + extra_used
      principal_regular <- payment_used - interest_period
      total_principal <- principal_regular + extra_used
      ob_new <- total_due - cashflow_period
    } else {
      payment_used <- min(payment_requested, ob)
      remaining_balance <- ob - payment_used
      extra_used <- min(extra_requested, remaining_balance)

      balance_after_payment <-
        remaining_balance - extra_used

      if (balance_after_payment <= tol) {
        interest_period <- 0
        ob_new <- 0
      } else {
        interest_period <-
          balance_after_payment * rate_period
        ob_new <-
          balance_after_payment + interest_period
      }

      cashflow_period <- payment_used + extra_used
      principal_regular <- payment_used
      total_principal <- payment_used + extra_used
    }

    if (abs(ob_new) <= tol) {
      ob_new <- 0
    }

    negative_amortization <-
      total_principal < -tol

    n_realized <- n_realized + 1L

    out_period[[n_realized]] <- period_idx
    out_ob_start[[n_realized]] <- ob
    out_i_annual[[n_realized]] <-
      i_effective_annual[[period_idx]]
    out_i_period[[n_realized]] <- rate_period
    out_interest[[n_realized]] <- interest_period
    out_payment[[n_realized]] <- payment_used
    out_extra[[n_realized]] <- extra_used
    out_principal[[n_realized]] <- principal_regular
    out_total_principal[[n_realized]] <-
      total_principal
    out_cashflow[[n_realized]] <- cashflow_period
    out_ob_end[[n_realized]] <- ob_new
    out_negative[[n_realized]] <-
      negative_amortization

    ob <- ob_new
  }

  idx <- seq_len(n_realized)

  schedule <- tibble::tibble(
    period = out_period[idx],
    ob_start = out_ob_start[idx],
    i_effective_annual = out_i_annual[idx],
    i_effective_period = out_i_period[idx],
    interest = out_interest[idx],
    payment = out_payment[idx],
    extra_principal = out_extra[idx],
    principal = out_principal[idx],
    total_principal = out_total_principal[idx],
    cashflow = out_cashflow[idx],
    ob_end = out_ob_end[idx],
    negative_amortization = out_negative[idx],
    timing = timing,
    k = k,
    payment_source = payment_source
  )

  if (identical(output, "schedule")) {
    return(schedule)
  }

  tibble::tibble(
    principal = principal,
    periods_realized = nrow(schedule),
    total_interest = sum(schedule$interest),
    total_paid = sum(schedule$cashflow),
    ending_balance = if (nrow(schedule) > 0L) {
      schedule$ob_end[[nrow(schedule)]]
    } else {
      principal
    },
    negative_amortization = any(
      schedule$negative_amortization
    )
  )
}


.amort_general_recycle_nonnegative <- function(
    x,
    n,
    argument
) {
  if (!is.numeric(x) ||
      !(length(x) %in% c(1L, n)) ||
      any(is.na(x)) ||
      any(!is.finite(x)) ||
      any(x < 0)) {
    stop(
      sprintf(
        "`%s` must be a finite nonnegative numeric vector of length 1 or `n`.",
        argument
      ),
      call. = FALSE
    )
  }

  rep_len(as.numeric(x), n)
}
