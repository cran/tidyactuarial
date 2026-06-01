#' Internal rate of return for a cash flow
#'
#' Computes the internal rate of return (IRR) of a cash flow by finding the
#' annual effective rate that makes its present value equal to zero, using
#' compact actuarial notation.
#'
#' The cash flow is supplied explicitly through \code{cf}. Its timing is given
#' either through \code{t} (in years) or \code{date} (calendar dates). If
#' \code{date} is supplied, the earliest date is treated as time 0.
#'
#' The IRR returned is interpreted as an annual effective rate.
#'
#' @param cf Numeric vector of cash flows.
#' @param t Optional numeric vector of cash-flow times in years.
#' @param date Optional vector of cash-flow dates. If supplied, the earliest
#'   date is treated as time 0.
#' @param m Positive integer used only to report an equivalent nominal annual
#'   interest rate convertible \code{m} times per year.
#' @param interval Numeric vector of length 2 giving the search interval for the
#'   annual effective IRR. Default is \code{c(-0.99, 10)}.
#' @param tol Numeric tolerance passed to \code{\link[stats]{uniroot}}.
#' @param maxiter Maximum number of iterations passed to
#'   \code{\link[stats]{uniroot}}.
#' @param day_count Day-count convention used when \code{date} is supplied.
#'   One of \code{"act/365"} or \code{"act/360"}.
#'
#' @return A one-row tibble with:
#' \describe{
#'   \item{irr}{Estimated IRR as an annual effective rate.}
#'   \item{i_effective_annual}{Same as \code{irr}, reported explicitly.}
#'   \item{j_nominal_interest}{Equivalent nominal annual interest rate
#'     convertible \code{m} times.}
#'   \item{delta}{Equivalent force of interest.}
#'   \item{npv}{Present value at the estimated IRR, close to zero.}
#'   \item{interval_left}{Left endpoint of the search interval.}
#'   \item{interval_right}{Right endpoint of the search interval.}
#'   \item{converged}{Logical flag indicating whether a root was found.}
#'   \item{n_iter}{Number of iterations used by \code{uniroot}.}
#'   \item{n_cashflows}{Length of \code{cf}.}
#'   \item{has_both_signs}{Whether the cash flow has at least one positive
#'     and one negative value.}
#'   \item{n_sign_changes}{Number of sign changes in the nonzero cash-flow
#'     sequence.}
#' }
#'
#' If no sign change is present in the cash flow, or if the NPV does not change
#' sign over \code{interval}, the function returns \code{converged = FALSE} and
#' \code{irr = NA_real_}.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{cf} denotes cash flows, \code{t} denotes time,
#' and \code{m} denotes the conversion frequency used to report the equivalent
#' nominal annual interest rate.
#'
#' The IRR is defined as the rate \eqn{r} satisfying
#' \deqn{\sum_{k} C_k (1+r)^{-t_k} = 0}{sum(C_k * (1+r)^(-t_k)) = 0}
#' where \eqn{C_k} are the cash flows and \eqn{t_k} the corresponding times
#' in years.
#'
#' The root is found using \code{\link[stats]{uniroot}} over the specified
#' \code{interval}. If the NPV does not change sign over the interval, no root
#' can be bracketed and the function returns gracefully with
#' \code{converged = FALSE}.
#'
#' The number of sign changes in the nonzero cash-flow sequence is reported as
#' a diagnostic. If there is exactly one sign change, the IRR is usually unique
#' under the usual ordered cash-flow setting.
#'
#' @seealso \code{\link{pv_flow}}, \code{\link{irr_flow_multi}},
#'   \code{\link{bond_ytm}}
#'
#' @family time-value
#'
#' @examples
#' irr_flow(
#'   cf = c(-1000, 300, 400, 500),
#'   t = c(0, 1, 2, 3)
#' )
#'
#' irr_flow(
#'   cf = c(-1000, 300, 400, 500),
#'   date = as.Date(c("2026-01-01", "2027-01-01", "2028-01-01", "2029-01-01"))
#' )
#'
#' @export
irr_flow <- function(
    cf,
    t = NULL,
    date = NULL,
    m = 1L,
    interval = c(-0.99, 10),
    tol = 1e-10,
    maxiter = 1000,
    day_count = c("act/365", "act/360")
) {
  day_count <- match.arg(day_count)

  if (!is.numeric(cf) || length(cf) < 2L) {
    stop("`cf` must be a numeric vector of length at least 2.", call. = FALSE)
  }
  if (any(is.na(cf)) || any(!is.finite(cf))) {
    stop("`cf` must contain only finite numeric values.", call. = FALSE)
  }

  if (!is.null(t) && !is.null(date)) {
    stop("Provide only one of `t` or `date`, not both.", call. = FALSE)
  }
  if (is.null(t) && is.null(date)) {
    stop("You must provide either `t` or `date`.", call. = FALSE)
  }

  if (!is.null(t)) {
    if (!is.numeric(t) || length(t) != length(cf)) {
      stop("`t` must be numeric and have the same length as `cf`.", call. = FALSE)
    }
    if (any(is.na(t)) || any(!is.finite(t)) || any(t < 0)) {
      stop("`t` must contain only finite values >= 0.", call. = FALSE)
    }
  }

  if (!is.null(date)) {
    date <- as.Date(date)
    if (length(date) != length(cf) || any(is.na(date))) {
      stop("`date` must contain valid dates and have the same length as `cf`.", call. = FALSE)
    }
  }

  if (!is.numeric(interval) || length(interval) != 2L || anyNA(interval)) {
    stop("`interval` must be a numeric vector of length 2.", call. = FALSE)
  }
  interval <- sort(interval)
  if (interval[1] <= -1) {
    stop("The lower bound of `interval` must be greater than -1.", call. = FALSE)
  }

  if (!is.numeric(m) || length(m) != 1L || is.na(m) ||
      !is.finite(m) || m <= 0 || m != floor(m)) {
    stop("`m` must be a positive integer.", call. = FALSE)
  }
  m <- as.integer(m)

  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("`tol` must be a single positive numeric value.", call. = FALSE)
  }

  if (!is.numeric(maxiter) || length(maxiter) != 1L || is.na(maxiter) ||
      !is.finite(maxiter) || maxiter <= 0 || maxiter != floor(maxiter)) {
    stop("`maxiter` must be a positive integer.", call. = FALSE)
  }
  maxiter <- as.integer(maxiter)

  # --- Sign pattern diagnostics ---
  has_pos <- any(cf > 0)
  has_neg <- any(cf < 0)
  has_both <- has_pos && has_neg

  cf_nz <- cf[cf != 0]
  n_sign_changes <- if (length(cf_nz) <= 1L) {
    0L
  } else {
    sum(diff(sign(cf_nz)) != 0)
  }

  # Helper: build the non-converged result tibble
  make_na_result <- function() {
    tibble::tibble(
      irr = NA_real_,
      i_effective_annual = NA_real_,
      j_nominal_interest = NA_real_,
      delta = NA_real_,
      npv = NA_real_,
      interval_left = interval[1L],
      interval_right = interval[2L],
      converged = FALSE,
      n_iter = NA_integer_,
      n_cashflows = length(cf),
      has_both_signs = has_both,
      n_sign_changes = n_sign_changes
    )
  }

  if (!has_both) {
    return(make_na_result())
  }

  # --- NPV as a function of annual effective rate ---
  npv_fun <- function(r) {
    pv_flow(
      cf = cf,
      i = r,
      i_type = "effective",
      t = t,
      date = date,
      day_count = day_count
    )
  }

  f_left <- npv_fun(interval[1L])
  f_right <- npv_fun(interval[2L])

  if (is.na(f_left) || is.na(f_right) || f_left * f_right > 0) {
    return(make_na_result())
  }

  # --- Find root ---
  root <- stats::uniroot(
    f = npv_fun,
    interval = interval,
    tol = tol,
    maxiter = maxiter
  )

  irr_est <- root$root
  npv_est <- root$f.root

  i_annual <- irr_est
  delta <- log1p(i_annual)
  j_nom <- m * ((1 + i_annual)^(1 / m) - 1)

  tibble::tibble(
    irr = irr_est,
    i_effective_annual = i_annual,
    j_nominal_interest = j_nom,
    delta = delta,
    npv = npv_est,
    interval_left = interval[1L],
    interval_right = interval[2L],
    converged = TRUE,
    n_iter = root$iter,
    n_cashflows = length(cf),
    has_both_signs = has_both,
    n_sign_changes = n_sign_changes
  )
}
