#' Present value of a single payment
#'
#' Computes the present value of a future payment \code{C} due at time
#' \code{t}, using the annual effective interest rate implied by the
#' supplied rate specification.
#'
#' The present value is computed as
#' \deqn{PV = C v^t = \frac{C}{(1+i)^t}}{PV = C * v^t = C / (1+i)^t}
#' where \eqn{i} is the annual effective interest rate and
#' \eqn{v = (1+i)^{-1}}{v = 1/(1+i)} is the annual discount factor.
#'
#' The input interest rate may be supplied as:
#' \itemize{
#'   \item annual effective interest rate \eqn{i},
#'   \item nominal annual interest rate \eqn{j^{(m)}}{j(m)},
#'   \item nominal annual discount rate \eqn{d^{(m)}}{d(m)},
#'   \item force of interest \eqn{\delta}{delta}.
#' }
#'
#' Internally, all rate specifications are first converted to the
#' equivalent annual effective interest rate using
#' \code{\link{standardize_interest}}.
#'
#' This is the core numeric version of the calculation. It is designed
#' to work naturally with vectors and with \code{dplyr::mutate()}.
#'
#' @param C Numeric vector of future payment amounts.
#' @param rate Numeric vector of rate values.
#' @param type Character vector indicating the rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the compounding frequency
#'   for nominal rates. Ignored for \code{"effective"} and \code{"force"}.
#' @param t Numeric vector of times in years until payment.
#'
#' @return Numeric vector of present values.
#'
#' @details
#' Input vectors must have length 1 or a common length. Standard recycling
#' is supported only under that rule.
#'
#' Missing values are propagated. This function does not accept dates.
#' If you need a tabular output with actuarial fields, use
#' \code{\link{present_value_tbl}}.
#'
#' @seealso \code{\link{standardize_interest}}, \code{\link{future_value}},
#'   \code{\link{present_value_tbl}}
#'
#' @family time-value
#'
#' @examples
#' # Simple scalar example
#' present_value(C = 1000, rate = 0.08, type = "effective", t = 3)
#'
#' # Medium vectorized example
#' present_value(
#'   C = c(1000, 2500, 4000),
#'   rate = c(0.08, 0.10, 0.12),
#'   type = c("effective", "nominal_interest", "force"),
#'   m = c(1, 12, 1),
#'   t = c(3, 5, 2)
#' )
#'
#' # Use inside a data pipeline
#' if (requireNamespace("dplyr", quietly = TRUE) &&
#'     requireNamespace("tibble", quietly = TRUE)) {
#'   cashflows <- tibble::tibble(
#'     amount = c(1000, 1500, 2000),
#'     rate   = c(0.08, 0.12, 0.09),
#'     type   = c("effective", "force", "nominal_interest"),
#'     m      = c(1, 1, 4),
#'     t      = c(2, 3, 5)
#'   )
#'
#'   dplyr::mutate(
#'     cashflows,
#'     pv = present_value(C = amount, rate = rate, type = type, m = m, t = t)
#'   )
#' }
#'
#' @export
present_value <- function(
    C,
    rate,
    type = "effective",
    m = 1,
    t
) {
  if (missing(C)) {
    stop("`C` must be provided.", call. = FALSE)
  }
  if (missing(rate)) {
    stop("`rate` must be provided.", call. = FALSE)
  }
  if (missing(t)) {
    stop("`t` must be provided.", call. = FALSE)
  }

  # --- Early type validation ---
  if (!is.numeric(C)) {
    stop("`C` must be a numeric vector.", call. = FALSE)
  }
  if (!is.numeric(rate)) {
    stop("`rate` must be a numeric vector.", call. = FALSE)
  }
  if (!is.character(type)) {
    stop("`type` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(m)) {
    stop("`m` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(t)) {
    stop("`t` must be a numeric vector.", call. = FALSE)
  }

  # --- Determine common size and validate lengths ---
  size <- max(length(C), length(rate), length(type), length(m), length(t))

  valid_size <- function(x) length(x) %in% c(1L, size)

  if (!valid_size(C) || !valid_size(rate) || !valid_size(type) ||
      !valid_size(m) || !valid_size(t)) {
    stop(
      "`C`, `rate`, `type`, `m`, and `t` must have length 1 or a common length.",
      call. = FALSE
    )
  }

  # --- Recycle ---
  C    <- rep_len(C, size)
  rate <- rep_len(rate, size)
  type <- rep_len(type, size)
  m    <- rep_len(m, size)
  t    <- rep_len(t, size)

  # --- Value-level validation ---
  bad_C <- !is.na(C) & !is.finite(C)
  if (any(bad_C)) {
    stop("`C` must contain only finite numeric values or NA.", call. = FALSE)
  }

  bad_t <- !is.na(t) & (!is.finite(t) | t < 0)
  if (any(bad_t)) {
    stop(
      "`t` must contain only finite values greater than or equal to 0.",
      call. = FALSE
    )
  }

  # --- Compute ---
  i <- standardize_interest(type = type, rate = rate, m = m)

  out <- rep(NA_real_, size)
  ok <- !is.na(C) & !is.na(t) & !is.na(i)

  out[ok] <- C[ok] / (1 + i[ok])^t[ok]

  out
}


#' Present value details in tibble form
#'
#' Computes the present value of a future payment and returns a tibble
#' containing both the inputs and the actuarial quantities used in the
#' calculation.
#'
#' This is a reporting wrapper around \code{\link{present_value}}. It is useful
#' when you want a tabular summary including the input rate specification,
#' the annual effective rate, the annual discount factor, and the final
#' present value.
#'
#' @param C Numeric vector of future payment amounts.
#' @param rate Numeric vector of rate values.
#' @param type Character vector indicating the rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the compounding frequency
#'   for nominal rates. Ignored for \code{"effective"} and \code{"force"}.
#' @param t Numeric vector of times in years until payment.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{C}{Future payment amount.}
#'   \item{t}{Time in years until payment.}
#'   \item{rate_input}{Original supplied rate.}
#'   \item{rate_type}{Type of supplied rate.}
#'   \item{m}{Compounding frequency.}
#'   \item{i_effective}{Equivalent annual effective interest rate.}
#'   \item{v}{Annual discount factor \eqn{(1+i)^{-1}}{1/(1+i)}.}
#'   \item{present_value}{Computed present value.}
#' }
#'
#' @details
#' This function follows the same recycling and validation rules as
#' \code{\link{present_value}}. Input vectors must have length 1 or a
#' common length.
#'
#' Missing values are propagated.
#'
#' @seealso \code{\link{present_value}}, \code{\link{standardize_interest}},
#'   \code{\link{future_value}}
#'
#' @family time-value
#'
#' @examples
#' # Simple scalar example
#' present_value_tbl(C = 1000, rate = 0.08, type = "effective", t = 3)
#'
#' # Medium vectorized example
#' present_value_tbl(
#'   C = c(1000, 2500, 4000),
#'   rate = c(0.08, 0.10, 0.12),
#'   type = c("effective", "nominal_interest", "force"),
#'   m = c(1, 12, 1),
#'   t = c(3, 5, 2)
#' )
#'
#' # Combine with dplyr
#' if (requireNamespace("dplyr", quietly = TRUE) &&
#'     requireNamespace("tibble", quietly = TRUE)) {
#'   cashflows <- tibble::tibble(
#'     amount = c(1000, 1500, 2000),
#'     rate   = c(0.08, 0.12, 0.09),
#'     type   = c("effective", "force", "nominal_interest"),
#'     m      = c(1, 1, 4),
#'     t      = c(2, 3, 5)
#'   )
#'
#'   present_value_tbl(
#'     C = cashflows$amount,
#'     rate = cashflows$rate,
#'     type = cashflows$type,
#'     m = cashflows$m,
#'     t = cashflows$t
#'   )
#' }
#'
#' @export
present_value_tbl <- function(
    C,
    rate,
    type = "effective",
    m = 1,
    t
) {
  if (missing(C)) {
    stop("`C` must be provided.", call. = FALSE)
  }
  if (missing(rate)) {
    stop("`rate` must be provided.", call. = FALSE)
  }
  if (missing(t)) {
    stop("`t` must be provided.", call. = FALSE)
  }

  # Delegate computation (includes all validation)
  pv <- present_value(C = C, rate = rate, type = type, m = m, t = t)

  # Recycle for tibble construction
  size <- max(length(C), length(rate), length(type), length(m), length(t))
  C    <- rep_len(C, size)
  rate <- rep_len(rate, size)
  type <- rep_len(type, size)
  m    <- rep_len(m, size)
  t    <- rep_len(t, size)

  i_effective <- standardize_interest(type = type, rate = rate, m = m)
  v <- 1 / (1 + i_effective)

  tibble::tibble(
    C = C,
    t = t,
    rate_input = rate,
    rate_type = type,
    m = m,
    i_effective = i_effective,
    v = v,
    present_value = pv
  )
}
