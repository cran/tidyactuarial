#' Future value of a single payment
#'
#' Computes the future value of a payment \code{C} invested at time 0
#' and accumulated to time \code{t}, using the annual effective interest
#' rate implied by the supplied rate specification.
#'
#' The future value is computed as
#' \deqn{FV = C (1+i)^t}{FV = C * (1+i)^t}
#' where \eqn{i} is the annual effective interest rate.
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
#' @param C Numeric vector of initial payment amounts.
#' @param rate Numeric vector of rate values.
#' @param type Character vector indicating the rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the compounding frequency
#'   for nominal rates. Ignored for \code{"effective"} and \code{"force"}.
#' @param t Numeric vector of times in years from valuation to accumulation.
#'
#' @return Numeric vector of future values.
#'
#' @details
#' Input vectors must have length 1 or a common length. Standard recycling
#' is supported only under that rule.
#'
#' Missing values are propagated. This function does not accept dates.
#' If you need a tabular output with actuarial fields, use
#' \code{\link{future_value_tbl}}.
#'
#' @seealso \code{\link{standardize_interest}}, \code{\link{present_value}},
#'   \code{\link{future_value_tbl}}
#'
#' @family time-value
#'
#' @examples
#' # Simple scalar example
#' future_value(C = 1000, rate = 0.08, type = "effective", t = 3)
#'
#' # Medium vectorized example
#' future_value(
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
#'   investments <- tibble::tibble(
#'     deposit = c(1000, 1500, 2000),
#'     rate    = c(0.08, 0.12, 0.09),
#'     type    = c("effective", "force", "nominal_interest"),
#'     m       = c(1, 1, 4),
#'     t       = c(2, 3, 5)
#'   )
#'
#'   dplyr::mutate(
#'     investments,
#'     fv = future_value(C = deposit, rate = rate, type = type, m = m, t = t)
#'   )
#' }
#'
#' @export
future_value <- function(
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

  out[ok] <- C[ok] * (1 + i[ok])^t[ok]

  out
}


#' Future value details in tibble form
#'
#' Computes the future value of a payment and returns a tibble containing
#' both the inputs and the actuarial quantities used in the calculation.
#'
#' This is a reporting wrapper around \code{\link{future_value}}. It is useful
#' when you want a tabular summary including the input rate specification,
#' the annual effective rate, and the final future value.
#'
#' @param C Numeric vector of initial payment amounts.
#' @param rate Numeric vector of rate values.
#' @param type Character vector indicating the rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the compounding frequency
#'   for nominal rates. Ignored for \code{"effective"} and \code{"force"}.
#' @param t Numeric vector of times in years from valuation to accumulation.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{C}{Initial payment amount.}
#'   \item{t}{Time in years to accumulation.}
#'   \item{rate_input}{Original supplied rate.}
#'   \item{rate_type}{Type of supplied rate.}
#'   \item{m}{Compounding frequency.}
#'   \item{i_effective}{Equivalent annual effective interest rate.}
#'   \item{future_value}{Computed future value.}
#' }
#'
#' @details
#' This function follows the same recycling and validation rules as
#' \code{\link{future_value}}. Input vectors must have length 1 or a
#' common length.
#'
#' Missing values are propagated.
#'
#' @seealso \code{\link{future_value}}, \code{\link{standardize_interest}},
#'   \code{\link{present_value}}
#'
#' @family time-value
#'
#' @examples
#' # Simple scalar example
#' future_value_tbl(C = 1000, rate = 0.08, type = "effective", t = 3)
#'
#' # Medium vectorized example
#' future_value_tbl(
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
#'   investments <- tibble::tibble(
#'     deposit = c(1000, 1500, 2000),
#'     rate    = c(0.08, 0.12, 0.09),
#'     type    = c("effective", "force", "nominal_interest"),
#'     m       = c(1, 1, 4),
#'     t       = c(2, 3, 5)
#'   )
#'
#'   future_value_tbl(
#'     C = investments$deposit,
#'     rate = investments$rate,
#'     type = investments$type,
#'     m = investments$m,
#'     t = investments$t
#'   )
#' }
#'
#' @export
future_value_tbl <- function(
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
  fv <- future_value(C = C, rate = rate, type = type, m = m, t = t)

  # Recycle for tibble construction
  size <- max(length(C), length(rate), length(type), length(m), length(t))
  C    <- rep_len(C, size)
  rate <- rep_len(rate, size)
  type <- rep_len(type, size)
  m    <- rep_len(m, size)
  t    <- rep_len(t, size)

  i_effective <- standardize_interest(type = type, rate = rate, m = m)

  tibble::tibble(
    C = C,
    t = t,
    rate_input = rate,
    rate_type = type,
    m = m,
    i_effective = i_effective,
    future_value = fv
  )
}
