#' Future value of a single payment
#'
#' Computes the future value of a payment invested at time 0 and accumulated to
#' a given time, using the annual effective interest rate implied by the
#' supplied interest-rate specification and compact actuarial notation.
#'
#' The future value is computed as
#' \deqn{FV = C(1+i)^t}
#' where \eqn{i} is the annual effective interest rate.
#'
#' The input interest rate may be supplied as:
#' \itemize{
#'   \item annual effective interest rate,
#'   \item nominal annual interest rate,
#'   \item nominal annual discount rate,
#'   \item force of interest.
#' }
#'
#' Internally, all rate specifications are first converted to the equivalent
#' annual effective interest rate using \code{\link{standardize_interest}}.
#'
#' @param C Numeric vector of initial payment amounts or capitals.
#' @param i Numeric vector of interest-rate values.
#' @param i_type Character vector indicating the interest-rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the conversion frequency for nominal
#'   rates. Ignored for \code{"effective"} and \code{"force"}.
#' @param t Numeric vector of times in years from valuation to accumulation.
#' @param tidy Logical scalar. If \code{FALSE}, returns a numeric future value.
#'   If \code{TRUE}, returns a tibble with intermediate calculations.
#'
#' @return
#' If \code{tidy = FALSE}, a numeric vector of future values.
#'
#' If \code{tidy = TRUE}, a tibble with input values, equivalent rates,
#' accumulation factors, and future values.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{C} denotes the initial payment amount or capital,
#' \code{t} denotes time, \code{i} denotes the interest-rate input,
#' \code{i_type} denotes the interest-rate type, and \code{m} denotes the
#' conversion frequency for nominal rates.
#'
#' Input vectors must have length 1 or a common length. Missing values are
#' propagated. This function does not accept dates; use \code{\link{fv_flow}}
#' for dated cash flows.
#'
#' @seealso \code{\link{standardize_interest}}, \code{\link{present_value}},
#'   \code{\link{fv_flow}}
#'
#' @family time-value
#'
#' @examples
#' # Numeric future value
#' future_value(C = 1000, i = 0.08, t = 3)
#'
#' # Nominal interest converted monthly
#' future_value(
#'   C = 1000,
#'   i = 0.12,
#'   i_type = "nominal_interest",
#'   m = 12,
#'   t = 5
#' )
#'
#' # Tibble output for teaching or auditing
#' future_value(
#'   C = 1000,
#'   i = 0.08,
#'   t = 3,
#'   tidy = TRUE
#' )
#'
#' # Vectorized example
#' future_value(
#'   C = c(1000, 2500, 4000),
#'   i = c(0.08, 0.10, 0.12),
#'   i_type = c("effective", "nominal_interest", "force"),
#'   m = c(1, 12, 1),
#'   t = c(3, 5, 2)
#' )
#'
#' @export
future_value <- function(
    C,
    i,
    i_type = "effective",
    m = 1,
    t,
    tidy = FALSE
) {
  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  if (missing(C)) {
    stop("`C` must be provided.", call. = FALSE)
  }
  if (missing(i)) {
    stop("`i` must be provided.", call. = FALSE)
  }
  if (missing(t)) {
    stop("`t` must be provided.", call. = FALSE)
  }

  # --- Early type validation ---
  if (!is.numeric(C)) {
    stop("`C` must be a numeric vector.", call. = FALSE)
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
  if (!is.numeric(t)) {
    stop("`t` must be a numeric vector.", call. = FALSE)
  }

  # --- Determine common size and validate lengths ---
  size <- max(
    length(C),
    length(i),
    length(i_type),
    length(m),
    length(t),
    1L
  )

  valid_size <- function(x) length(x) %in% c(1L, size)

  if (!valid_size(C) ||
      !valid_size(i) ||
      !valid_size(i_type) ||
      !valid_size(m) ||
      !valid_size(t)) {
    stop(
      "`C`, `i`, `i_type`, `m`, and `t` must have length 1 ",
      "or a common length.",
      call. = FALSE
    )
  }

  # --- Recycle ---
  C      <- rep_len(C, size)
  i      <- rep_len(i, size)
  i_type <- rep_len(i_type, size)
  m      <- rep_len(m, size)
  t      <- rep_len(t, size)

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

  bad_m <- !is.na(m) & (!is.finite(m) | m < 1 | abs(m - round(m)) > 1e-10)
  if (any(bad_m)) {
    stop("`m` must contain positive integer values or NA.", call. = FALSE)
  }

  # --- Compute equivalent rates ---
  i_effective <- standardize_interest(i_type = i_type, i = i, m = m)

  if (any(!is.na(i_effective) & (!is.finite(i_effective) | i_effective <= -1))) {
    stop(
      "The standardized annual effective interest rates must be greater than -1.",
      call. = FALSE
    )
  }

  accumulation_factor <- (1 + i_effective)^t

  future_value_out <- rep(NA_real_, size)

  ok <- !is.na(C) &
    !is.na(t) &
    !is.na(i_effective)

  future_value_out[ok] <- C[ok] * accumulation_factor[ok]

  if (!tidy) {
    return(future_value_out)
  }

  tibble::tibble(
    C = C,
    t = t,
    i_input = i,
    i_type = i_type,
    m = m,
    i_effective = i_effective,
    accumulation_factor = accumulation_factor,
    future_value = future_value_out
  )
}
