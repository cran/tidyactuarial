#' Present value of a single payment
#'
#' Computes the present value of a future payment due at a given time, using
#' the annual effective interest rate implied by the supplied rate specification.
#'
#' The present value is computed as
#' \deqn{PV = C v^t = \frac{C}{(1+i)^t}}
#' where \eqn{i} is the annual effective interest rate and
#' \eqn{v = (1+i)^{-1}} is the annual discount factor.
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
#' @param amount Numeric vector of future payment amounts.
#' @param rate Numeric vector of rate values.
#' @param rate_type Character vector indicating the rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the compounding frequency
#'   for nominal rates. Ignored for \code{"effective"} and \code{"force"}.
#' @param time Numeric vector of times in years until payment.
#' @param output Character string. Use \code{"value"} to return a numeric
#'   present value, or \code{"table"} to return a tibble with intermediate
#'   calculations.
#'
#' @return
#' If \code{output = "value"}, a numeric vector of present values.
#'
#' If \code{output = "table"}, a tibble with input values, equivalent rates,
#' discount factors, and present values.
#'
#' @details
#' Input vectors must have length 1 or a common length. Missing values are
#' propagated. This function does not accept dates.
#'
#' @seealso \code{\link{standardize_interest}}, \code{\link{future_value}}
#'
#' @family time-value
#'
#' @examples
#' # Numeric present value
#' present_value(amount = 1000, rate = 0.08, time = 3)
#'
#' # Nominal interest converted monthly
#' present_value(
#'   amount = 1000,
#'   rate = 0.12,
#'   rate_type = "nominal_interest",
#'   m = 12,
#'   time = 5
#' )
#'
#' # Tibble output for teaching or auditing
#' present_value(
#'   amount = 1000,
#'   rate = 0.08,
#'   time = 3,
#'   output = "table"
#' )
#'
#' # Vectorized example
#' present_value(
#'   amount = c(1000, 2500, 4000),
#'   rate = c(0.08, 0.10, 0.12),
#'   rate_type = c("effective", "nominal_interest", "force"),
#'   m = c(1, 12, 1),
#'   time = c(3, 5, 2)
#' )
#'
#' @export
present_value <- function(
    amount,
    rate,
    rate_type = "effective",
    m = 1,
    time,
    output = c("value", "table")
) {
  output <- match.arg(output)

  if (missing(amount)) {
    stop("`amount` must be provided.", call. = FALSE)
  }
  if (missing(rate)) {
    stop("`rate` must be provided.", call. = FALSE)
  }
  if (missing(time)) {
    stop("`time` must be provided.", call. = FALSE)
  }

  # --- Early type validation ---
  if (!is.numeric(amount)) {
    stop("`amount` must be a numeric vector.", call. = FALSE)
  }
  if (!is.numeric(rate)) {
    stop("`rate` must be a numeric vector.", call. = FALSE)
  }
  if (!is.character(rate_type)) {
    stop("`rate_type` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(m)) {
    stop("`m` must be numeric.", call. = FALSE)
  }
  if (!is.numeric(time)) {
    stop("`time` must be a numeric vector.", call. = FALSE)
  }

  # --- Determine common size and validate lengths ---
  size <- max(
    length(amount),
    length(rate),
    length(rate_type),
    length(m),
    length(time),
    1L
  )

  valid_size <- function(x) length(x) %in% c(1L, size)

  if (!valid_size(amount) ||
      !valid_size(rate) ||
      !valid_size(rate_type) ||
      !valid_size(m) ||
      !valid_size(time)) {
    stop(
      "`amount`, `rate`, `rate_type`, `m`, and `time` must have length 1 ",
      "or a common length.",
      call. = FALSE
    )
  }

  # --- Recycle ---
  amount    <- rep_len(amount, size)
  rate      <- rep_len(rate, size)
  rate_type <- rep_len(rate_type, size)
  m         <- rep_len(m, size)
  time      <- rep_len(time, size)

  # --- Value-level validation ---
  bad_amount <- !is.na(amount) & !is.finite(amount)
  if (any(bad_amount)) {
    stop("`amount` must contain only finite numeric values or NA.", call. = FALSE)
  }

  bad_time <- !is.na(time) & (!is.finite(time) | time < 0)
  if (any(bad_time)) {
    stop(
      "`time` must contain only finite values greater than or equal to 0.",
      call. = FALSE
    )
  }

  # --- Compute equivalent rates ---
  i_effective <- standardize_interest(type = rate_type, rate = rate, m = m)
  v <- 1 / (1 + i_effective)

  present_value_out <- rep(NA_real_, size)

  ok <- !is.na(amount) &
    !is.na(time) &
    !is.na(i_effective)

  present_value_out[ok] <- amount[ok] / (1 + i_effective[ok])^time[ok]

  if (output == "value") {
    return(present_value_out)
  }

  tibble::tibble(
    amount = amount,
    time = time,
    rate_input = rate,
    rate_type = rate_type,
    m = m,
    i_effective = i_effective,
    v = v,
    present_value = present_value_out
  )
}
