#' Spot discount factor
#'
#' Computes the discount factor implied by a spot rate for a given time.
#'
#' The spot rate may be supplied in FM-style notation:
#' \itemize{
#'   \item annual effective rate,
#'   \item nominal annual interest rate,
#'   \item nominal annual discount rate,
#'   \item force of interest.
#' }
#'
#' Internally, the supplied spot rate is first converted to the equivalent
#' annual effective rate using \code{\link{standardize_interest}}. The discount
#' factor is then computed as
#' \deqn{v(t) = (1+i)^{-t}}
#'
#' @param time Numeric vector of times in years.
#'   Each value must be greater than or equal to 0.
#' @param rate Numeric vector of spot-rate values.
#' @param rate_type Character vector indicating the spot-rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the compounding frequency for
#'   nominal spot-rate inputs.
#' @param output Character string. Use \code{"value"} to return a numeric
#'   discount factor, or \code{"table"} to return a tibble with intermediate
#'   calculations.
#'
#' @return
#' If \code{output = "value"}, a numeric vector of discount factors.
#'
#' If \code{output = "table"}, a tibble with input values, standardized rates,
#' and discount factors.
#'
#' @details
#' If \eqn{t = 0}, the discount factor is 1.
#'
#' Input vectors must have length 1 or a common length.
#' Missing values are propagated.
#'
#' @seealso \code{\link{standardize_interest}}, \code{\link{present_value}}
#'
#' @family interest
#'
#' @examples
#' # Numeric discount factor
#' discount_factor_spot(
#'   time = 3,
#'   rate = 0.05
#' )
#'
#' # Vectorized example
#' discount_factor_spot(
#'   time = c(1, 2, 3),
#'   rate = c(0.05, 0.055, 0.06)
#' )
#'
#' # FM-style input with nominal annual interest
#' discount_factor_spot(
#'   time = 2,
#'   rate = 0.08,
#'   rate_type = "nominal_interest",
#'   m = 2
#' )
#'
#' # Tibble output for teaching or auditing
#' discount_factor_spot(
#'   time = c(1, 2, 3),
#'   rate = c(0.05, 0.055, 0.06),
#'   output = "table"
#' )
#'
#' @export
discount_factor_spot <- function(
    time,
    rate,
    rate_type = "effective",
    m = 1L,
    output = c("value", "table")
) {
  output <- match.arg(output)

  if (missing(time)) {
    stop("`time` must be provided.", call. = FALSE)
  }
  if (missing(rate)) {
    stop("`rate` must be provided.", call. = FALSE)
  }

  # --- Early type validation ---
  if (!is.numeric(time)) {
    stop("`time` must be a numeric vector.", call. = FALSE)
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

  # --- Determine common size and validate lengths ---
  size <- max(
    length(time),
    length(rate),
    length(rate_type),
    length(m),
    1L
  )

  valid_size <- function(x) length(x) %in% c(1L, size)

  if (!valid_size(time) ||
      !valid_size(rate) ||
      !valid_size(rate_type) ||
      !valid_size(m)) {
    stop(
      "`time`, `rate`, `rate_type`, and `m` must have length 1 ",
      "or a common length.",
      call. = FALSE
    )
  }

  # --- Recycle ---
  time      <- rep_len(time, size)
  rate      <- rep_len(rate, size)
  rate_type <- rep_len(rate_type, size)
  m         <- rep_len(m, size)

  # --- Value-level validation ---
  bad_time <- !is.na(time) & (!is.finite(time) | time < 0)
  if (any(bad_time)) {
    stop(
      "`time` must contain only finite values greater than or equal to 0 or NA.",
      call. = FALSE
    )
  }

  bad_m <- !is.na(m) & (!is.finite(m) | m < 1 | m != floor(m))
  if (any(bad_m)) {
    stop("`m` must contain only positive integers or NA.", call. = FALSE)
  }

  # --- Compute ---
  i_effective <- standardize_interest(
    type = rate_type,
    rate = rate,
    m = m
  )

  discount_factor <- (1 + i_effective)^(-time)

  if (output == "value") {
    return(discount_factor)
  }

  tibble::tibble(
    time = time,
    rate_input = rate,
    rate_type = rate_type,
    m = as.integer(m),
    i_effective = i_effective,
    discount_factor = discount_factor
  )
}
