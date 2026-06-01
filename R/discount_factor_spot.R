#' Spot discount factor
#'
#' Computes the discount factor implied by a spot rate for a given time, using
#' compact actuarial notation.
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
#' \deqn{v(t) = (1+i)^{-t}.}
#'
#' @param t Numeric vector of times in years. Each value must be greater than
#'   or equal to 0.
#' @param i Numeric vector of spot-rate values.
#' @param i_type Character vector indicating the spot-rate type. Allowed values
#'   are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param m Positive integer vector giving the conversion frequency for nominal
#'   spot-rate inputs.
#' @param tidy Logical scalar. If \code{FALSE}, returns a numeric discount
#'   factor. If \code{TRUE}, returns a tibble with intermediate calculations.
#'
#' @return
#' If \code{tidy = FALSE}, a numeric vector of discount factors.
#'
#' If \code{tidy = TRUE}, a tibble with input values, standardized rates, and
#' discount factors.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{t} denotes time, \code{i} denotes the spot-rate
#' input, \code{i_type} denotes the interest-rate type, and \code{m} denotes the
#' conversion frequency for nominal rates.
#'
#' If \eqn{t = 0}, the discount factor is 1.
#'
#' Input vectors must have length 1 or a common length. Missing values are
#' propagated.
#'
#' @seealso \code{\link{standardize_interest}}, \code{\link{present_value}},
#'   \code{\link{pv_flow}}
#'
#' @family interest
#'
#' @examples
#' # Numeric discount factor
#' discount_factor_spot(
#'   t = 3,
#'   i = 0.05
#' )
#'
#' # Vectorized example
#' discount_factor_spot(
#'   t = c(1, 2, 3),
#'   i = c(0.05, 0.055, 0.06)
#' )
#'
#' # FM-style input with nominal annual interest
#' discount_factor_spot(
#'   t = 2,
#'   i = 0.08,
#'   i_type = "nominal_interest",
#'   m = 2
#' )
#'
#' # Tibble output for teaching or auditing
#' discount_factor_spot(
#'   t = c(1, 2, 3),
#'   i = c(0.05, 0.055, 0.06),
#'   tidy = TRUE
#' )
#'
#' @export
discount_factor_spot <- function(
    t,
    i,
    i_type = "effective",
    m = 1L,
    tidy = FALSE
) {
  if (!is.logical(tidy) || length(tidy) != 1L || is.na(tidy)) {
    stop("`tidy` must be a logical scalar.", call. = FALSE)
  }

  if (missing(t)) {
    stop("`t` must be provided.", call. = FALSE)
  }
  if (missing(i)) {
    stop("`i` must be provided.", call. = FALSE)
  }

  # --- Early type validation ---
  if (!is.numeric(t)) {
    stop("`t` must be a numeric vector.", call. = FALSE)
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

  # --- Determine common size and validate lengths ---
  size <- max(
    length(t),
    length(i),
    length(i_type),
    length(m),
    1L
  )

  valid_size <- function(x) length(x) %in% c(1L, size)

  if (!valid_size(t) ||
      !valid_size(i) ||
      !valid_size(i_type) ||
      !valid_size(m)) {
    stop(
      "`t`, `i`, `i_type`, and `m` must have length 1 or a common length.",
      call. = FALSE
    )
  }

  # --- Recycle ---
  t      <- rep_len(t, size)
  i      <- rep_len(i, size)
  i_type <- rep_len(i_type, size)
  m      <- rep_len(m, size)

  # --- Value-level validation ---
  bad_t <- !is.na(t) & (!is.finite(t) | t < 0)
  if (any(bad_t)) {
    stop(
      "`t` must contain only finite values greater than or equal to 0 or NA.",
      call. = FALSE
    )
  }

  bad_m <- !is.na(m) & (!is.finite(m) | m < 1 | abs(m - round(m)) > 1e-10)
  if (any(bad_m)) {
    stop("`m` must contain only positive integers or NA.", call. = FALSE)
  }

  # --- Compute ---
  i_effective <- standardize_interest(
    i_type = i_type,
    i = i,
    m = m
  )

  if (any(!is.na(i_effective) & (!is.finite(i_effective) | i_effective <= -1))) {
    stop(
      "The standardized annual effective spot rates must be greater than -1.",
      call. = FALSE
    )
  }

  discount_factor <- (1 + i_effective)^(-t)

  if (!tidy) {
    return(discount_factor)
  }

  tibble::tibble(
    t = t,
    i_input = i,
    i_type = i_type,
    m = as.integer(round(m)),
    i_effective = i_effective,
    discount_factor = discount_factor
  )
}
