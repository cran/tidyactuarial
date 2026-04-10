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
#' \deqn{v(t) = (1+i)^{-t}}{v(t) = (1+i)^(-t)}
#'
#' @param time Numeric vector of times in years.
#'   Each value must be greater than or equal to 0.
#' @param spot_rate Numeric vector of spot-rate values.
#' @param spot_type Character vector indicating the spot-rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param spot_m Positive integer vector giving the compounding frequency for
#'   nominal spot-rate inputs.
#'
#' @return Numeric vector of discount factors.
#'
#' @details
#' If \eqn{t = 0}, the discount factor is 1.
#'
#' Input vectors must have length 1 or a common length.
#' Missing values are propagated.
#'
#' @seealso \code{\link{discount_factor_spot_tbl}},
#'   \code{\link{standardize_interest}}, \code{\link{present_value}}
#'
#' @family interest
#'
#' @examples
#' # Simple scalar example
#' discount_factor_spot(
#'   time = 3,
#'   spot_rate = 0.05,
#'   spot_type = "effective"
#' )
#'
#' # Vectorized example
#' discount_factor_spot(
#'   time = c(1, 2, 3),
#'   spot_rate = c(0.05, 0.055, 0.06),
#'   spot_type = "effective"
#' )
#'
#' # FM-style input with nominal annual interest
#' discount_factor_spot(
#'   time = 2,
#'   spot_rate = 0.08,
#'   spot_type = "nominal_interest",
#'   spot_m = 2
#' )
#'
#' @export
discount_factor_spot <- function(
    time,
    spot_rate,
    spot_type = "effective",
    spot_m = 1L
) {
  # --- Early type validation ---
  if (!is.numeric(time)) {
    stop("`time` must be a numeric vector.", call. = FALSE)
  }
  if (!is.numeric(spot_rate)) {
    stop("`spot_rate` must be a numeric vector.", call. = FALSE)
  }
  if (!is.character(spot_type)) {
    stop("`spot_type` must be a character vector.", call. = FALSE)
  }
  if (!is.numeric(spot_m)) {
    stop("`spot_m` must be numeric.", call. = FALSE)
  }

  # --- Determine common size and validate lengths ---
  size <- max(length(time), length(spot_rate), length(spot_type), length(spot_m))

  valid_size <- function(x) length(x) %in% c(1L, size)

  if (!valid_size(time) || !valid_size(spot_rate) ||
      !valid_size(spot_type) || !valid_size(spot_m)) {
    stop(
      "`time`, `spot_rate`, `spot_type`, and `spot_m` must have length 1 ",
      "or a common length.",
      call. = FALSE
    )
  }

  # --- Recycle ---
  time      <- rep_len(time, size)
  spot_rate <- rep_len(spot_rate, size)
  spot_type <- rep_len(spot_type, size)
  spot_m    <- rep_len(spot_m, size)

  # --- Value-level validation ---
  bad_time <- !is.na(time) & (!is.finite(time) | time < 0)
  if (any(bad_time)) {
    stop(
      "`time` must contain only finite values greater than or equal to 0 or NA.",
      call. = FALSE
    )
  }

  bad_m <- !is.na(spot_m) & (!is.finite(spot_m) | spot_m < 1 | spot_m != floor(spot_m))
  if (any(bad_m)) {
    stop(
      "`spot_m` must contain only positive integers or NA.",
      call. = FALSE
    )
  }

  i_effective <- standardize_interest(
    type = spot_type,
    rate = spot_rate,
    m = spot_m
  )

  (1 + i_effective)^(-time)
}


#' Spot discount factor in tibble form
#'
#' Computes the discount factor implied by a spot rate for a given time and
#' returns a tibble with the main input values, the standardized annual
#' effective rate, and the discount factor.
#'
#' This is a reporting wrapper around \code{\link{discount_factor_spot}}.
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
#' \deqn{v(t) = (1+i)^{-t}}{v(t) = (1+i)^(-t)}
#'
#' @param time Numeric vector of times in years.
#'   Each value must be greater than or equal to 0.
#' @param spot_rate Numeric vector of spot-rate values.
#' @param spot_type Character vector indicating the spot-rate type.
#'   Allowed values are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}.
#' @param spot_m Positive integer vector giving the compounding frequency for
#'   nominal spot-rate inputs.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{time}{Time in years.}
#'   \item{spot_rate_input}{Original supplied spot rate.}
#'   \item{spot_type}{Type of supplied spot rate.}
#'   \item{spot_m}{Compounding frequency for nominal spot-rate inputs.}
#'   \item{i_effective}{Equivalent annual effective spot rate.}
#'   \item{discount_factor}{Computed discount factor.}
#' }
#'
#' @details
#' If \eqn{t = 0}, the discount factor is 1.
#'
#' Input vectors must have length 1 or a common length.
#' Missing values are propagated.
#'
#' @seealso \code{\link{discount_factor_spot}},
#'   \code{\link{standardize_interest}}
#'
#' @family interest
#'
#' @examples
#' # Simple scalar example
#' discount_factor_spot_tbl(
#'   time = 3,
#'   spot_rate = 0.05,
#'   spot_type = "effective"
#' )
#'
#' # Vectorized example
#' discount_factor_spot_tbl(
#'   time = c(1, 2, 3),
#'   spot_rate = c(0.05, 0.055, 0.06),
#'   spot_type = "effective"
#' )
#'
#' # FM-style input with nominal annual discount
#' discount_factor_spot_tbl(
#'   time = c(1, 2),
#'   spot_rate = c(0.04, 0.05),
#'   spot_type = "nominal_discount",
#'   spot_m = 2
#' )
#'
#' @export
discount_factor_spot_tbl <- function(
    time,
    spot_rate,
    spot_type = "effective",
    spot_m = 1L
) {
  # Delegate validation and computation to the core function
  discount_factor <- discount_factor_spot(
    time = time,
    spot_rate = spot_rate,
    spot_type = spot_type,
    spot_m = spot_m
  )

  # Recycle for tibble construction (must match what core function used)
  size <- max(length(time), length(spot_rate), length(spot_type), length(spot_m))
  time      <- rep_len(time, size)
  spot_rate <- rep_len(spot_rate, size)
  spot_type <- rep_len(spot_type, size)
  spot_m    <- rep_len(spot_m, size)

  i_effective <- standardize_interest(
    type = spot_type,
    rate = spot_rate,
    m = spot_m
  )

  tibble::tibble(
    time = time,
    spot_rate_input = spot_rate,
    spot_type = spot_type,
    spot_m = as.integer(spot_m),
    i_effective = i_effective,
    discount_factor = discount_factor
  )
}
