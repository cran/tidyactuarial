#' Standardize an interest rate to the annual effective rate i
#'
#' Converts common interest-rate specifications to the equivalent
#' annual effective interest rate \code{i}.
#'
#' @param type Character vector indicating the rate type.
#'   Must be one of \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, or \code{"force"}.
#' @param rate Numeric vector of rate values.
#' @param m Positive integer vector giving the compounding frequency
#'   for nominal rates. Ignored for \code{"effective"} and \code{"force"}.
#'
#' @return Numeric vector of annual effective rates.
#'   Missing values are propagated.
#'
#' @details
#' The conversion formulas are:
#' \describe{
#'   \item{effective:}{\eqn{i = \text{rate}}{i = rate} (identity)}
#'   \item{nominal_interest:}{\eqn{i = (1 + j^{(m)}/m)^m - 1}{i = (1 + j(m)/m)^m - 1}}
#'   \item{nominal_discount:}{\eqn{i = (1 - d^{(m)}/m)^{-m} - 1}{i = (1 - d(m)/m)^(-m) - 1}}
#'   \item{force:}{\eqn{i = e^{\delta} - 1}{i = exp(delta) - 1}}
#' }
#'
#' Input vectors must have length 1 or a common length.
#'
#' @seealso \code{\link{interest_equivalents}},
#'   \code{\link{discount_factor_spot}}
#'
#' @family interest
#'
#' @examples
#' # Scalar cases
#' standardize_interest("nominal_interest", rate = 0.18, m = 4)
#' standardize_interest("nominal_discount", rate = 0.10, m = 12)
#' standardize_interest("force", rate = 0.12)
#'
#' # Vectorized case
#' standardize_interest(
#'   type = c("nominal_interest", "force", "effective"),
#'   rate = c(0.06, 0.05, 0.04),
#'   m = c(12, 1, 1)
#' )
#'
#' # Use inside a data pipeline
#' if (requireNamespace("dplyr", quietly = TRUE) &&
#'     requireNamespace("tibble", quietly = TRUE)) {
#'   portfolio <- tibble::tibble(
#'     policy_id  = 1:3,
#'     gross_rate = c(0.05, 0.08, 0.10),
#'     rate_type  = c("force", "nominal_interest", "nominal_discount"),
#'     frequency  = c(NA, 4, 12)
#'   )
#'
#'   dplyr::mutate(
#'     portfolio,
#'     effective_rate = standardize_interest(
#'       type = rate_type,
#'       rate = gross_rate,
#'       m = frequency
#'     )
#'   )
#' }
#'
#' @export
standardize_interest <- function(
    type = "effective",
    rate,
    m = 1
) {
  if (missing(rate)) {
    stop("`rate` must be provided.", call. = FALSE)
  }

  if (!is.character(type)) {
    stop("`type` must be a character vector.", call. = FALSE)
  }

  if (!is.numeric(rate)) {
    stop("`rate` must be a numeric vector.", call. = FALSE)
  }

  if (!is.numeric(m)) {
    stop("`m` must be a numeric vector.", call. = FALSE)
  }

  size <- max(length(type), length(rate), length(m))

  valid_size <- function(x) length(x) %in% c(1L, size)

  if (!valid_size(type) || !valid_size(rate) || !valid_size(m)) {
    stop(
      "`type`, `rate`, and `m` must have length 1 or a common length.",
      call. = FALSE
    )
  }

  type <- rep_len(tolower(type), size)
  rate <- rep_len(rate, size)
  m    <- rep_len(m, size)

  valid_types <- c(
    "effective",
    "nominal_interest",
    "nominal_discount",
    "force"
  )

  bad_type <- !is.na(type) & !(type %in% valid_types)
  if (any(bad_type)) {
    stop(
      "`type` must contain only 'effective', 'nominal_interest', ",
      "'nominal_discount', or 'force'.",
      call. = FALSE
    )
  }

  bad_rate <- !is.na(rate) & !is.finite(rate)
  if (any(bad_rate)) {
    stop("`rate` must contain only finite numeric values or NA.", call. = FALSE)
  }

  idx_nom_i <- type == "nominal_interest"
  idx_nom_d <- type == "nominal_discount"
  idx_eff   <- type == "effective"
  idx_force <- type == "force"
  idx_nom   <- idx_nom_i | idx_nom_d

  missing_m <- idx_nom & is.na(m)
  if (any(missing_m)) {
    stop("`m` must be provided for nominal rates.", call. = FALSE)
  }

  bad_m <- idx_nom & !(is.finite(m) & m > 0 & m == floor(m))
  if (any(bad_m)) {
    stop(
      "`m` must be a positive integer for nominal rates.",
      call. = FALSE
    )
  }

  bad_eff <- idx_eff & !is.na(rate) & rate <= -1
  if (any(bad_eff)) {
    stop(
      "For `type = 'effective'`, `rate` must be greater than -1.",
      call. = FALSE
    )
  }

  bad_nom_i <- idx_nom_i & !is.na(rate) & (1 + rate / m) <= 0
  if (any(bad_nom_i)) {
    stop(
      "For `type = 'nominal_interest'`, `1 + rate / m` must be positive.",
      call. = FALSE
    )
  }

  bad_nom_d <- idx_nom_d & !is.na(rate) & (1 - rate / m) <= 0
  if (any(bad_nom_d)) {
    stop(
      "For `type = 'nominal_discount'`, `1 - rate / m` must be positive.",
      call. = FALSE
    )
  }

  out <- rep(NA_real_, size)

  ok_eff <- idx_eff & !is.na(rate)
  out[ok_eff] <- rate[ok_eff]

  ok_nom_i <- idx_nom_i & !is.na(rate)
  out[ok_nom_i] <- (1 + rate[ok_nom_i] / m[ok_nom_i])^m[ok_nom_i] - 1

  ok_nom_d <- idx_nom_d & !is.na(rate)
  out[ok_nom_d] <- (1 - rate[ok_nom_d] / m[ok_nom_d])^(-m[ok_nom_d]) - 1

  ok_force <- idx_force & !is.na(rate)
  out[ok_force] <- exp(rate[ok_force]) - 1

  out
}
