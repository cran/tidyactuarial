#' Equivalent interest rates in FM actuarial notation
#'
#' Converts a single interest-rate specification into equivalent actuarial
#' rates for the same conversion frequency \code{m}.
#'
#' Internally, the supplied rate is first converted to the annual effective
#' interest rate \eqn{i} using \code{\link{standardize_interest}}.
#'
#' @param i_type Character string indicating the input interest-rate type.
#'   Must be one of \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, or \code{"force"}.
#' @param i Numeric scalar giving the interest-rate value.
#' @param m Positive integer scalar giving the conversion frequency for nominal
#'   rates.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{family}{Rate family: \code{"effective"}, \code{"discount"},
#'     \code{"force"}, \code{"nominal_interest"}, or
#'     \code{"nominal_discount"}.}
#'   \item{notation}{Actuarial notation for the equivalent rate.}
#'   \item{m}{Conversion frequency used for nominal rates.}
#'   \item{description}{Human-readable description.}
#'   \item{value}{Equivalent rate value.}
#' }
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{i} denotes the interest-rate value,
#' \code{i_type} denotes the type of interest rate, and \code{m} denotes the
#' conversion frequency for nominal rates.
#'
#' Given the annual effective interest rate \eqn{i}, the equivalents are:
#' \describe{
#'   \item{Effective discount rate:}{\eqn{d = i/(1+i)}{d = i/(1+i)}}
#'   \item{Discount factor:}{\eqn{v = 1/(1+i)}{v = 1/(1+i)}}
#'   \item{Force of interest:}{\eqn{\delta = \ln(1+i)}{delta = ln(1+i)}}
#'   \item{Nominal interest:}{\eqn{j^{(m)} = m[(1+i)^{1/m} - 1]}{j(m) = m*((1+i)^(1/m) - 1)}}
#'   \item{Nominal discount:}{\eqn{d^{(m)} = m[1 - (1+i)^{-1/m}]}{d(m) = m*(1 - (1+i)^(-1/m))}}
#' }
#'
#' @seealso \code{\link{standardize_interest}},
#'   \code{\link{discount_factor_spot}}
#'
#' @family interest
#'
#' @examples
#' interest_equivalents(i_type = "nominal_interest", i = 0.18, m = 4)
#' interest_equivalents(i_type = "nominal_discount", i = 0.10, m = 12)
#' interest_equivalents(i_type = "force", i = 0.12)
#' interest_equivalents(i_type = "effective", i = 0.08)
#'
#' # Batch use with purrr
#' if (requireNamespace("purrr", quietly = TRUE) &&
#'     requireNamespace("tibble", quietly = TRUE)) {
#'   cases <- tibble::tibble(
#'     i_type = c("effective", "force"),
#'     i = c(0.08, 0.12),
#'     m = c(1, 1)
#'   )
#'
#'   purrr::pmap(cases, interest_equivalents)
#' }
#'
#' @export
interest_equivalents <- function(
    i_type = c("effective", "nominal_interest", "nominal_discount", "force"),
    i,
    m = 1L
) {
  i_type <- match.arg(i_type)

  if (missing(i)) {
    stop("`i` must be provided.", call. = FALSE)
  }

  if (!is.numeric(i) || length(i) != 1L || is.na(i) || !is.finite(i)) {
    stop("`i` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.numeric(m) || length(m) != 1L || is.na(m) || !is.finite(m) ||
      m <= 0 || m != floor(m)) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }

  m <- as.integer(m)

  i_eff <- standardize_interest(type = i_type, rate = i, m = m)

  d     <- i_eff / (1 + i_eff)
  v     <- 1 / (1 + i_eff)
  delta <- log1p(i_eff)
  j_m   <- m * expm1(log1p(i_eff) / m)
  d_m   <- m * (1 - exp(-log1p(i_eff) / m))

  freq_label <- switch(
    as.character(m),
    "1"  = "annual",
    "2"  = "semiannual",
    "4"  = "quarterly",
    "12" = "monthly",
    paste0(m, "-thly")
  )

  tibble::tibble(
    family = c(
      "effective",
      "discount",
      "discount_factor",
      "force",
      "nominal_interest",
      "nominal_discount"
    ),
    notation = c(
      "i",
      "d",
      "v",
      "delta",
      paste0("j^(", m, ")"),
      paste0("d^(", m, ")")
    ),
    m = c(NA_integer_, NA_integer_, NA_integer_, NA_integer_, m, m),
    description = c(
      "effective annual interest rate",
      "effective annual discount rate",
      "annual discount factor",
      "force of interest",
      paste0("nominal annual interest rate convertible ", freq_label),
      paste0("nominal annual discount rate convertible ", freq_label)
    ),
    value = c(i_eff, d, v, delta, j_m, d_m)
  )
}
