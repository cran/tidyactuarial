#' Equivalent interest rates in FM actuarial notation
#'
#' Converts a single interest-rate specification into equivalent
#' actuarial rates for the same compounding frequency \code{m}.
#'
#' Internally, the supplied rate is first converted to the annual
#' effective interest rate \eqn{i} using \code{\link{standardize_interest}}.
#'
#' @param type Character string indicating the input rate type.
#'   Must be one of \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, or \code{"force"}.
#' @param rate Numeric scalar giving the rate value.
#' @param m Positive integer scalar giving the compounding frequency
#'   for nominal rates.
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{family}{Rate family: \code{"effective"}, \code{"discount"},
#'     \code{"force"}, \code{"nominal_interest"}, or
#'     \code{"nominal_discount"}.}
#'   \item{notation}{Actuarial notation for the equivalent rate.}
#'   \item{m}{Compounding frequency used for nominal rates.}
#'   \item{description}{Human-readable description.}
#'   \item{value}{Equivalent rate value.}
#' }
#'
#' @details
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
#' interest_equivalents("nominal_interest", rate = 0.18, m = 4)
#' interest_equivalents("nominal_discount", rate = 0.10, m = 12)
#' interest_equivalents("force", rate = 0.12)
#' interest_equivalents("effective", rate = 0.08)
#'
#' # Batch use with purrr
#' if (requireNamespace("purrr", quietly = TRUE) &&
#'     requireNamespace("tibble", quietly = TRUE)) {
#'   cases <- tibble::tibble(
#'     type = c("effective", "force"),
#'     rate = c(0.08, 0.12),
#'     m = c(1, 1)
#'   )
#'
#'   purrr::pmap(cases, interest_equivalents)
#' }
#'
#' @export
interest_equivalents <- function(
    type = c("effective", "nominal_interest", "nominal_discount", "force"),
    rate,
    m = 1L
) {
  type <- match.arg(type)

  if (missing(rate)) {
    stop("`rate` must be provided.", call. = FALSE)
  }

  if (!is.numeric(rate) || length(rate) != 1L || is.na(rate) || !is.finite(rate)) {
    stop("`rate` must be a single finite numeric value.", call. = FALSE)
  }

  if (!is.numeric(m) || length(m) != 1L || is.na(m) || !is.finite(m) ||
      m <= 0 || m != floor(m)) {
    stop("`m` must be a single positive integer.", call. = FALSE)
  }

  m <- as.integer(m)

  i <- standardize_interest(type = type, rate = rate, m = m)

  d     <- i / (1 + i)
  v     <- 1 / (1 + i)
  delta <- log1p(i)
  j_m   <- m * expm1(log1p(i) / m)
  d_m   <- m * (1 - exp(-log1p(i) / m))

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
    value = c(i, d, v, delta, j_m, d_m)
  )
}
