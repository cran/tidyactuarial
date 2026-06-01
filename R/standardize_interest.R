#' Standardize an interest rate to the annual effective rate i
#'
#' Converts common interest-rate specifications to the equivalent annual
#' effective interest rate \code{i}, using compact actuarial notation.
#'
#' @param i_type Character vector indicating the interest-rate type. Must be one
#'   of \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, or \code{"force"}.
#' @param i Numeric vector of interest-rate values.
#' @param m Positive integer vector giving the conversion frequency for nominal
#'   rates. Ignored for \code{"effective"} and \code{"force"}.
#' @param ... Transitional compatibility for older internal calls using
#'   \code{type = } and \code{rate = }. These names are accepted and mapped to
#'   \code{i_type} and \code{i}.
#'
#' @return Numeric vector of annual effective rates. Missing values are
#' propagated.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{i} denotes the interest-rate value,
#' \code{i_type} denotes the interest-rate type, and \code{m} denotes the
#' conversion frequency for nominal rates.
#'
#' The conversion formulas are:
#' \describe{
#'   \item{effective:}{\eqn{i = i}{i = i} (identity)}
#'   \item{nominal_interest:}{\eqn{i_e = (1 + j^{(m)}/m)^m - 1}{i_e = (1 + j(m)/m)^m - 1}}
#'   \item{nominal_discount:}{\eqn{i_e = (1 - d^{(m)}/m)^{-m} - 1}{i_e = (1 - d(m)/m)^(-m) - 1}}
#'   \item{force:}{\eqn{i_e = e^{\delta} - 1}{i_e = exp(delta) - 1}}
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
#' standardize_interest(i_type = "nominal_interest", i = 0.18, m = 4)
#' standardize_interest(i_type = "nominal_discount", i = 0.10, m = 12)
#' standardize_interest(i_type = "force", i = 0.12)
#'
#' # Vectorized case
#' standardize_interest(
#'   i_type = c("nominal_interest", "force", "effective"),
#'   i = c(0.06, 0.05, 0.04),
#'   m = c(12, 1, 1)
#' )
#'
#' # Use inside a data pipeline
#' if (requireNamespace("dplyr", quietly = TRUE) &&
#'     requireNamespace("tibble", quietly = TRUE)) {
#'   portfolio <- tibble::tibble(
#'     policy_id = 1:3,
#'     i = c(0.05, 0.08, 0.10),
#'     i_type = c("force", "nominal_interest", "nominal_discount"),
#'     m = c(1, 4, 12)
#'   )
#'
#'   dplyr::mutate(
#'     portfolio,
#'     i_effective = standardize_interest(
#'       i_type = i_type,
#'       i = i,
#'       m = m
#'     )
#'   )
#' }
#'
#' @export
standardize_interest <- function(
    i_type = "effective",
    i,
    m = 1,
    ...
) {
  dots <- list(...)

  if (length(dots) > 0L) {
    allowed_old <- c("type", "rate")
    bad_dots <- setdiff(names(dots), allowed_old)

    if (length(bad_dots) > 0L) {
      stop(
        "Unused argument(s): ",
        paste(sprintf("`%s`", bad_dots), collapse = ", "),
        ".",
        call. = FALSE
      )
    }

    # Transitional compatibility: many internal calls still use
    # standardize_interest(type = ..., rate = ..., m = ...).
    if (!is.null(dots$type)) {
      i_type <- dots$type
    }

    if (missing(i) && !is.null(dots$rate)) {
      i <- dots$rate
    } else if (!missing(i) && !is.null(dots$rate)) {
      stop("Provide only one of `i` or deprecated `rate`.", call. = FALSE)
    }
  }

  if (missing(i)) {
    stop("`i` must be provided.", call. = FALSE)
  }

  if (!is.character(i_type)) {
    stop("`i_type` must be a character vector.", call. = FALSE)
  }

  if (!is.numeric(i)) {
    stop("`i` must be a numeric vector.", call. = FALSE)
  }

  if (!is.numeric(m)) {
    stop("`m` must be a numeric vector.", call. = FALSE)
  }

  size <- max(length(i_type), length(i), length(m))

  valid_size <- function(x) length(x) %in% c(1L, size)

  if (!valid_size(i_type) || !valid_size(i) || !valid_size(m)) {
    stop(
      "`i_type`, `i`, and `m` must have length 1 or a common length.",
      call. = FALSE
    )
  }

  i_type <- rep_len(tolower(i_type), size)
  i      <- rep_len(i, size)
  m      <- rep_len(m, size)

  valid_types <- c(
    "effective",
    "nominal_interest",
    "nominal_discount",
    "force"
  )

  bad_type <- !is.na(i_type) & !(i_type %in% valid_types)
  if (any(bad_type)) {
    stop(
      "`i_type` must contain only 'effective', 'nominal_interest', ",
      "'nominal_discount', or 'force'.",
      call. = FALSE
    )
  }

  bad_i <- !is.na(i) & !is.finite(i)
  if (any(bad_i)) {
    stop("`i` must contain only finite numeric values or NA.", call. = FALSE)
  }

  idx_nom_i <- i_type == "nominal_interest"
  idx_nom_d <- i_type == "nominal_discount"
  idx_eff   <- i_type == "effective"
  idx_force <- i_type == "force"
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

  bad_eff <- idx_eff & !is.na(i) & i <= -1
  if (any(bad_eff)) {
    stop(
      "For `i_type = 'effective'`, `i` must be greater than -1.",
      call. = FALSE
    )
  }

  bad_nom_i <- idx_nom_i & !is.na(i) & (1 + i / m) <= 0
  if (any(bad_nom_i)) {
    stop(
      "For `i_type = 'nominal_interest'`, `1 + i / m` must be positive.",
      call. = FALSE
    )
  }

  bad_nom_d <- idx_nom_d & !is.na(i) & (1 - i / m) <= 0
  if (any(bad_nom_d)) {
    stop(
      "For `i_type = 'nominal_discount'`, `1 - i / m` must be positive.",
      call. = FALSE
    )
  }

  out <- rep(NA_real_, size)

  ok_eff <- idx_eff & !is.na(i)
  out[ok_eff] <- i[ok_eff]

  ok_nom_i <- idx_nom_i & !is.na(i)
  out[ok_nom_i] <- (1 + i[ok_nom_i] / m[ok_nom_i])^m[ok_nom_i] - 1

  ok_nom_d <- idx_nom_d & !is.na(i)
  out[ok_nom_d] <- (1 - i[ok_nom_d] / m[ok_nom_d])^(-m[ok_nom_d]) - 1

  ok_force <- idx_force & !is.na(i)
  out[ok_force] <- exp(i[ok_force]) - 1

  out
}
