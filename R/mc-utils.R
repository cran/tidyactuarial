#' Internal helper: validate numeric scalar
#'
#' @param x Object to validate.
#' @param arg Character string with the argument name.
#' @param min Minimum allowed value. Default is `-Inf`.
#' @param max Maximum allowed value. Default is `Inf`.
#' @param strict_min Logical. Should the lower bound be strict?
#' @param strict_max Logical. Should the upper bound be strict?
#' @param allow_null Logical. Should `NULL` be allowed?
#'
#' @return Invisibly returns `TRUE` if validation is successful.
#'
#' @keywords internal
.mc_assert_numeric_scalar <- function(x,
                                      arg,
                                      min = -Inf,
                                      max = Inf,
                                      strict_min = FALSE,
                                      strict_max = FALSE,
                                      allow_null = FALSE) {
  if (is.null(x) && allow_null) {
    return(invisible(TRUE))
  }

  if (is.null(x) && !allow_null) {
    stop("`", arg, "` must not be NULL.", call. = FALSE)
  }

  if (!is.numeric(x) || length(x) != 1 || is.na(x) || !is.finite(x)) {
    stop("`", arg, "` must be a finite numeric scalar.", call. = FALSE)
  }

  lower_ok <- if (strict_min) x > min else x >= min
  upper_ok <- if (strict_max) x < max else x <= max

  if (!lower_ok || !upper_ok) {
    stop(
      "`", arg, "` is outside the allowed range.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Internal helper: validate positive integer-like scalar
#'
#' @param x Object to validate.
#' @param arg Character string with the argument name.
#'
#' @return Invisibly returns `TRUE` if validation is successful.
#'
#' @keywords internal
.mc_assert_positive_integer <- function(x, arg) {
  .mc_assert_numeric_scalar(
    x = x,
    arg = arg,
    min = 0,
    strict_min = TRUE
  )

  if (abs(x - round(x)) > sqrt(.Machine$double.eps)) {
    stop("`", arg, "` must be a positive integer.", call. = FALSE)
  }

  invisible(TRUE)
}


#' Internal helper: validate character scalar
#'
#' @param x Object to validate.
#' @param arg Character string with the argument name.
#'
#' @return Invisibly returns `TRUE` if validation is successful.
#'
#' @keywords internal
.mc_assert_character_scalar <- function(x, arg) {
  if (!is.character(x) || length(x) != 1 || is.na(x)) {
    stop("`", arg, "` must be a character string.", call. = FALSE)
  }

  invisible(TRUE)
}


#' Internal helper: validate column existence
#'
#' @param data A data frame or tibble.
#' @param col Character string with the column name.
#' @param arg Character string with the argument name.
#'
#' @return Invisibly returns `TRUE` if validation is successful.
#'
#' @keywords internal
.mc_assert_column <- function(data, col, arg) {
  .mc_assert_character_scalar(col, arg)

  if (!col %in% names(data)) {
    stop(
      "`", arg, "` must identify a column in `data`.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Internal helper: validate numeric column
#'
#' @param data A data frame or tibble.
#' @param col Character string with the column name.
#' @param arg Character string with the argument name.
#'
#' @return Invisibly returns `TRUE` if validation is successful.
#'
#' @keywords internal
.mc_assert_numeric_column <- function(data, col, arg) {
  .mc_assert_column(data = data, col = col, arg = arg)

  if (!is.numeric(data[[col]])) {
    stop(
      "`", arg, "` must identify a numeric column in `data`.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Internal helper: convert interest rate to annual effective rate
#'
#' Converts an interest rate supplied under a standard actuarial convention
#' into an equivalent annual effective interest rate.
#'
#' @param rate Numeric scalar. Interest rate.
#' @param interest_type Character string. Interest rate convention. One of
#'   `"effective"`, `"nominal"`, or `"force"`.
#' @param m Numeric scalar. Number of interest conversion periods per year when
#'   `interest_type = "nominal"`. Default is `1`.
#'
#' @details
#' This helper separates the financial interest-rate convention from the
#' actuarial cash-flow structure. The argument `m` is used only to convert a
#' nominal annual rate into an equivalent annual effective rate. It does not
#' represent annuity payment frequency.
#'
#' The conversion rules are:
#'
#' * Effective annual rate:
#'
#' \deqn{
#'   i_{\mathrm{eff}} = i.
#' }
#'
#' * Nominal annual rate convertible `m` times per year:
#'
#' \deqn{
#'   i_{\mathrm{eff}} = \left(1 + \frac{i^{(m)}}{m}\right)^m - 1.
#' }
#'
#' * Constant force of interest:
#'
#' \deqn{
#'   i_{\mathrm{eff}} = e^\delta - 1.
#' }
#'
#' @return Numeric scalar with the equivalent annual effective interest rate.
#'
#' @keywords internal
.mc_effective_rate <- function(rate,
                               interest_type = c("effective", "nominal", "force"),
                               m = 1) {
  interest_type <- match.arg(interest_type)

  .mc_assert_numeric_scalar(rate, "rate")
  .mc_assert_numeric_scalar(m, "m", min = 0, strict_min = TRUE)

  if (interest_type == "effective" && rate <= -1) {
    stop(
      "`rate` must be greater than -1 when `interest_type = 'effective'`.",
      call. = FALSE
    )
  }

  if (interest_type == "nominal" && (1 + rate / m) <= 0) {
    stop(
      "`1 + rate / m` must be positive when `interest_type = 'nominal'`.",
      call. = FALSE
    )
  }

  switch(
    interest_type,
    effective = rate,
    nominal = (1 + rate / m)^m - 1,
    force = exp(rate) - 1
  )
}


#' Internal helper: compute annual discount factor
#'
#' Computes the annual discount factor from an interest rate convention.
#'
#' @param rate Numeric scalar. Interest rate.
#' @param interest_type Character string. Interest rate convention. One of
#'   `"effective"`, `"nominal"`, or `"force"`.
#' @param m Numeric scalar. Number of interest conversion periods per year when
#'   `interest_type = "nominal"`. Default is `1`.
#'
#' @return Numeric scalar with the annual discount factor.
#'
#' @keywords internal
.mc_discount_factor <- function(rate,
                                interest_type = c("effective", "nominal", "force"),
                                m = 1) {
  interest_type <- match.arg(interest_type)

  effective_rate <- .mc_effective_rate(
    rate = rate,
    interest_type = interest_type,
    m = m
  )

  1 / (1 + effective_rate)
}


#' Internal helper: safe payment times
#'
#' Creates a sequence of payment times. If the upper bound is smaller than the
#' lower bound, it returns `numeric(0)`.
#'
#' @param from Numeric scalar. First payment time.
#' @param to Numeric scalar. Last payment time.
#' @param by Numeric scalar. Step between payment times. Default is `1`.
#'
#' @return Numeric vector.
#'
#' @keywords internal
.mc_payment_times <- function(from, to, by = 1) {
  if (is.null(from) || is.null(to)) {
    return(numeric(0))
  }

  if (length(from) != 1 || length(to) != 1 || length(by) != 1) {
    stop("`from`, `to`, and `by` must be numeric scalars.", call. = FALSE)
  }

  if (is.na(from) || is.na(to) || is.na(by)) {
    return(numeric(0))
  }

  if (!is.finite(from) || !is.finite(to) || !is.finite(by)) {
    return(numeric(0))
  }

  if (by <= 0) {
    stop("`by` must be positive.", call. = FALSE)
  }

  if (to < from) {
    return(numeric(0))
  }

  seq(from = from, to = to, by = by)
}


#' Internal helper: create quantile column names
#'
#' @param probs Numeric vector of probabilities.
#'
#' @return Character vector with syntactically valid quantile names.
#'
#' @keywords internal
.mc_quantile_names <- function(probs) {
  raw <- 100 * probs

  names <- vapply(
    raw,
    function(x) {
      label <- formatC(x, format = "f", digits = 3)
      label <- sub("0+$", "", label)
      label <- sub("\\.$", "", label)
      label <- gsub("\\.", "_", label)
      paste0("q", label)
    },
    character(1)
  )

  names
}
