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
  if (!is.character(x) ||
      length(x) != 1L ||
      is.na(x) ||
      !nzchar(trimws(x))) {
    stop(
      "`", arg, "` must be a single nonempty character string.",
      call. = FALSE
    )
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
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame or tibble.", call. = FALSE)
  }

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


#' Internal helper: normalize Monte Carlo interest-rate type
#'
#' Converts legacy Monte Carlo labels to the package-wide actuarial labels.
#'
#' @param i_type Character scalar. Interest-rate type.
#'
#' @return Character scalar.
#'
#' @keywords internal
.mc_normalize_i_type <- function(i_type) {
  .mc_assert_character_scalar(i_type, "i_type")

  i_type <- tolower(i_type)

  if (identical(i_type, "nominal")) {
    i_type <- "nominal_interest"
  }

  valid_i_type <- c(
    "effective",
    "nominal_interest",
    "nominal_discount",
    "force"
  )

  if (!i_type %in% valid_i_type) {
    stop(
      "`i_type` must be one of: ",
      paste(sprintf("'%s'", valid_i_type), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  i_type
}


#' Internal helper: convert interest rate to annual effective rate
#'
#' Converts an interest rate supplied under a standard actuarial convention into
#' an equivalent annual effective interest rate.
#'
#' @param i Numeric scalar. Interest-rate input.
#' @param i_type Character string. Interest-rate convention. One of
#'   `"effective"`, `"nominal_interest"`, `"nominal_discount"`, or `"force"`.
#'   The legacy value `"nominal"` is accepted internally and treated as
#'   `"nominal_interest"`.
#' @param m Numeric scalar. Number of interest conversion periods per year for
#'   nominal annual rates. Default is `1`.
#' @param rate Deprecated internal alias for `i`.
#' @param interest_type Deprecated internal alias for `i_type`.
#'
#' @details
#' This helper follows the compact actuarial notation used throughout
#' `tidyactuarial`: `i` is the interest-rate input, `i_type` is the
#' interest-rate type, and `m` is the conversion frequency for nominal rates.
#'
#' The argument `m` is used only to convert nominal annual rates into equivalent
#' annual effective rates. It does not represent annuity payment frequency.
#'
#' Transitional compatibility is intentionally retained for older internal Monte
#' Carlo calls using `rate =` and `interest_type =`.
#'
#' @return Numeric scalar with the equivalent annual effective interest rate.
#'
#' @keywords internal
.mc_effective_rate <- function(i,
                               i_type = c(
                                 "effective",
                                 "nominal_interest",
                                 "nominal_discount",
                                 "force"
                               ),
                               m = 1,
                               rate = NULL,
                               interest_type = NULL) {
  i_missing <- missing(i)
  i_type_missing <- missing(i_type)

  if (!is.null(rate)) {
    if (!i_missing) {
      stop(
        "Provide only one of `i` or deprecated `rate`.",
        call. = FALSE
      )
    }

    i <- rate
    i_missing <- FALSE
  }

  if (!is.null(interest_type)) {
    if (!i_type_missing) {
      stop(
        "Provide only one of `i_type` or deprecated `interest_type`.",
        call. = FALSE
      )
    }

    i_type <- interest_type
    i_type_missing <- FALSE
  }

  if (i_missing) {
    stop("`i` must be provided.", call. = FALSE)
  }

  if (i_type_missing) {
    i_type <- "effective"
  }

  i_type <- .mc_normalize_i_type(i_type)

  .mc_assert_numeric_scalar(i, "i")
  .mc_assert_positive_integer(m, "m")

  m <- as.integer(round(m))

  i_effective <- standardize_interest(
    i_type = i_type,
    i = i,
    m = m
  )

  if (!is.numeric(i_effective) ||
      length(i_effective) != 1L ||
      is.na(i_effective) ||
      !is.finite(i_effective) ||
      i_effective <= -1) {
    stop(
      "The supplied interest convention must imply a finite annual ",
      "effective rate greater than -1.",
      call. = FALSE
    )
  }

  i_effective
}


#' Internal helper: compute annual discount factor
#'
#' Computes the annual discount factor from an interest-rate convention.
#'
#' @param i Numeric scalar. Interest-rate input.
#' @param i_type Character string. Interest-rate convention. One of
#'   `"effective"`, `"nominal_interest"`, `"nominal_discount"`, or `"force"`.
#'   The legacy value `"nominal"` is accepted internally and treated as
#'   `"nominal_interest"`.
#' @param m Numeric scalar. Number of interest conversion periods per year for
#'   nominal annual rates. Default is `1`.
#' @param rate Deprecated internal alias for `i`.
#' @param interest_type Deprecated internal alias for `i_type`.
#'
#' @return Numeric scalar with the annual discount factor.
#'
#' @keywords internal
.mc_discount_factor <- function(i,
                                i_type = c(
                                  "effective",
                                  "nominal_interest",
                                  "nominal_discount",
                                  "force"
                                ),
                                m = 1,
                                rate = NULL,
                                interest_type = NULL) {
  args <- list(m = m)

  if (!missing(i)) {
    args$i <- i
  }

  if (!missing(i_type)) {
    args$i_type <- i_type
  }

  if (!is.null(rate)) {
    args$rate <- rate
  }

  if (!is.null(interest_type)) {
    args$interest_type <- interest_type
  }

  effective_rate <- do.call(
    .mc_effective_rate,
    args
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

  if (!is.numeric(from) ||
      !is.numeric(to) ||
      !is.numeric(by) ||
      length(from) != 1L ||
      length(to) != 1L ||
      length(by) != 1L) {
    stop(
      "`from`, `to`, and `by` must be numeric scalars.",
      call. = FALSE
    )
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

  tol <- sqrt(.Machine$double.eps) *
    max(1, abs(from), abs(to), abs(by))

  if (to < from - tol) {
    return(numeric(0))
  }

  n_steps <- floor((to - from) / by + tol)

  if (!is.finite(n_steps) || n_steps < 0) {
    return(numeric(0))
  }

  times <- from + seq.int(0, n_steps) * by
  times[times <= to + tol]
}


#' Internal helper: create quantile column names
#'
#' @param probs Numeric vector of probabilities.
#'
#' @return Character vector with syntactically valid quantile names.
#'
#' @keywords internal
.mc_quantile_names <- function(probs) {
  if (!is.numeric(probs) ||
      length(probs) == 0L ||
      anyNA(probs) ||
      any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1)) {
    stop(
      "`probs` must be a nonempty numeric vector with values in [0, 1].",
      call. = FALSE
    )
  }

  raw <- 100 * probs

  quantile_names <- vapply(
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

  if (anyDuplicated(quantile_names)) {
    stop(
      "`probs` produces duplicated quantile column names at the supported ",
      "display precision.",
      call. = FALSE
    )
  }

  quantile_names
}

