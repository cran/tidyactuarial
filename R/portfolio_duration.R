#' Compute portfolio duration as a market-value-weighted average
#'
#' Computes portfolio duration from individual position durations using present
#' values, prices, or market values as weights, using compact actuarial
#' notation.
#'
#' This is a summarise-style tibble-first function. Each input row represents
#' one position, and each output row represents one portfolio.
#'
#' The function does not compute individual durations from bond terms or yields.
#' Instead, it assumes that the input duration column already contains valid
#' duration measures on a common basis within each portfolio.
#'
#' The portfolio duration is computed as:
#' \deqn{D_P = \frac{\sum_{j=1}^r P_j D_j}{\sum_{j=1}^r P_j}}
#' where \eqn{P_j} is the present value, price, or market value of position
#' \eqn{j}, and \eqn{D_j} is its duration.
#'
#' @param .data A data.frame or tibble. If \code{NULL}, \code{P} and
#'   \code{D} must be supplied as vectors.
#' @param portfolio_id Optional vector of portfolio identifiers when
#'   \code{.data = NULL}. If omitted, all positions are treated as belonging to
#'   a single portfolio.
#' @param P Numeric vector of present values, prices, or market values when
#'   \code{.data = NULL}.
#' @param D Numeric vector of individual durations when \code{.data = NULL}.
#' @param col_portfolio Name of the portfolio identifier column. If \code{NULL},
#'   all rows are treated as one portfolio.
#' @param col_P Name of the numeric column containing present values, prices, or
#'   market values.
#' @param col_D Name of the numeric column containing individual durations.
#' @param .out Name of the output column containing portfolio duration.
#' @param .out_value Name of the output column containing total portfolio value.
#' @param .out_n Name of the output column containing the number of positions
#'   used in the calculation.
#' @param .na NA handling policy: \code{"propagate"}, \code{"error"}, or
#'   \code{"drop"}.
#' @param ... Transitional compatibility for older calls using
#'   \code{market_value}, \code{duration}, \code{col_market_value}, and
#'   \code{col_duration}. These names are mapped to \code{P}, \code{D},
#'   \code{col_P}, and \code{col_D}.
#'
#' @return A tibble with one row per portfolio and columns for portfolio
#'   duration, total portfolio value, and number of positions used.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{P} denotes price, present value, or market value,
#' and \code{D} denotes duration.
#'
#' The function is deliberately agnostic about the duration convention, but all
#' individual durations must be expressed on the same basis within each
#' portfolio. For example, do not mix Macaulay durations in years with
#' durations measured in coupon periods.
#'
#' @seealso \code{\link{portfolio_convexity}},
#'   \code{\link{bond_duration}}, \code{\link{bond_convexity}}
#'
#' @family bonds
#'
#' @examples
#' # Simple example: one portfolio
#' portfolio_duration(
#'   P = c(1000, 2000, 500),
#'   D = c(7, 5, 10)
#' )
#'
#' # Medium example: two portfolios
#' positions <- tibble::tibble(
#'   portfolio_id = c("A", "A", "B", "B"),
#'   P = c(1000, 2000, 1000 / 1.08^2, 1000 / 1.08^4),
#'   D = c(7, 5, 2, 4)
#' )
#'
#' portfolio_duration(
#'   positions,
#'   col_portfolio = "portfolio_id",
#'   col_P = "P",
#'   col_D = "D"
#' )
#'
#' @references
#' Marcel B. Finan, \emph{A Basic Course in the Theory of Interest and
#' Derivatives Markets: A Preparation for the Actuarial Exam FM/2},
#' Section 54: Macaulay and Modified Durations.
#'
#' Kellison, S. G. \emph{The Theory of Interest}, Chapter 11:
#' Duration, Convexity and Immunization.
#'
#' @export
portfolio_duration <- function(
    .data = NULL,
    portfolio_id = NULL,
    P = NULL,
    D = NULL,
    col_portfolio = "portfolio_id",
    col_P = "P",
    col_D = "D",
    .out = "D_P",
    .out_value = "P_total",
    .out_n = "n_positions",
    .na = c("propagate", "error", "drop"),
    ...
) {
  .na <- match.arg(.na)
  dots <- list(...)

  abort <- function(msg) rlang::abort(msg)

  # --- Transitional compatibility with old argument names ---
  allowed_old <- c(
    "market_value",
    "duration",
    "col_market_value",
    "col_duration"
  )
  bad_dots <- setdiff(names(dots), allowed_old)

  if (length(bad_dots) > 0L) {
    abort(paste0(
      "Unused argument(s): ",
      paste(sprintf("`%s`", bad_dots), collapse = ", "),
      "."
    ))
  }

  if (!is.null(dots$market_value)) {
    if (!is.null(P)) {
      abort("Provide only one of `P` or deprecated `market_value`.")
    }
    P <- dots$market_value
  }

  if (!is.null(dots$duration)) {
    if (!is.null(D)) {
      abort("Provide only one of `D` or deprecated `duration`.")
    }
    D <- dots$duration
  }

  if (!is.null(dots$col_market_value)) {
    if (!identical(col_P, "P")) {
      abort("Provide only one of `col_P` or deprecated `col_market_value`.")
    }
    col_P <- dots$col_market_value
  }

  if (!is.null(dots$col_duration)) {
    if (!identical(col_D, "D")) {
      abort("Provide only one of `col_D` or deprecated `col_duration`.")
    }
    col_D <- dots$col_duration
  }

  # --- Validate names ---
  if (!is.null(col_portfolio) &&
      (!is.character(col_portfolio) || length(col_portfolio) != 1L ||
       is.na(col_portfolio) || !nzchar(col_portfolio))) {
    abort("`col_portfolio` must be NULL or a single non-empty string.")
  }

  if (!is.character(col_P) || length(col_P) != 1L ||
      is.na(col_P) || !nzchar(col_P)) {
    abort("`col_P` must be a single non-empty string.")
  }

  if (!is.character(col_D) || length(col_D) != 1L ||
      is.na(col_D) || !nzchar(col_D)) {
    abort("`col_D` must be a single non-empty string.")
  }

  for (nm in c(".out", ".out_value", ".out_n")) {
    val <- get(nm, inherits = FALSE)

    if (!is.character(val) || length(val) != 1L || is.na(val) || !nzchar(val)) {
      abort(paste0("`", nm, "` must be a single non-empty string."))
    }
  }

  out_names <- c(.out, .out_value, .out_n)

  if (length(unique(out_names)) != 3L) {
    abort("`.out`, `.out_value`, and `.out_n` must have different names.")
  }

  # --- Build internal working data ---
  if (is.null(.data)) {
    if (is.null(P) || is.null(D)) {
      abort("When `.data = NULL`, `P` and `D` must be supplied.")
    }

    if (!is.numeric(P)) {
      abort("`P` must be a numeric vector when `.data = NULL`.")
    }

    if (!is.numeric(D)) {
      abort("`D` must be a numeric vector when `.data = NULL`.")
    }

    n <- length(P)

    if (length(D) != n) {
      abort("`P` and `D` must have the same length.")
    }

    if (is.null(portfolio_id)) {
      portfolio_id <- rep(1L, n)
      has_portfolio_output <- FALSE
    } else {
      if (length(portfolio_id) != n) {
        abort("`portfolio_id` must have the same length as `P` and `D`.")
      }

      has_portfolio_output <- TRUE
    }

    data_in <- tibble::tibble(
      .portfolio = portfolio_id,
      .P = P,
      .D = D,
      .row_id = seq_len(n)
    )
  } else {
    if (!inherits(.data, "data.frame")) {
      abort("`.data` must be a data.frame or tibble.")
    }

    data_in <- tibble::as_tibble(.data)
    data_in$.row_id <- seq_len(nrow(data_in))

    if (!col_P %in% names(data_in)) {
      abort(paste0(
        "Column `", col_P, "` was not found. ",
        "Pass the correct column name via `col_P`."
      ))
    }

    if (!col_D %in% names(data_in)) {
      abort(paste0(
        "Column `", col_D, "` was not found. ",
        "Pass the correct column name via `col_D`."
      ))
    }

    if (is.null(col_portfolio)) {
      data_in$.portfolio <- 1L
      has_portfolio_output <- FALSE
    } else {
      if (!col_portfolio %in% names(data_in)) {
        abort(paste0(
          "Column `", col_portfolio, "` was not found. ",
          "Pass the correct column name via `col_portfolio`, or set ",
          "`col_portfolio = NULL` to treat all rows as one portfolio."
        ))
      }

      data_in$.portfolio <- data_in[[col_portfolio]]
      has_portfolio_output <- TRUE
    }

    data_in$.P <- data_in[[col_P]]
    data_in$.D <- data_in[[col_D]]

    if (!is.numeric(data_in$.P)) {
      abort(paste0("`", col_P, "` must be numeric."))
    }

    if (!is.numeric(data_in$.D)) {
      abort(paste0("`", col_D, "` must be numeric."))
    }
  }

  # --- Missing-value handling ---
  data_in$.bad_na <- is.na(data_in$.P) |
    is.na(data_in$.D) |
    is.na(data_in$.portfolio)

  if (.na == "error" && any(data_in$.bad_na)) {
    bad_rows <- data_in$.row_id[data_in$.bad_na]

    abort(paste0(
      "Missing values found in required inputs at row(s): ",
      paste(bad_rows, collapse = ", "), "."
    ))
  }

  if (.na == "drop" && any(data_in$.bad_na)) {
    bad_rows <- data_in$.row_id[data_in$.bad_na]

    warning(
      paste0(
        "Dropping row(s) with missing values in required inputs: ",
        paste(bad_rows, collapse = ", "), "."
      ),
      call. = FALSE
    )

    data_work <- data_in[!data_in$.bad_na, , drop = FALSE]
  } else {
    data_work <- data_in
  }

  data_valid <- data_work[!data_work$.bad_na, , drop = FALSE]

  if (nrow(data_valid) > 0L) {
    bad_P <- !is.finite(data_valid$.P) |
      data_valid$.P < 0

    if (any(bad_P)) {
      bad_rows <- data_valid$.row_id[bad_P]

      abort(paste0(
        "`", if (is.null(.data)) "P" else col_P,
        "` must be finite and >= 0. Problem at row(s): ",
        paste(bad_rows, collapse = ", "), "."
      ))
    }

    bad_D <- !is.finite(data_valid$.D) |
      data_valid$.D < 0

    if (any(bad_D)) {
      bad_rows <- data_valid$.row_id[bad_D]

      abort(paste0(
        "`", if (is.null(.data)) "D" else col_D,
        "` must be finite and >= 0. Problem at row(s): ",
        paste(bad_rows, collapse = ", "), "."
      ))
    }
  }

  # --- Empty output ---
  if (nrow(data_work) == 0L) {
    out <- tibble::tibble()

    if (has_portfolio_output) {
      out[[if (is.null(.data)) "portfolio_id" else col_portfolio]] <-
        vector(mode = "logical", length = 0L)
    }

    out[[.out]] <- numeric(0)
    out[[.out_value]] <- numeric(0)
    out[[.out_n]] <- integer(0)

    return(out)
  }

  # --- Summarise by portfolio ---
  split_all <- split(data_work, data_work$.portfolio, drop = TRUE)

  res_list <- lapply(split_all, function(df_port) {
    has_na_port <- any(df_port$.bad_na)

    if (.na == "propagate" && has_na_port) {
      tibble::tibble(
        .portfolio = df_port$.portfolio[[1]],
        .D_P = NA_real_,
        .P_total = NA_real_,
        .n_positions = NA_integer_
      )
    } else {
      df_ok <- df_port[!df_port$.bad_na, , drop = FALSE]
      P_total <- sum(df_ok$.P)

      if (!is.finite(P_total) || P_total <= 0) {
        abort(paste0(
          "Each portfolio must have strictly positive total value. ",
          "Problem detected for portfolio `",
          as.character(df_port$.portfolio[[1]]),
          "`."
        ))
      }

      tibble::tibble(
        .portfolio = df_port$.portfolio[[1]],
        .D_P = sum(df_ok$.P * df_ok$.D) / P_total,
        .P_total = P_total,
        .n_positions = nrow(df_ok)
      )
    }
  })

  out <- dplyr::bind_rows(res_list)

  names(out)[names(out) == ".D_P"] <- .out
  names(out)[names(out) == ".P_total"] <- .out_value
  names(out)[names(out) == ".n_positions"] <- .out_n

  if (has_portfolio_output) {
    portfolio_name_out <- if (is.null(.data)) "portfolio_id" else col_portfolio
    names(out)[names(out) == ".portfolio"] <- portfolio_name_out
  } else {
    out$.portfolio <- NULL
  }

  tibble::as_tibble(out)
}
