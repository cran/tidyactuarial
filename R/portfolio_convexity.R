#' Compute portfolio convexity as a market-value-weighted average
#'
#' Computes portfolio convexity from individual position convexities using
#' market values as weights.
#'
#' This is a summarise-style tibble-first function. Each input row represents
#' one position, and each output row represents one portfolio.
#'
#' The function does not compute individual convexities from bond terms
#' or yields. Instead, it assumes that the input convexity column already
#' contains valid convexity measures on a common basis within each portfolio.
#'
#' The portfolio convexity is computed as:
#' \deqn{C_P = \frac{\sum_{k=1}^n P_k C_k}{\sum_{k=1}^n P_k}}{C_P = sum(P_k * C_k) / sum(P_k)}
#' where \eqn{P_k}{P_k} is the market value of position \eqn{k} and
#' \eqn{C_k}{C_k} is its convexity.
#'
#' @param .data A data.frame or tibble. If \code{NULL}, \code{market_value} and
#'   \code{convexity} must be supplied as vectors.
#' @param portfolio_id Optional vector of portfolio identifiers when
#'   \code{.data = NULL}. If omitted, all positions are treated as belonging to
#'   a single portfolio.
#' @param market_value Numeric vector of market values when \code{.data = NULL}.
#' @param convexity Numeric vector of individual convexities when
#'   \code{.data = NULL}.
#' @param col_portfolio Name of the portfolio identifier column. If \code{NULL},
#'   all rows are treated as one portfolio.
#' @param col_market_value Name of the numeric column containing market values.
#' @param col_convexity Name of the numeric column containing individual
#'   convexities.
#' @param .out Name of the output column containing portfolio convexity.
#' @param .out_value Name of the output column containing total portfolio
#'   market value.
#' @param .out_n Name of the output column containing the number of positions
#'   used in the calculation.
#' @param .na NA handling policy: \code{"propagate"}, \code{"error"}, or
#'   \code{"drop"}.
#'
#' @return A tibble with one row per portfolio and columns for portfolio
#'   convexity, total market value, and number of positions used.
#'
#' @seealso \code{\link{portfolio_duration}},
#'   \code{\link{bond_convexity}}, \code{\link{bond_duration}}
#'
#' @family bonds
#'
#' @examples
#' # Simple example: one portfolio
#' portfolio_convexity(
#'   market_value = c(1000, 2000, 500),
#'   convexity = c(20, 12, 35)
#' )
#'
#' # Medium example: two portfolios
#' positions <- tibble::tibble(
#'   portfolio_id = c("A", "A", "B", "B"),
#'   market_value = c(1000, 2000, 1000 / 1.08^2, 1000 / 1.08^4),
#'   convexity = c(20, 12, 6, 18)
#' )
#'
#' portfolio_convexity(
#'   positions,
#'   col_portfolio = "portfolio_id",
#'   col_market_value = "market_value",
#'   col_convexity = "convexity"
#' )
#'
#' @references
#' Marcel B. Finan, *A Basic Course in the Theory of Interest and
#' Derivatives Markets: A Preparation for the Actuarial Exam FM/2*,
#' Section 55: Redington Immunization and Convexity.
#'
#' Kellison, S. G. *The Theory of Interest*, Chapter 11:
#' Duration, Convexity and Immunization.
#'
#' @export
portfolio_convexity <- function(
    .data = NULL,
    portfolio_id = NULL,
    market_value = NULL,
    convexity = NULL,
    col_portfolio = "portfolio_id",
    col_market_value = "market_value",
    col_convexity = "convexity",
    .out = "portfolio_convexity",
    .out_value = "portfolio_market_value",
    .out_n = "n_positions",
    .na = c("propagate", "error", "drop")
) {
  .na <- match.arg(.na)

  abort <- function(msg) rlang::abort(msg)

  if (!is.null(col_portfolio) &&
      (!is.character(col_portfolio) || length(col_portfolio) != 1L ||
       is.na(col_portfolio) || !nzchar(col_portfolio))) {
    abort("`col_portfolio` must be NULL or a single non-empty string.")
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

  if (is.null(.data)) {
    if (is.null(market_value) || is.null(convexity)) {
      abort("When `.data = NULL`, `market_value` and `convexity` must be supplied.")
    }

    if (!is.numeric(market_value)) {
      abort("`market_value` must be a numeric vector when `.data = NULL`.")
    }

    if (!is.numeric(convexity)) {
      abort("`convexity` must be a numeric vector when `.data = NULL`.")
    }

    n <- length(market_value)

    if (length(convexity) != n) {
      abort("`market_value` and `convexity` must have the same length.")
    }

    if (is.null(portfolio_id)) {
      portfolio_id <- rep(1L, n)
      has_portfolio_output <- FALSE
    } else {
      if (length(portfolio_id) != n) {
        abort("`portfolio_id` must have the same length as `market_value` and `convexity`.")
      }

      has_portfolio_output <- TRUE
    }

    data_in <- tibble::tibble(
      .portfolio = portfolio_id,
      .market_value = market_value,
      .convexity = convexity,
      .row_id = seq_len(n)
    )
  } else {
    if (!inherits(.data, "data.frame")) {
      abort("`.data` must be a data.frame or tibble.")
    }

    data_in <- tibble::as_tibble(.data)
    data_in$.row_id <- seq_len(nrow(data_in))

    if (!col_market_value %in% names(data_in)) {
      abort(paste0(
        "Column `", col_market_value, "` was not found. ",
        "Pass the correct column name via `col_market_value`."
      ))
    }

    if (!col_convexity %in% names(data_in)) {
      abort(paste0(
        "Column `", col_convexity, "` was not found. ",
        "Pass the correct column name via `col_convexity`."
      ))
    }

    if (is.null(col_portfolio)) {
      data_in$.portfolio <- 1L
      has_portfolio_output <- FALSE
    } else {
      if (!col_portfolio %in% names(data_in)) {
        abort(paste0(
          "Column `", col_portfolio, "` was not found. ",
          "Pass the correct column name via `col_portfolio`, or set `col_portfolio = NULL` ",
          "to treat all rows as one portfolio."
        ))
      }

      data_in$.portfolio <- data_in[[col_portfolio]]
      has_portfolio_output <- TRUE
    }

    data_in$.market_value <- data_in[[col_market_value]]
    data_in$.convexity <- data_in[[col_convexity]]

    if (!is.numeric(data_in$.market_value)) {
      abort(paste0("`", col_market_value, "` must be numeric."))
    }

    if (!is.numeric(data_in$.convexity)) {
      abort(paste0("`", col_convexity, "` must be numeric."))
    }
  }

  data_in$.bad_na <- is.na(data_in$.market_value) |
    is.na(data_in$.convexity) |
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
    bad_mv <- !is.finite(data_valid$.market_value) |
      data_valid$.market_value < 0

    if (any(bad_mv)) {
      bad_rows <- data_valid$.row_id[bad_mv]

      abort(paste0(
        "`", if (is.null(.data)) "market_value" else col_market_value,
        "` must be finite and >= 0. Problem at row(s): ",
        paste(bad_rows, collapse = ", "), "."
      ))
    }

    bad_cx <- !is.finite(data_valid$.convexity) |
      data_valid$.convexity < 0

    if (any(bad_cx)) {
      bad_rows <- data_valid$.row_id[bad_cx]

      abort(paste0(
        "`", if (is.null(.data)) "convexity" else col_convexity,
        "` must be finite and >= 0. Problem at row(s): ",
        paste(bad_rows, collapse = ", "), "."
      ))
    }
  }

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

  split_all <- split(data_work, data_work$.portfolio, drop = TRUE)

  res_list <- lapply(split_all, function(df_port) {
    has_na_port <- any(df_port$.bad_na)

    if (.na == "propagate" && has_na_port) {
      tibble::tibble(
        .portfolio = df_port$.portfolio[[1]],
        .portfolio_convexity = NA_real_,
        .portfolio_market_value = NA_real_,
        .n_positions = NA_integer_
      )
    } else {
      df_ok <- df_port[!df_port$.bad_na, , drop = FALSE]
      total_mv <- sum(df_ok$.market_value)

      if (!is.finite(total_mv) || total_mv <= 0) {
        abort(paste0(
          "Each portfolio must have strictly positive total market value. ",
          "Problem detected for portfolio `",
          as.character(df_port$.portfolio[[1]]),
          "`."
        ))
      }

      tibble::tibble(
        .portfolio = df_port$.portfolio[[1]],
        .portfolio_convexity = sum(df_ok$.market_value * df_ok$.convexity) / total_mv,
        .portfolio_market_value = total_mv,
        .n_positions = nrow(df_ok)
      )
    }
  })

  out <- dplyr::bind_rows(res_list)

  names(out)[names(out) == ".portfolio_convexity"] <- .out
  names(out)[names(out) == ".portfolio_market_value"] <- .out_value
  names(out)[names(out) == ".n_positions"] <- .out_n

  if (has_portfolio_output) {
    portfolio_name_out <- if (is.null(.data)) "portfolio_id" else col_portfolio
    names(out)[names(out) == ".portfolio"] <- portfolio_name_out
  } else {
    out$.portfolio <- NULL
  }

  tibble::as_tibble(out)
}
