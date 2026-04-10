#' Validate a yield curve, compute discount factors, and optionally create a plot
#'
#' Builds a tibble-first representation of a discrete yield curve using
#' annual effective spot rates and computes the corresponding discount factors.
#'
#' Each row is treated as one curve (one case). For tibble input,
#' \code{col_term} and \code{col_spot} must be list-columns of equal-length
#' numeric vectors. When \code{.data = NULL}, \code{term} and \code{spot}
#' must be numeric vectors and a one-row tibble is returned.
#'
#' The discount factors are computed as:
#' \deqn{v_t = (1 + i_t)^{-t}}{v_t = (1 + i_t)^(-t)}
#' where \eqn{i_t}{i_t} is the annual effective spot rate for maturity \eqn{t}.
#'
#' If \code{plot = TRUE}, the function also returns a list-column of
#' \code{ggplot2} objects showing the spot yield curve for each row.
#'
#' @param .data A data.frame or tibble. If \code{NULL}, \code{term} and
#'   \code{spot} must be supplied as numeric vectors.
#' @param term Numeric vector of maturities when \code{.data = NULL}.
#' @param spot Numeric vector of annual effective spot rates when
#'   \code{.data = NULL}.
#' @param col_term Name of the list-column containing maturities.
#' @param col_spot Name of the list-column containing spot rates.
#' @param plot Logical; if \code{TRUE}, adds a list-column of
#'   \code{ggplot2} objects.
#' @param .out Name of the output list-column containing discount factors.
#' @param .out_plot Name of the output list-column containing \code{ggplot2}
#'   objects. Used only if \code{plot = TRUE}.
#' @param .keep One of \code{"all"}, \code{"used"}, or \code{"none"}.
#' @param .na NA handling policy: \code{"propagate"}, \code{"error"}, or
#'   \code{"drop"}.
#'
#' @return A tibble. By default it returns the original columns plus a new
#'   list-column named by \code{.out} containing discount-factor vectors.
#'   If \code{plot = TRUE}, it also adds a list-column named by \code{.out_plot}
#'   containing \code{ggplot2} objects.
#'
#' @seealso \code{\link{forward_rate_tbl}},
#'   \code{\link{discount_factor_spot}}, \code{\link{standardize_interest}}
#'
#' @family interest
#'
#' @examples
#' # Simple example
#' res <- yield_curve_tbl(
#'   term = c(1, 2, 3, 4, 5),
#'   spot = c(0.040, 0.045, 0.048, 0.050, 0.051),
#'   plot = TRUE
#' )
#'
#' res$yield_curve_plot[[1]]
#'
#' # Medium example
#' curves <- tibble::tibble(
#'   curve_id = c("A", "B"),
#'   term = list(c(1, 2, 3), c(1, 3, 5)),
#'   spot = list(c(0.04, 0.05, 0.06), c(0.03, 0.035, 0.04))
#' )
#'
#' res2 <- yield_curve_tbl(
#'   curves,
#'   col_term = "term",
#'   col_spot = "spot",
#'   plot = TRUE,
#'   .out = "v",
#'   .out_plot = "curve_plot"
#' )
#'
#' res2$curve_plot[[2]]
#'
#' @references
#' Marcel B. Finan, *A Basic Course in the Theory of Interest and Derivatives
#' Markets: A Preparation for the Actuarial Exam FM/2*, Section 53:
#' The Term Structure of Interest Rates and Yield Curves.
#'
#' Kellison, S. G. *The Theory of Interest*.
#'
#' @export
yield_curve_tbl <- function(
    .data = NULL,
    term = NULL,
    spot = NULL,
    col_term = "term",
    col_spot = "spot",
    plot = FALSE,
    .out = "discount",
    .out_plot = "yield_curve_plot",
    .keep = c("all", "used", "none"),
    .na = c("propagate", "error", "drop")
) {
  .keep <- match.arg(.keep)
  .na <- match.arg(.na)

  abort <- function(msg) rlang::abort(msg)

  if (!is.logical(plot) || length(plot) != 1L || is.na(plot)) {
    abort("`plot` must be TRUE or FALSE.")
  }

  if (!is.character(.out) || length(.out) != 1L || is.na(.out) || !nzchar(.out)) {
    abort("`.out` must be a single non-empty string.")
  }

  if (!is.character(.out_plot) || length(.out_plot) != 1L || is.na(.out_plot) || !nzchar(.out_plot)) {
    abort("`.out_plot` must be a single non-empty string.")
  }

  if (.out == .out_plot) {
    abort("`.out` and `.out_plot` must have different names.")
  }

  # --- Build data_in from scalar or tibble input ---
  if (is.null(.data)) {
    if (is.null(term) || is.null(spot)) {
      abort("When `.data = NULL`, both `term` and `spot` must be supplied.")
    }

    if (!is.numeric(term)) {
      abort("`term` must be a numeric vector when `.data = NULL`.")
    }

    if (!is.numeric(spot)) {
      abort("`spot` must be a numeric vector when `.data = NULL`.")
    }

    data_in <- tibble::tibble(
      term = list(term),
      spot = list(spot)
    )

    term_name <- "term"
    spot_name <- "spot"
  } else {
    if (!inherits(.data, "data.frame")) {
      abort("`.data` must be a data.frame or tibble.")
    }

    data_in <- tibble::as_tibble(.data)

    if (!col_term %in% names(data_in)) {
      abort(paste0(
        "Column `", col_term, "` was not found. ",
        "Pass the correct column name via `col_term`."
      ))
    }

    if (!col_spot %in% names(data_in)) {
      abort(paste0(
        "Column `", col_spot, "` was not found. ",
        "Pass the correct column name via `col_spot`."
      ))
    }

    term_name <- col_term
    spot_name <- col_spot
  }

  # --- Extract and validate column types ---
  term_col <- data_in[[term_name]]
  spot_col <- data_in[[spot_name]]

  if (!is.list(term_col)) {
    abort(paste0(
      "`", term_name, "` must be a list-column ",
      "(one numeric vector of maturities per row)."
    ))
  }

  if (!is.list(spot_col)) {
    abort(paste0(
      "`", spot_name, "` must be a list-column ",
      "(one numeric vector of spot rates per row)."
    ))
  }

  n <- nrow(data_in)

  # --- NA detection ---
  bad_na <- vapply(seq_len(n), function(i) {
    anyNA(term_col[[i]]) || anyNA(spot_col[[i]])
  }, logical(1))

  if (.na == "error" && any(bad_na)) {
    idx <- which(bad_na)
    abort(paste0(
      "Missing values found in required inputs at row(s): ",
      paste(idx, collapse = ", "), "."
    ))
  }

  if (.na == "drop" && any(bad_na)) {
    warning(
      paste0(
        "Dropping row(s) with missing values in required inputs: ",
        paste(which(bad_na), collapse = ", "), "."
      ),
      call. = FALSE
    )
    data_in <- data_in[!bad_na, , drop = FALSE]
    term_col <- term_col[!bad_na]
    spot_col <- spot_col[!bad_na]
    n <- nrow(data_in)
    bad_na <- rep(FALSE, n)
  }

  # --- Per-row validation ---
  validate_curve <- function(ti, si, row_id) {
    if (!is.numeric(ti)) {
      abort(paste0("`", term_name, "` must contain numeric vectors. Problem at row ", row_id, "."))
    }
    if (!is.numeric(si)) {
      abort(paste0("`", spot_name, "` must contain numeric vectors. Problem at row ", row_id, "."))
    }
    if (length(ti) == 0L) {
      abort(paste0("Empty curve detected at row ", row_id, "."))
    }
    if (length(ti) != length(si)) {
      abort(paste0(
        "`", term_name, "` and `", spot_name,
        "` must have the same length at row ", row_id, "."
      ))
    }
    if (!anyNA(ti) && any(!is.finite(ti))) {
      abort(paste0("`", term_name, "` must be finite at row ", row_id, "."))
    }
    if (!anyNA(si) && any(!is.finite(si))) {
      abort(paste0("`", spot_name, "` must be finite at row ", row_id, "."))
    }
    if (!anyNA(ti) && any(ti < 0)) {
      abort(paste0("All maturities in `", term_name, "` must satisfy t >= 0. Problem at row ", row_id, "."))
    }
    if (!anyNA(si) && any(si <= -1)) {
      abort(paste0("All spot rates in `", spot_name, "` must satisfy i > -1. Problem at row ", row_id, "."))
    }
    if (!anyNA(ti) && length(ti) >= 2L && any(diff(ti) <= 0)) {
      abort(paste0("`", term_name, "` must be strictly increasing without duplicates at row ", row_id, "."))
    }
    invisible(TRUE)
  }

  # --- Plot helper (uses .data pronoun to avoid NSE NOTEs) ---
  make_plot <- function(ti, si) {
    plot_df <- tibble::tibble(term = ti, spot = si)

    ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$term, y = .data$spot)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(
        title = "Yield curve",
        x = "Term",
        y = "Spot rate"
      ) +
      ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 0.01)) +
      ggplot2::theme_minimal()
  }

  # --- Main computation ---
  out_discount <- vector("list", n)
  out_plot_list <- if (plot) vector("list", n) else NULL

  for (i in seq_len(n)) {
    ti <- term_col[[i]]
    si <- spot_col[[i]]

    validate_curve(ti, si, i)

    if (.na == "propagate" && (anyNA(ti) || anyNA(si))) {
      out_discount[[i]] <- rep(NA_real_, length(ti))
      if (plot) out_plot_list[[i]] <- NULL
    } else {
      out_discount[[i]] <- (1 + si)^(-ti)
      if (plot) out_plot_list[[i]] <- make_plot(ti, si)
    }
  }

  # --- Build output ---
  if (.keep == "all") {
    out <- data_in
  } else if (.keep == "used") {
    out <- data_in[, c(term_name, spot_name), drop = FALSE]
  } else {
    out <- tibble::tibble()
  }

  out[[.out]] <- out_discount

  if (plot) {
    out[[.out_plot]] <- out_plot_list
  }

  tibble::as_tibble(out)
}
