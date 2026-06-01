#' Validate a yield curve and compute discount factors
#'
#' Builds a tibble-first representation of a discrete yield curve and computes
#' the corresponding spot discount factors, using compact actuarial notation.
#'
#' Each row is treated as one curve. For tibble input, \code{col_t} and
#' \code{col_i} must identify list-columns of equal-length numeric vectors.
#' When \code{.data = NULL}, \code{t} and \code{i} must be numeric vectors and
#' a one-row tibble is returned.
#'
#' The discount factors are computed as:
#' \deqn{v_t = (1+i_t)^{-t}}
#' where \eqn{i_t} is the annual effective spot rate for maturity \eqn{t}.
#'
#' If \code{plot = TRUE}, the function also returns a list-column of
#' \code{ggplot2} objects showing the spot yield curve for each row.
#'
#' @param .data A data frame or tibble. If \code{NULL}, \code{t} and \code{i}
#'   must be supplied as numeric vectors.
#' @param t Numeric vector of maturities in years when \code{.data = NULL}.
#' @param i Numeric vector of spot-rate values when \code{.data = NULL}.
#' @param col_t Name of the list-column containing maturities.
#' @param col_i Name of the list-column containing spot rates.
#' @param i_type Character vector indicating the spot-rate type. Allowed values
#'   are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}. May have length 1 or the
#'   same length as each curve.
#' @param m Positive integer vector giving the conversion frequency for nominal
#'   spot rates. May have length 1 or the same length as each curve.
#' @param plot Logical. If \code{TRUE}, adds a list-column of \code{ggplot2}
#'   objects.
#' @param .out Name of the output list-column containing discount factors.
#' @param .out_plot Name of the output list-column containing \code{ggplot2}
#'   objects. Used only if \code{plot = TRUE}.
#' @param .keep One of \code{"all"}, \code{"used"}, or \code{"none"}.
#' @param .na Missing-value handling policy: \code{"propagate"},
#'   \code{"error"}, or \code{"drop"}.
#'
#' @return A tibble. By default, it returns the original columns plus a new
#'   list-column named by \code{.out} containing discount-factor vectors. If
#'   \code{plot = TRUE}, it also adds a list-column named by \code{.out_plot}
#'   containing \code{ggplot2} objects.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{t} denotes maturity, \code{i} denotes the spot
#' rate, \code{i_type} denotes the interest-rate type, and \code{m} denotes the
#' conversion frequency for nominal spot rates.
#'
#' The spot-rate input is converted to annual effective form through
#' \code{\link{standardize_interest}} before discount factors are computed.
#'
#' @seealso \code{\link{forward_rate}},
#'   \code{\link{discount_factor_spot}}, \code{\link{standardize_interest}}
#'
#' @family interest
#'
#' @examples
#' # Simple example
#' res <- yield_curve(
#'   t = c(1, 2, 3, 4, 5),
#'   i = c(0.040, 0.045, 0.048, 0.050, 0.051),
#'   plot = TRUE
#' )
#'
#' res$yield_curve_plot[[1]]
#'
#' # Multiple curves in a tibble
#' curves <- tibble::tibble(
#'   curve_id = c("A", "B"),
#'   t = list(c(1, 2, 3), c(1, 3, 5)),
#'   i = list(c(0.04, 0.05, 0.06), c(0.03, 0.035, 0.04))
#' )
#'
#' res2 <- yield_curve(
#'   curves,
#'   col_t = "t",
#'   col_i = "i",
#'   plot = TRUE,
#'   .out = "v",
#'   .out_plot = "curve_plot"
#' )
#'
#' res2$curve_plot[[2]]
#'
#' # Nominal annual spot rates convertible semiannually
#' yield_curve(
#'   t = c(1, 2, 3),
#'   i = c(0.05, 0.055, 0.06),
#'   i_type = "nominal_interest",
#'   m = 2
#' )
#'
#' @references
#' Marcel B. Finan, \emph{A Basic Course in the Theory of Interest and
#' Derivatives Markets: A Preparation for the Actuarial Exam FM/2}, Section 53:
#' The Term Structure of Interest Rates and Yield Curves.
#'
#' Kellison, S. G. \emph{The Theory of Interest}.
#'
#' @export
yield_curve <- function(
    .data = NULL,
    t = NULL,
    i = NULL,
    col_t = "t",
    col_i = "i",
    i_type = "effective",
    m = 1L,
    plot = FALSE,
    .out = "v",
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

  if (!is.character(i_type)) {
    abort("`i_type` must be a character vector.")
  }

  if (!is.numeric(m)) {
    abort("`m` must be numeric.")
  }

  if (any(is.na(m)) || any(!is.finite(m)) ||
      any(m < 1) || any(abs(m - round(m)) > 1e-10)) {
    abort("`m` must contain positive integer values.")
  }

  # --- Build data_in from vector or tibble input ---
  if (is.null(.data)) {
    if (is.null(t) || is.null(i)) {
      abort("When `.data = NULL`, both `t` and `i` must be supplied.")
    }

    if (!is.numeric(t)) {
      abort("`t` must be a numeric vector when `.data = NULL`.")
    }

    if (!is.numeric(i)) {
      abort("`i` must be a numeric vector when `.data = NULL`.")
    }

    data_in <- tibble::tibble(
      t = list(t),
      i = list(i)
    )

    t_name <- "t"
    i_name <- "i"
  } else {
    if (!inherits(.data, "data.frame")) {
      abort("`.data` must be a data.frame or tibble.")
    }

    data_in <- tibble::as_tibble(.data)

    if (!is.character(col_t) || length(col_t) != 1L || is.na(col_t) || !nzchar(col_t)) {
      abort("`col_t` must be a single non-empty string.")
    }

    if (!is.character(col_i) || length(col_i) != 1L || is.na(col_i) || !nzchar(col_i)) {
      abort("`col_i` must be a single non-empty string.")
    }

    if (!col_t %in% names(data_in)) {
      abort(paste0(
        "Column `", col_t, "` was not found. ",
        "Pass the correct column name via `col_t`."
      ))
    }

    if (!col_i %in% names(data_in)) {
      abort(paste0(
        "Column `", col_i, "` was not found. ",
        "Pass the correct column name via `col_i`."
      ))
    }

    t_name <- col_t
    i_name <- col_i
  }

  # --- Extract and validate column types ---
  t_col <- data_in[[t_name]]
  i_col <- data_in[[i_name]]

  if (!is.list(t_col)) {
    abort(paste0(
      "`", t_name, "` must be a list-column ",
      "(one numeric vector of maturities per row)."
    ))
  }

  if (!is.list(i_col)) {
    abort(paste0(
      "`", i_name, "` must be a list-column ",
      "(one numeric vector of spot rates per row)."
    ))
  }

  n_rows <- nrow(data_in)

  # --- NA detection ---
  bad_na <- vapply(seq_len(n_rows), function(row_idx) {
    anyNA(t_col[[row_idx]]) || anyNA(i_col[[row_idx]])
  }, logical(1L))

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
    t_col <- t_col[!bad_na]
    i_col <- i_col[!bad_na]
    n_rows <- nrow(data_in)
    bad_na <- rep(FALSE, n_rows)
  }

  # --- Per-row validation ---
  validate_curve <- function(tt, ii, row_id) {
    if (!is.numeric(tt)) {
      abort(paste0(
        "`", t_name, "` must contain numeric vectors. ",
        "Problem at row ", row_id, "."
      ))
    }

    if (!is.numeric(ii)) {
      abort(paste0(
        "`", i_name, "` must contain numeric vectors. ",
        "Problem at row ", row_id, "."
      ))
    }

    if (length(tt) == 0L) {
      abort(paste0("Empty curve detected at row ", row_id, "."))
    }

    if (length(tt) != length(ii)) {
      abort(paste0(
        "`", t_name, "` and `", i_name,
        "` must have the same length at row ", row_id, "."
      ))
    }

    if (!anyNA(tt) && any(!is.finite(tt))) {
      abort(paste0("`", t_name, "` must be finite at row ", row_id, "."))
    }

    if (!anyNA(ii) && any(!is.finite(ii))) {
      abort(paste0("`", i_name, "` must be finite at row ", row_id, "."))
    }

    if (!anyNA(tt) && any(tt < 0)) {
      abort(paste0(
        "All maturities in `", t_name,
        "` must satisfy t >= 0. Problem at row ", row_id, "."
      ))
    }

    if (!anyNA(tt) && length(tt) >= 2L && any(diff(tt) <= 0)) {
      abort(paste0(
        "`", t_name,
        "` must be strictly increasing without duplicates at row ",
        row_id, "."
      ))
    }

    invisible(TRUE)
  }

  validate_rate_inputs <- function(ii, row_id) {
    if (!length(i_type) %in% c(1L, length(ii))) {
      abort(paste0(
        "`i_type` must have length 1 or the same length as the spot-rate vector. ",
        "Problem at row ", row_id, "."
      ))
    }

    if (!length(m) %in% c(1L, length(ii))) {
      abort(paste0(
        "`m` must have length 1 or the same length as the spot-rate vector. ",
        "Problem at row ", row_id, "."
      ))
    }

    invisible(TRUE)
  }

  # --- Plot helper ---
  make_plot <- function(tt, ii_eff) {
    plot_df <- tibble::tibble(t = tt, i = ii_eff)

    ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[["t"]], y = .data[["i"]])) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(
        title = "Yield curve",
        x = "Term t",
        y = "Annual effective spot rate i"
      ) +
      ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 0.01)) +
      ggplot2::theme_minimal()
  }

  # --- Main computation ---
  out_discount <- vector("list", n_rows)
  out_i_effective <- vector("list", n_rows)
  out_plot_list <- if (plot) vector("list", n_rows) else NULL

  for (row_idx in seq_len(n_rows)) {
    tt <- t_col[[row_idx]]
    ii <- i_col[[row_idx]]

    validate_curve(tt, ii, row_idx)
    validate_rate_inputs(ii, row_idx)

    if (.na == "propagate" && (anyNA(tt) || anyNA(ii))) {
      out_discount[[row_idx]] <- rep(NA_real_, length(tt))
      out_i_effective[[row_idx]] <- rep(NA_real_, length(tt))

      if (plot) {
        out_plot_list[[row_idx]] <- NULL
      }
    } else {
      ii_eff <- standardize_interest(
        i_type = rep_len(i_type, length(ii)),
        i = ii,
        m = rep_len(m, length(ii))
      )

      if (any(!is.finite(ii_eff) | ii_eff <= -1)) {
        abort(paste0(
          "The standardized annual effective spot rates must be greater than -1. ",
          "Problem at row ", row_idx, "."
        ))
      }

      out_i_effective[[row_idx]] <- ii_eff
      out_discount[[row_idx]] <- (1 + ii_eff)^(-tt)

      if (plot) {
        out_plot_list[[row_idx]] <- make_plot(tt, ii_eff)
      }
    }
  }

  # --- Build output ---
  if (.keep == "all") {
    out <- data_in
  } else if (.keep == "used") {
    out <- data_in[, c(t_name, i_name), drop = FALSE]
  } else {
    out <- tibble::tibble()
  }

  out[[.out]] <- out_discount
  out$i_effective <- out_i_effective

  if (plot) {
    out[[.out_plot]] <- out_plot_list
  }

  tibble::as_tibble(out)
}
