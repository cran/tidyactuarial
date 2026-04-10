#' Compute an implied forward rate from a discrete spot curve
#'
#' Returns the annual effective forward rate implied between two maturities
#' from a discrete yield curve stored in tibble-first format.
#'
#' Each row is treated as one curve (one case). For tibble input,
#' \code{col_term} and \code{col_spot} must be list-columns of equal-length
#' numeric vectors, and \code{col_t_start} and \code{col_t_end} must be
#' numeric columns giving the forward interval for each row.
#'
#' The implied forward rate is computed from the spot curve through:
#' \deqn{(1+i_1)^{t_1}(1+f)^{t_2-t_1}=(1+i_2)^{t_2}}{(1+i1)^t1 * (1+f)^(t2-t1) = (1+i2)^t2}
#' so that
#' \deqn{f_{t_1,t_2} = \left(\frac{(1+i_2)^{t_2}}{(1+i_1)^{t_1}}\right)^{1/(t_2-t_1)} - 1}{f = ((1+i2)^t2 / (1+i1)^t1)^(1/(t2-t1)) - 1}
#'
#' Two extraction methods are supported for the spot rates:
#' \itemize{
#'   \item \code{"exact"}: requires that \code{t_start} and \code{t_end}
#'     match curve nodes.
#'   \item \code{"linear"}: uses linear interpolation between adjacent nodes.
#' }
#'
#' No extrapolation is performed outside the observed maturity range.
#'
#' @param .data A data.frame or tibble. If \code{NULL}, \code{term},
#'   \code{spot}, \code{t_start}, and \code{t_end} must be supplied.
#' @param term Numeric vector of maturities when \code{.data = NULL}.
#' @param spot Numeric vector of annual effective spot rates when
#'   \code{.data = NULL}.
#' @param t_start Numeric scalar giving the start maturity when
#'   \code{.data = NULL}.
#' @param t_end Numeric scalar giving the end maturity when
#'   \code{.data = NULL}.
#' @param col_term Name of the list-column containing maturities.
#' @param col_spot Name of the list-column containing spot rates.
#' @param col_t_start Name of the numeric column containing the start maturity.
#' @param col_t_end Name of the numeric column containing the end maturity.
#' @param method Spot extraction method: \code{"exact"} or \code{"linear"}.
#' @param plot Logical; if \code{TRUE}, adds a list-column of
#'   \code{ggplot2} objects.
#' @param .out Name of the output column containing the forward rate.
#' @param .out_plot Name of the output list-column containing \code{ggplot2}
#'   objects. Used only if \code{plot = TRUE}.
#' @param .keep One of \code{"all"}, \code{"used"}, or \code{"none"}.
#' @param .na NA handling policy: \code{"propagate"}, \code{"error"}, or
#'   \code{"drop"}.
#'
#' @return A tibble. By default it returns the original columns plus a new
#'   numeric column named by \code{.out}. If \code{plot = TRUE}, it also adds
#'   a list-column named by \code{.out_plot} containing \code{ggplot2} objects.
#'
#' @seealso \code{yield_curve}, \code{\link{discount_factor_spot}},
#'   \code{\link{standardize_interest}}
#'
#' @family interest
#'
#' @examples
#' # Simple example: exact forward rate
#' forward_rate_tbl(
#'   term = c(1, 2, 3, 4, 5),
#'   spot = c(0.040, 0.045, 0.048, 0.050, 0.051),
#'   t_start = 2,
#'   t_end = 5
#' )
#'
#' # Medium example: interpolated forward rates for multiple curves
#' curves <- tibble::tibble(
#'   curve_id = c("A", "B"),
#'   term = list(c(1, 2, 3, 5), c(1, 3, 5, 7)),
#'   spot = list(c(0.04, 0.05, 0.055, 0.06),
#'               c(0.03, 0.035, 0.04, 0.045)),
#'   t_start = c(2, 2),
#'   t_end = c(4, 6)
#' )
#'
#' forward_rate_tbl(
#'   curves,
#'   method = "linear",
#'   plot = TRUE
#' )
#'
#' @references
#' Marcel B. Finan, *A Basic Course in the Theory of Interest and Derivatives
#' Markets: A Preparation for the Actuarial Exam FM/2*, Section 53:
#' The Term Structure of Interest Rates and Yield Curves.
#'
#' Kellison, S. G. *The Theory of Interest*, Chapter 10:
#' The Term Structure of Interest Rates.
#'
#' @export
forward_rate_tbl <- function(
    .data = NULL,
    term = NULL,
    spot = NULL,
    t_start = NULL,
    t_end = NULL,
    col_term = "term",
    col_spot = "spot",
    col_t_start = "t_start",
    col_t_end = "t_end",
    method = c("exact", "linear"),
    plot = FALSE,
    .out = "forward_rate",
    .out_plot = "forward_rate_plot",
    .keep = c("all", "used", "none"),
    .na = c("propagate", "error", "drop")
) {
  method <- match.arg(method)
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

  if (plot && identical(.out, .out_plot)) {
    abort("`.out` and `.out_plot` must have different names.")
  }

  # --- Build data_in from scalar or tibble input ---
  if (is.null(.data)) {
    if (is.null(term) || is.null(spot) || is.null(t_start) || is.null(t_end)) {
      abort("When `.data = NULL`, `term`, `spot`, `t_start`, and `t_end` must all be supplied.")
    }

    if (!is.numeric(term)) {
      abort("`term` must be a numeric vector when `.data = NULL`.")
    }

    if (!is.numeric(spot)) {
      abort("`spot` must be a numeric vector when `.data = NULL`.")
    }

    if (!is.numeric(t_start) || length(t_start) != 1L) {
      abort("`t_start` must be a numeric scalar when `.data = NULL`.")
    }

    if (!is.numeric(t_end) || length(t_end) != 1L) {
      abort("`t_end` must be a numeric scalar when `.data = NULL`.")
    }

    data_in <- tibble::tibble(
      term = list(term),
      spot = list(spot),
      t_start = t_start,
      t_end = t_end
    )

    term_name <- "term"
    spot_name <- "spot"
    t_start_name <- "t_start"
    t_end_name <- "t_end"
  } else {
    if (!inherits(.data, "data.frame")) {
      abort("`.data` must be a data.frame or tibble.")
    }

    data_in <- tibble::as_tibble(.data)

    required <- c(col_term, col_spot, col_t_start, col_t_end)
    missing_cols <- setdiff(required, names(data_in))
    if (length(missing_cols) > 0L) {
      abort(paste0(
        "Missing required column(s): ",
        paste(sprintf("`%s`", missing_cols), collapse = ", "),
        ". Pass the correct names via the corresponding `col_*` arguments."
      ))
    }

    term_name <- col_term
    spot_name <- col_spot
    t_start_name <- col_t_start
    t_end_name <- col_t_end
  }

  # --- Extract and validate column types ---
  term_col <- data_in[[term_name]]
  spot_col <- data_in[[spot_name]]
  t_start_col <- data_in[[t_start_name]]
  t_end_col <- data_in[[t_end_name]]

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

  if (!is.numeric(t_start_col)) {
    abort(paste0("`", t_start_name, "` must be a numeric column (one start maturity per row)."))
  }

  if (!is.numeric(t_end_col)) {
    abort(paste0("`", t_end_name, "` must be a numeric column (one end maturity per row)."))
  }

  n <- nrow(data_in)

  # --- NA detection ---
  bad_na <- vapply(seq_len(n), function(i) {
    anyNA(term_col[[i]]) || anyNA(spot_col[[i]]) ||
      is.na(t_start_col[[i]]) || is.na(t_end_col[[i]])
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
    t_start_col <- t_start_col[!bad_na]
    t_end_col <- t_end_col[!bad_na]
    n <- nrow(data_in)
    bad_na <- rep(FALSE, n)
  }

  # --- Per-row validation helper ---
  validate_curve <- function(ti, si, s0, s1, row_id) {
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
    if (!is.numeric(s0) || length(s0) != 1L || (!is.na(s0) && !is.finite(s0))) {
      abort(paste0("`", t_start_name, "` must be a finite numeric scalar at row ", row_id, "."))
    }
    if (!is.numeric(s1) || length(s1) != 1L || (!is.na(s1) && !is.finite(s1))) {
      abort(paste0("`", t_end_name, "` must be a finite numeric scalar at row ", row_id, "."))
    }
    if (!anyNA(ti) && any(ti < 0)) {
      abort(paste0("All maturities in `", term_name, "` must satisfy t >= 0. Problem at row ", row_id, "."))
    }
    if (!is.na(s0) && s0 < 0) {
      abort(paste0("Start maturity `", t_start_name, "` must satisfy t >= 0. Problem at row ", row_id, "."))
    }
    if (!is.na(s1) && s1 <= s0) {
      abort(paste0("`", t_end_name, "` must be strictly greater than `", t_start_name, "` at row ", row_id, "."))
    }
    if (!anyNA(si) && any(si <= -1)) {
      abort(paste0("All spot rates in `", spot_name, "` must satisfy i > -1. Problem at row ", row_id, "."))
    }
    if (!anyNA(ti) && length(ti) >= 2L && any(diff(ti) <= 0)) {
      abort(paste0("`", term_name, "` must be strictly increasing without duplicates at row ", row_id, "."))
    }
    invisible(TRUE)
  }

  # --- Spot extraction helpers ---
  extract_exact <- function(ti, si, target, row_id, target_name, tol = 1e-12) {
    idx <- which(abs(ti - target) <= tol)
    if (length(idx) == 1L) return(si[idx])
    abort(paste0(
      "Target maturity `", target_name, "` = ", target,
      " does not match any curve node at row ", row_id,
      ". Use `method = \"linear\"` if interpolation is intended."
    ))
  }

  extract_linear <- function(ti, si, target, row_id, target_name, tol = 1e-12) {
    idx <- which(abs(ti - target) <= tol)
    if (length(idx) == 1L) return(si[idx])

    if (target < min(ti) || target > max(ti)) {
      abort(paste0(
        "Target maturity `", target_name, "` = ", target,
        " is outside the observed maturity range at row ", row_id,
        ". No extrapolation is performed."
      ))
    }

    j <- max(which(ti < target))
    k <- min(which(ti > target))
    t_left <- ti[j]; t_right <- ti[k]
    i_left <- si[j]; i_right <- si[k]
    i_left + (target - t_left) * (i_right - i_left) / (t_right - t_left)
  }

  # --- Plot helper (uses .data pronoun to avoid NSE NOTEs) ---
  make_plot <- function(ti, si, s0, s1, i0, i1) {
    plot_df <- tibble::tibble(term = ti, spot = si)
    point_df <- tibble::tibble(
      term = c(s0, s1),
      spot = c(i0, i1),
      point_type = c("start", "end")
    )

    ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$term, y = .data$spot)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_point(
        data = point_df,
        ggplot2::aes(x = .data$term, y = .data$spot),
        inherit.aes = FALSE,
        size = 3,
        colour = "red"
      ) +
      ggplot2::geom_segment(
        data = point_df[1, , drop = FALSE],
        ggplot2::aes(
          x = s0, xend = s1,
          y = i0, yend = i1
        ),
        inherit.aes = FALSE,
        linetype = "dashed"
      ) +
      ggplot2::labs(
        title = "Forward rate interval on spot curve",
        x = "Term",
        y = "Spot rate"
      ) +
      ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 0.01)) +
      ggplot2::theme_minimal()
  }

  # --- Main computation loop ---
  out_forward <- rep(NA_real_, n)
  out_plot_list <- if (plot) vector("list", n) else NULL

  for (i in seq_len(n)) {
    ti <- term_col[[i]]
    si <- spot_col[[i]]
    s0 <- t_start_col[[i]]
    s1 <- t_end_col[[i]]

    validate_curve(ti, si, s0, s1, i)

    if (.na == "propagate" && (anyNA(ti) || anyNA(si) || is.na(s0) || is.na(s1))) {
      out_forward[i] <- NA_real_
      if (plot) out_plot_list[[i]] <- NULL
      next
    }

    i_start <- switch(
      method,
      exact = extract_exact(ti, si, s0, i, t_start_name),
      linear = extract_linear(ti, si, s0, i, t_start_name)
    )

    i_end <- switch(
      method,
      exact = extract_exact(ti, si, s1, i, t_end_name),
      linear = extract_linear(ti, si, s1, i, t_end_name)
    )

    fwd <- ((1 + i_end)^s1 / (1 + i_start)^s0)^(1 / (s1 - s0)) - 1
    out_forward[i] <- fwd

    if (plot) {
      out_plot_list[[i]] <- make_plot(ti, si, s0, s1, i_start, i_end)
    }
  }

  # --- Build output ---
  if (.keep == "all") {
    out <- data_in
  } else if (.keep == "used") {
    out <- data_in[, c(term_name, spot_name, t_start_name, t_end_name), drop = FALSE]
  } else {
    out <- tibble::tibble()
  }

  out[[.out]] <- out_forward

  if (plot) {
    out[[.out_plot]] <- out_plot_list
  }

  tibble::as_tibble(out)
}
