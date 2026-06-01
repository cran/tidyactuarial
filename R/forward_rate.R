#' Compute an implied forward rate from a discrete spot curve
#'
#' Returns the annual effective forward rate implied between two maturities from
#' a discrete yield curve stored in tibble-first format, using compact
#' actuarial notation.
#'
#' Each row is treated as one curve. For tibble input, \code{col_t} and
#' \code{col_i} must be list-columns of equal-length numeric vectors, and
#' \code{col_t_start} and \code{col_t_end} must be numeric columns giving the
#' forward interval for each row.
#'
#' The implied forward rate is computed from the standardized annual effective
#' spot curve through:
#' \deqn{(1+i_1)^{t_1}(1+f)^{t_2-t_1}=(1+i_2)^{t_2}}
#' so that
#' \deqn{f_{t_1,t_2} =
#' \left(\frac{(1+i_2)^{t_2}}{(1+i_1)^{t_1}}\right)^{1/(t_2-t_1)} - 1.}
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
#' @param .data A data.frame or tibble. If \code{NULL}, \code{t}, \code{i},
#'   \code{t_start}, and \code{t_end} must be supplied.
#' @param t Numeric vector of maturities in years when \code{.data = NULL}.
#' @param i Numeric vector of spot-rate values when \code{.data = NULL}.
#' @param t_start Numeric scalar giving the start maturity when
#'   \code{.data = NULL}.
#' @param t_end Numeric scalar giving the end maturity when
#'   \code{.data = NULL}.
#' @param col_t Name of the list-column containing maturities.
#' @param col_i Name of the list-column containing spot rates.
#' @param col_t_start Name of the numeric column containing the start maturity.
#' @param col_t_end Name of the numeric column containing the end maturity.
#' @param i_type Character vector indicating the spot-rate type. Allowed values
#'   are \code{"effective"}, \code{"nominal_interest"},
#'   \code{"nominal_discount"}, and \code{"force"}. May have length 1 or the
#'   same length as each curve.
#' @param m Positive integer vector giving the conversion frequency for nominal
#'   spot rates. May have length 1 or the same length as each curve.
#' @param method Spot extraction method: \code{"exact"} or \code{"linear"}.
#' @param plot Logical; if \code{TRUE}, adds a list-column of \code{ggplot2}
#'   objects.
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
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{t} denotes maturity, \code{i} denotes the spot
#' rate, \code{i_type} denotes the interest-rate type, and \code{m} denotes the
#' conversion frequency for nominal spot rates. The output column \code{f}
#' denotes the implied annual effective forward rate.
#'
#' Spot-rate inputs are converted to annual effective form using
#' \code{\link{standardize_interest}} before interpolation and forward-rate
#' calculation.
#'
#' @seealso \code{\link{yield_curve}}, \code{\link{discount_factor_spot}},
#'   \code{\link{standardize_interest}}
#'
#' @family interest
#'
#' @examples
#' # Simple example: exact forward rate
#' forward_rate(
#'   t = c(1, 2, 3, 4, 5),
#'   i = c(0.040, 0.045, 0.048, 0.050, 0.051),
#'   t_start = 2,
#'   t_end = 5
#' )
#'
#' # Interpolated forward rates for multiple curves
#' curves <- tibble::tibble(
#'   curve_id = c("A", "B"),
#'   t = list(c(1, 2, 3, 5), c(1, 3, 5, 7)),
#'   i = list(c(0.04, 0.05, 0.055, 0.06),
#'            c(0.03, 0.035, 0.04, 0.045)),
#'   t_start = c(2, 2),
#'   t_end = c(4, 6)
#' )
#'
#' forward_rate(
#'   curves,
#'   method = "linear",
#'   plot = TRUE
#' )
#'
#' # Nominal annual spot rates convertible semiannually
#' forward_rate(
#'   t = c(1, 2, 3),
#'   i = c(0.05, 0.055, 0.06),
#'   i_type = "nominal_interest",
#'   m = 2,
#'   t_start = 1,
#'   t_end = 3
#' )
#'
#' @references
#' Marcel B. Finan, \emph{A Basic Course in the Theory of Interest and
#' Derivatives Markets: A Preparation for the Actuarial Exam FM/2}, Section 53:
#' The Term Structure of Interest Rates and Yield Curves.
#'
#' Kellison, S. G. \emph{The Theory of Interest}, Chapter 10:
#' The Term Structure of Interest Rates.
#'
#' @export
forward_rate <- function(
    .data = NULL,
    t = NULL,
    i = NULL,
    t_start = NULL,
    t_end = NULL,
    col_t = "t",
    col_i = "i",
    col_t_start = "t_start",
    col_t_end = "t_end",
    i_type = "effective",
    m = 1L,
    method = c("exact", "linear"),
    plot = FALSE,
    .out = "f",
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

  # --- Build data_in from scalar or tibble input ---
  if (is.null(.data)) {
    if (is.null(t) || is.null(i) || is.null(t_start) || is.null(t_end)) {
      abort("When `.data = NULL`, `t`, `i`, `t_start`, and `t_end` must all be supplied.")
    }

    if (!is.numeric(t)) {
      abort("`t` must be a numeric vector when `.data = NULL`.")
    }

    if (!is.numeric(i)) {
      abort("`i` must be a numeric vector when `.data = NULL`.")
    }

    if (!is.numeric(t_start) || length(t_start) != 1L) {
      abort("`t_start` must be a numeric scalar when `.data = NULL`.")
    }

    if (!is.numeric(t_end) || length(t_end) != 1L) {
      abort("`t_end` must be a numeric scalar when `.data = NULL`.")
    }

    data_in <- tibble::tibble(
      t = list(t),
      i = list(i),
      t_start = t_start,
      t_end = t_end
    )

    t_name <- "t"
    i_name <- "i"
    t_start_name <- "t_start"
    t_end_name <- "t_end"
  } else {
    if (!inherits(.data, "data.frame")) {
      abort("`.data` must be a data.frame or tibble.")
    }

    data_in <- tibble::as_tibble(.data)

    required <- c(col_t, col_i, col_t_start, col_t_end)
    missing_cols <- setdiff(required, names(data_in))

    if (length(missing_cols) > 0L) {
      abort(paste0(
        "Missing required column(s): ",
        paste(sprintf("`%s`", missing_cols), collapse = ", "),
        ". Pass the correct names via the corresponding `col_*` arguments."
      ))
    }

    t_name <- col_t
    i_name <- col_i
    t_start_name <- col_t_start
    t_end_name <- col_t_end
  }

  # --- Extract and validate column types ---
  t_col <- data_in[[t_name]]
  i_col <- data_in[[i_name]]
  t_start_col <- data_in[[t_start_name]]
  t_end_col <- data_in[[t_end_name]]

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

  if (!is.numeric(t_start_col)) {
    abort(paste0("`", t_start_name, "` must be a numeric column (one start maturity per row)."))
  }

  if (!is.numeric(t_end_col)) {
    abort(paste0("`", t_end_name, "` must be a numeric column (one end maturity per row)."))
  }

  n_rows <- nrow(data_in)

  # --- NA detection ---
  bad_na <- vapply(seq_len(n_rows), function(row_idx) {
    anyNA(t_col[[row_idx]]) || anyNA(i_col[[row_idx]]) ||
      is.na(t_start_col[[row_idx]]) || is.na(t_end_col[[row_idx]])
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
    t_col <- t_col[!bad_na]
    i_col <- i_col[!bad_na]
    t_start_col <- t_start_col[!bad_na]
    t_end_col <- t_end_col[!bad_na]
    n_rows <- nrow(data_in)
    bad_na <- rep(FALSE, n_rows)
  }

  # --- Per-row validation helper ---
  validate_curve <- function(tt, ii, s0, s1, row_id) {
    if (!is.numeric(tt)) {
      abort(paste0(
        "`", t_name,
        "` must contain numeric vectors. Problem at row ",
        row_id, "."
      ))
    }

    if (!is.numeric(ii)) {
      abort(paste0(
        "`", i_name,
        "` must contain numeric vectors. Problem at row ",
        row_id, "."
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

    if (!is.numeric(s0) || length(s0) != 1L || (!is.na(s0) && !is.finite(s0))) {
      abort(paste0(
        "`", t_start_name,
        "` must be a finite numeric scalar at row ",
        row_id, "."
      ))
    }

    if (!is.numeric(s1) || length(s1) != 1L || (!is.na(s1) && !is.finite(s1))) {
      abort(paste0(
        "`", t_end_name,
        "` must be a finite numeric scalar at row ",
        row_id, "."
      ))
    }

    if (!anyNA(tt) && any(tt < 0)) {
      abort(paste0(
        "All maturities in `", t_name,
        "` must satisfy t >= 0. Problem at row ",
        row_id, "."
      ))
    }

    if (!is.na(s0) && s0 < 0) {
      abort(paste0(
        "Start maturity `", t_start_name,
        "` must satisfy t >= 0. Problem at row ",
        row_id, "."
      ))
    }

    if (!is.na(s1) && s1 <= s0) {
      abort(paste0(
        "`", t_end_name,
        "` must be strictly greater than `",
        t_start_name,
        "` at row ",
        row_id, "."
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

  # --- Spot extraction helpers ---
  extract_exact <- function(tt, ii, target, row_id, target_name, tol = 1e-12) {
    idx <- which(abs(tt - target) <= tol)

    if (length(idx) == 1L) {
      return(ii[idx])
    }

    abort(paste0(
      "Target maturity `", target_name, "` = ", target,
      " does not match any curve node at row ", row_id,
      ". Use `method = \"linear\"` if interpolation is intended."
    ))
  }

  extract_linear <- function(tt, ii, target, row_id, target_name, tol = 1e-12) {
    idx <- which(abs(tt - target) <= tol)

    if (length(idx) == 1L) {
      return(ii[idx])
    }

    if (target < min(tt) || target > max(tt)) {
      abort(paste0(
        "Target maturity `", target_name, "` = ", target,
        " is outside the observed maturity range at row ", row_id,
        ". No extrapolation is performed."
      ))
    }

    left_idx <- max(which(tt < target))
    right_idx <- min(which(tt > target))

    t_left <- tt[left_idx]
    t_right <- tt[right_idx]
    i_left <- ii[left_idx]
    i_right <- ii[right_idx]

    i_left + (target - t_left) * (i_right - i_left) / (t_right - t_left)
  }

  # --- Plot helper ---
  make_plot <- function(tt, ii, s0, s1, i0, i1) {
    plot_df <- tibble::tibble(t = tt, i = ii)

    point_df <- tibble::tibble(
      t = c(s0, s1),
      i = c(i0, i1),
      point_type = c("start", "end")
    )

    ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[["t"]], y = .data[["i"]])) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_point(
        data = point_df,
        ggplot2::aes(x = .data[["t"]], y = .data[["i"]]),
        inherit.aes = FALSE,
        size = 3,
        colour = "red"
      ) +
      ggplot2::geom_segment(
        data = point_df[1, , drop = FALSE],
        ggplot2::aes(
          x = s0,
          xend = s1,
          y = i0,
          yend = i1
        ),
        inherit.aes = FALSE,
        linetype = "dashed"
      ) +
      ggplot2::labs(
        title = "Forward rate interval on spot curve",
        x = "Term t",
        y = "Annual effective spot rate i"
      ) +
      ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 0.01)) +
      ggplot2::theme_minimal()
  }

  # --- Main computation loop ---
  out_forward <- rep(NA_real_, n_rows)
  out_i_start <- rep(NA_real_, n_rows)
  out_i_end <- rep(NA_real_, n_rows)
  out_plot_list <- if (plot) vector("list", n_rows) else NULL

  for (row_idx in seq_len(n_rows)) {
    tt <- t_col[[row_idx]]
    ii <- i_col[[row_idx]]
    s0 <- t_start_col[[row_idx]]
    s1 <- t_end_col[[row_idx]]

    validate_curve(tt, ii, s0, s1, row_idx)
    validate_rate_inputs(ii, row_idx)

    if (.na == "propagate" && (anyNA(tt) || anyNA(ii) || is.na(s0) || is.na(s1))) {
      out_forward[row_idx] <- NA_real_

      if (plot) {
        out_plot_list[[row_idx]] <- NULL
      }

      next
    }

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

    i_start <- switch(
      method,
      exact = extract_exact(tt, ii_eff, s0, row_idx, t_start_name),
      linear = extract_linear(tt, ii_eff, s0, row_idx, t_start_name)
    )

    i_end <- switch(
      method,
      exact = extract_exact(tt, ii_eff, s1, row_idx, t_end_name),
      linear = extract_linear(tt, ii_eff, s1, row_idx, t_end_name)
    )

    fwd <- ((1 + i_end)^s1 / (1 + i_start)^s0)^(1 / (s1 - s0)) - 1

    out_forward[row_idx] <- fwd
    out_i_start[row_idx] <- i_start
    out_i_end[row_idx] <- i_end

    if (plot) {
      out_plot_list[[row_idx]] <- make_plot(tt, ii_eff, s0, s1, i_start, i_end)
    }
  }

  # --- Build output ---
  if (.keep == "all") {
    out <- data_in
  } else if (.keep == "used") {
    out <- data_in[, c(t_name, i_name, t_start_name, t_end_name), drop = FALSE]
  } else {
    out <- tibble::tibble()
  }

  out[[.out]] <- out_forward
  out$i_start <- out_i_start
  out$i_end <- out_i_end

  if (plot) {
    out[[.out_plot]] <- out_plot_list
  }

  tibble::as_tibble(out)
}
