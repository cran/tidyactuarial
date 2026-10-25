#' Construct a forward-rate curve from spot rates
#'
#' @md
#' @description
#' Computes the annual effective forward rates for every consecutive interval
#' implied by an observed spot yield curve. Unlike [forward_rate()], which
#' evaluates a selected interval, `forward_curve()` returns the complete
#' adjacent-interval schedule in a tidy, long-format tibble.
#'
#' @details
#' Let \eqn{0=t_0<t_1<\cdots<t_n} denote maturities, and let \eqn{v(t_j)}
#' be their spot discount factors. For each consecutive interval,
#' \deqn{(1+f_{t_{j-1},t_j})^{t_j-t_{j-1}}=
#'       \frac{v(t_{j-1})}{v(t_j)},}
#' so the annual effective forward rate is
#' \deqn{f_{t_{j-1},t_j}=\left(
#'       \frac{v(t_{j-1})}{v(t_j)}\right)^{1/(t_j-t_{j-1})}-1.}
#' The function evaluates log discount-factor ratios from original log spot
#' growth rates, retaining small rates even when reported discounts round to 1.
#' Equal spot log rates are preserved across intervals with a large common
#' origin. Discount factors are obtained from the
#' existing [discount_factor_spot()] function, which standardizes the chosen
#' interest-rate convention. Negative spot and forward rates are supported
#' when the discount factors are positive and finite.
#'
#' With `include_origin = TRUE`, the first interval starts at time zero and
#' uses the exact identity \eqn{v(0)=1}. A supplied zero-maturity node is
#' retained without duplication. With `include_origin = FALSE`, only
#' intervals between supplied nodes are returned.
#'
#' Maturities are sorted within each curve. Repeated, nonfinite or negative
#' maturities, missing or nonfinite rates, and invalid discount factors are
#' rejected. No interpolation, extrapolation, curve fitting or forward-rate
#' forecast is performed. Intervals may have unequal lengths. Spot-to-forward
#' relationships represent no-arbitrage implications of the supplied
#' discount factors, not predictions of realized future rates.
#'
#' Accepted input shapes are: direct numeric `t` and `i` vectors; a data
#' frame with numeric `t` and `i` columns (one node per row, optionally
#' partitioned using `col_curve`); or a data frame with `t` and `i`
#' list-columns (one complete curve per row), such as [yield_curve()] output.
#' Each curve produces one or more rows. `curve` is a character identifier:
#' an ID from `col_curve`, or the input curve's sequential number if no ID
#' is supplied. If more than one curve is supplied, `i_type` and `m` must be
#' scalars; for a single curve they can instead have one value per node.
#'
#' @param .data Optional data frame or tibble containing spot maturities and
#'   rates, in numeric columns or corresponding list-columns.
#' @param t Numeric vector of maturities in years if `.data = NULL`.
#' @param i Numeric vector of spot rates if `.data = NULL`.
#' @param col_t Name of the maturity column in `.data`.
#' @param col_i Name of the spot-rate column in `.data`.
#' @param col_curve Optional name of a curve-identifier column in `.data`.
#'   Required to distinguish multiple curves stored as numeric node rows;
#'   optional for list-column input.
#' @param i_type Spot-rate convention: `"effective"`, `"nominal_interest"`,
#'   `"nominal_discount"` or `"force"`. Scalar for multiple curves, or
#'   scalar/nodewise for a single curve.
#' @param m Conversion frequency for nominal rates. A positive integer;
#'   scalar for multiple curves, or scalar/nodewise for one curve.
#' @param include_origin Logical scalar. If `TRUE` (default), include the
#'   initial interval from time zero when the first observed node is positive.
#'
#' @return A tibble, with one row per forward interval and columns `curve`,
#'   `t_start`, `t_end`, `duration`, `v_start`, `v_end`, and `f`. `f` is an
#'   annual effective forward rate and `duration` is measured in years.
#'
#' @seealso [forward_rate()], [yield_curve()], [discount_factor_spot()],
#'   [bootstrap_spot_curve()]
#' @family interest
#'
#' @examples
#' # A forward curve for consecutive annual spot maturities.
#' forward_curve(t = 1:3, i = c(0.04, 0.045, 0.05))
#'
#' # The output of yield_curve() can be used directly in a pipe.
#' yield_curve(t = c(1, 2, 3), i = c(0.04, 0.045, 0.05)) |>
#'   forward_curve()
#'
#' # Multiple curves, using tidyverse-friendly list-columns.
#' tibble::tibble(
#'   scenario = c("base", "stress"),
#'   t = list(c(1, 2, 3), c(1, 3, 5)),
#'   i = list(c(0.04, 0.045, 0.05), c(0.05, 0.055, 0.06))
#' ) |>
#'   forward_curve(col_curve = "scenario")
#'
#' @export
forward_curve <- function(
    .data = NULL,
    t = NULL,
    i = NULL,
    col_t = "t",
    col_i = "i",
    col_curve = NULL,
    i_type = "effective",
    m = 1L,
    include_origin = TRUE
) {
  abort <- function(message) stop(message, call. = FALSE)

  is_colname <- function(x) {
    is.character(x) && length(x) == 1L && !is.na(x) && nzchar(x)
  }
  if (!is_colname(col_t) || !is_colname(col_i) || identical(col_t, col_i)) {
    abort("`col_t` and `col_i` must be distinct nonempty column names.")
  }
  if (!is.null(col_curve) &&
      (!is_colname(col_curve) || col_curve %in% c(col_t, col_i))) {
    abort("`col_curve` must be NULL or a distinct nonempty column name.")
  }
  if (!is.logical(include_origin) || length(include_origin) != 1L ||
      is.na(include_origin)) {
    abort("`include_origin` must be TRUE or FALSE.")
  }
  if (!is.character(i_type) || !length(i_type) || anyNA(i_type)) {
    abort("`i_type` must be a nonmissing character vector.")
  }
  if (!is.numeric(m) || !length(m) || anyNA(m) ||
      any(!is.finite(m)) || any(m < 1) || any(m != floor(m))) {
    abort("`m` must contain positive integer conversion frequencies.")
  }

  if (is.null(.data)) {
    if (!is.null(col_curve)) {
      abort("`col_curve` requires `.data`.")
    }
    if (!is.numeric(t) || !is.numeric(i) || is.list(t) || is.list(i)) {
      abort("Supply numeric `t` and `i` vectors or a data frame in `.data`.")
    }
    curve_times <- list(t)
    curve_rates <- list(i)
    curve_ids <- "1"
  } else {
    if (!is.data.frame(.data) || nrow(.data) == 0L) {
      abort("`.data` must be a nonempty data frame.")
    }
    if (!is.null(t) || !is.null(i)) {
      abort("Do not supply `t` or `i` separately with `.data`.")
    }
    required <- c(col_t, col_i, col_curve)
    missing_columns <- setdiff(required, names(.data))
    if (length(missing_columns)) {
      abort(paste0("Missing required column(s): ",
                   paste(sprintf("`%s`", missing_columns), collapse = ", "), "."))
    }
    times_column <- .data[[col_t]]
    rates_column <- .data[[col_i]]
    ids <- if (is.null(col_curve)) NULL else .data[[col_curve]]
    if (!is.null(ids)) {
      if (!is.atomic(ids) || is.list(ids) || anyNA(ids) ||
          any(!nzchar(as.character(ids)))) {
        abort("`col_curve` must contain nonmissing scalar identifiers.")
      }
      ids <- as.character(ids)
    }

    if (is.list(times_column) && is.list(rates_column)) {
      curve_times <- times_column
      curve_rates <- rates_column
      curve_ids <- if (is.null(ids)) as.character(seq_len(nrow(.data))) else ids
      if (anyDuplicated(curve_ids)) {
        abort("Curve identifiers must be unique for list-column input.")
      }
    } else if (is.numeric(times_column) && is.numeric(rates_column) &&
               !is.list(times_column) && !is.list(rates_column)) {
      if (is.null(ids)) {
        curve_times <- list(times_column)
        curve_rates <- list(rates_column)
        curve_ids <- "1"
      } else {
        curve_ids <- unique(ids)
        groups <- lapply(curve_ids, function(id) which(ids == id))
        curve_times <- lapply(groups, function(idx) times_column[idx])
        curve_rates <- lapply(groups, function(idx) rates_column[idx])
      }
    } else {
      abort("`col_t` and `col_i` must both be numeric columns or both be list-columns.")
    }
  }

  n_curves <- length(curve_times)
  if (n_curves > 1L && (length(i_type) != 1L || length(m) != 1L)) {
    abort("For multiple curves, `i_type` and `m` must be scalars.")
  }

  results <- vector("list", n_curves)
  for (row in seq_len(n_curves)) {
    tt <- curve_times[[row]]
    ii <- curve_rates[[row]]
    if (!is.numeric(tt) || !is.null(dim(tt)) ||
        !is.numeric(ii) || !is.null(dim(ii)) ||
        !length(tt) || length(tt) != length(ii) ||
        anyNA(tt) || any(!is.finite(tt)) || any(tt < 0) ||
        anyNA(ii) || any(!is.finite(ii))) {
      abort(sprintf("Curve %s requires equal-length finite numeric maturity/rate vectors; maturities must be nonnegative.",
                    curve_ids[[row]]))
    }
    if (!(length(i_type) %in% c(1L, length(tt))) ||
        !(length(m) %in% c(1L, length(tt)))) {
      abort(sprintf("`i_type` and `m` must be scalars or have one value per node in curve %s.",
                    curve_ids[[row]]))
    }
    ord <- order(tt)
    times <- tt[ord]
    rates <- ii[ord]
    if (length(times) > 1L && any(diff(times) <= 0)) {
      abort(sprintf("Curve %s has duplicate maturity nodes.", curve_ids[[row]]))
    }
    typ <- if (length(i_type) == 1L) i_type else i_type[ord]
    frequency <- if (length(m) == 1L) m else m[ord]

    discounts <- discount_factor_spot(
      t = times, i = rates, i_type = typ, m = frequency
    )
    growth <- .interest_basis(i_type = typ, i = rates, m = frequency)$growth
    if (length(discounts) != length(times) || anyNA(discounts) ||
        any(!is.finite(discounts)) || any(discounts <= 0)) {
      abort(sprintf("Curve %s implies nonpositive or nonfinite discount factors.",
                    curve_ids[[row]]))
    }

    if (include_origin && times[[1L]] > 0) {
      nodes <- c(0, times)
      factors <- c(1, discounts)
      growth <- c(0, growth)
    } else {
      nodes <- times
      factors <- discounts
    }
    if (length(nodes) < 2L) {
      abort(sprintf("Curve %s has fewer than two nodes; no forward interval can be constructed.",
                    curve_ids[[row]]))
    }

    start <- head(nodes, -1L)
    end <- tail(nodes, -1L)
    start_v <- head(factors, -1L)
    end_v <- tail(factors, -1L)
    period <- end - start
    start_growth <- head(growth, -1L)
    end_growth <- tail(growth, -1L)
    forward_growth <- end_growth
    nonzero <- start != 0
    forward_growth[nonzero] <- end_growth[nonzero] + start[nonzero] *
      ((end_growth[nonzero] - start_growth[nonzero]) / period[nonzero])
    forward <- expm1(forward_growth)
    if (any(!is.finite(forward)) || any(forward <= -1)) {
      abort(sprintf("Curve %s produces nonfinite or invalid forward rates.",
                    curve_ids[[row]]))
    }
    results[[row]] <- tibble::tibble(
      curve = rep(curve_ids[[row]], length(period)),
      t_start = start,
      t_end = end,
      duration = period,
      v_start = start_v,
      v_end = end_v,
      f = forward
    )
  }
  tibble::as_tibble(do.call(rbind, results))
}
