#' Multiple internal rates of return for a cash flow
#'
#' Searches for multiple internal rates of return (IRRs) of a cash flow by
#' scanning a search interval and solving for all detectable roots of the net
#' present value (NPV) function, using compact actuarial notation.
#'
#' This function is intended for cash flows with multiple sign changes, where
#' more than one IRR may exist. It evaluates the NPV on a fine grid over the
#' search interval, identifies subintervals with sign changes (and grid points
#' where the NPV is approximately zero), and applies
#' \code{\link[stats]{uniroot}} to each candidate interval.
#'
#' The IRRs returned are interpreted as annual effective rates.
#'
#' Timing can be supplied either through \code{t} (in years) or \code{date}
#' (calendar dates). If \code{date} is supplied, the earliest date is treated
#' as time 0.
#'
#' @param cf Numeric vector of cash flows.
#' @param t Optional numeric vector of cash-flow times in years.
#' @param date Optional vector of cash-flow dates. If supplied, the earliest
#'   date is treated as time 0.
#' @param m Positive integer used only to report an equivalent nominal annual
#'   interest rate convertible \code{m} times per year.
#' @param search_interval Numeric vector of length 2 giving the search interval
#'   for annual effective IRRs. Default is \code{c(-0.99, 10)}.
#' @param grid_points Positive integer giving the number of grid points used
#'   to scan the interval. Larger values improve detection at the cost of speed.
#' @param tol Numeric tolerance passed to \code{\link[stats]{uniroot}}.
#' @param maxiter Positive integer passed to \code{\link[stats]{uniroot}}.
#' @param day_count Day-count convention used when \code{date} is supplied.
#'   One of \code{"act/365"} or \code{"act/360"}.
#'
#' @return A tibble with one row per detected IRR and columns:
#' \describe{
#'   \item{root_id}{Root index.}
#'   \item{irr}{Detected IRR as an annual effective rate.}
#'   \item{i_effective_annual}{Same as \code{irr}, reported explicitly.}
#'   \item{j_nominal_interest}{Equivalent nominal annual interest rate
#'     convertible \code{m} times.}
#'   \item{delta}{Equivalent force of interest.}
#'   \item{npv}{NPV evaluated at the detected root, approximately zero.}
#'   \item{interval_left}{Left endpoint of the local search bracket.}
#'   \item{interval_right}{Right endpoint of the local search bracket.}
#'   \item{n_cashflows}{Length of \code{cf}.}
#'   \item{has_both_signs}{Whether the cash flow has at least one positive
#'     and one negative value.}
#'   \item{n_sign_changes_cashflow}{Number of sign changes in the nonzero
#'     cash-flow sequence.}
#' }
#'
#' If no roots are detected, the function returns a tibble with zero rows.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{cf} denotes cash flows, \code{t} denotes time,
#' and \code{m} denotes the conversion frequency used to report the equivalent
#' nominal annual interest rate.
#'
#' This function detects roots numerically over a finite search interval.
#' It may miss roots if:
#' \itemize{
#'   \item the grid is too coarse,
#'   \item two roots are extremely close,
#'   \item the NPV touches zero without changing sign,
#'   \item or the root lies outside the search interval.
#' }
#'
#' For a single-IRR workflow, use \code{\link{irr_flow}}.
#'
#' @seealso \code{\link{irr_flow}}, \code{\link{pv_flow}}
#'
#' @family time-value
#'
#' @examples
#' # A standard single-IRR cash flow
#' irr_flow_multi(
#'   cf = c(-1000, 300, 400, 500),
#'   t = c(0, 1, 2, 3)
#' )
#'
#' # A cash flow with multiple sign changes
#' irr_flow_multi(
#'   cf = c(-1000, 5000, -4500, 200),
#'   t = c(0, 1, 2, 3),
#'   search_interval = c(-0.99, 5),
#'   grid_points = 5000
#' )
#'
#' # Date-based version
#' irr_flow_multi(
#'   cf = c(-1000, 300, 400, 500),
#'   date = as.Date(c("2026-01-01", "2027-01-01", "2028-01-01", "2029-01-01"))
#' )
#'
#' @export
irr_flow_multi <- function(
    cf,
    t = NULL,
    date = NULL,
    m = 1L,
    search_interval = c(-0.99, 10),
    grid_points = 2000L,
    tol = 1e-10,
    maxiter = 1000L,
    day_count = c("act/365", "act/360")
) {
  day_count <- match.arg(day_count)

  if (!is.numeric(cf) || length(cf) < 2L) {
    stop("`cf` must be a numeric vector of length at least 2.", call. = FALSE)
  }

  if (any(is.na(cf)) || any(!is.finite(cf))) {
    stop("`cf` must contain only finite numeric values.", call. = FALSE)
  }

  if (!is.null(t) && !is.null(date)) {
    stop("Provide only one of `t` or `date`, not both.", call. = FALSE)
  }

  if (is.null(t) && is.null(date)) {
    stop("You must provide either `t` or `date`.", call. = FALSE)
  }

  if (!is.null(t)) {
    if (!is.numeric(t) || length(t) != length(cf)) {
      stop("`t` must be numeric and have the same length as `cf`.", call. = FALSE)
    }

    if (any(is.na(t)) || any(!is.finite(t)) || any(t < 0)) {
      stop("`t` must contain only finite values >= 0.", call. = FALSE)
    }
  }

  if (!is.null(date)) {
    date <- as.Date(date)

    if (length(date) != length(cf) || any(is.na(date))) {
      stop("`date` must contain valid dates and have the same length as `cf`.", call. = FALSE)
    }
  }

  if (!is.numeric(search_interval) ||
      length(search_interval) != 2L ||
      anyNA(search_interval)) {
    stop("`search_interval` must be a numeric vector of length 2.", call. = FALSE)
  }

  search_interval <- sort(search_interval)

  if (search_interval[1] <= -1) {
    stop("The lower bound of `search_interval` must be greater than -1.", call. = FALSE)
  }

  if (!is.numeric(grid_points) ||
      length(grid_points) != 1L ||
      is.na(grid_points) ||
      !is.finite(grid_points) ||
      grid_points < 2 ||
      grid_points != floor(grid_points)) {
    stop("`grid_points` must be an integer greater than or equal to 2.", call. = FALSE)
  }

  grid_points <- as.integer(grid_points)

  if (!is.numeric(m) ||
      length(m) != 1L ||
      is.na(m) ||
      !is.finite(m) ||
      m <= 0 ||
      m != floor(m)) {
    stop("`m` must be a positive integer.", call. = FALSE)
  }

  m <- as.integer(m)

  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol <= 0) {
    stop("`tol` must be a single positive numeric value.", call. = FALSE)
  }

  if (!is.numeric(maxiter) ||
      length(maxiter) != 1L ||
      is.na(maxiter) ||
      !is.finite(maxiter) ||
      maxiter <= 0 ||
      maxiter != floor(maxiter)) {
    stop("`maxiter` must be a positive integer.", call. = FALSE)
  }

  maxiter <- as.integer(maxiter)

  # --- Cash-flow sign diagnostics ---
  has_pos <- any(cf > 0)
  has_neg <- any(cf < 0)
  has_both <- has_pos && has_neg

  cf_nz <- cf[cf != 0]
  n_sign_changes <- if (length(cf_nz) <= 1L) {
    0L
  } else {
    sum(diff(sign(cf_nz)) != 0)
  }

  # Empty result schema (used for early returns)
  empty_result <- tibble::tibble(
    root_id = integer(0),
    irr = numeric(0),
    i_effective_annual = numeric(0),
    j_nominal_interest = numeric(0),
    delta = numeric(0),
    npv = numeric(0),
    interval_left = numeric(0),
    interval_right = numeric(0),
    n_cashflows = integer(0),
    has_both_signs = logical(0),
    n_sign_changes_cashflow = integer(0)
  )

  if (!has_both) {
    return(empty_result)
  }

  # --- NPV function ---
  npv_fun <- function(r) {
    pv_flow(
      cf = cf,
      i = r,
      i_type = "effective",
      t = t,
      date = date,
      day_count = day_count
    )
  }

  # --- Evaluate NPV on grid ---
  grid <- seq(search_interval[1], search_interval[2], length.out = grid_points)

  f_vals <- vapply(
    grid,
    function(r) {
      out <- tryCatch(npv_fun(r), error = function(e) NA_real_)
      if (is.finite(out)) out else NA_real_
    },
    numeric(1L)
  )

  # --- Candidate roots from near-zero grid evaluations ---
  zero_idx <- which(is.finite(f_vals) & abs(f_vals) <= sqrt(tol))
  point_roots <- if (length(zero_idx) > 0L) grid[zero_idx] else numeric(0)

  # --- Candidate brackets from sign changes (vectorized) ---
  ok <- is.finite(f_vals[-length(f_vals)]) & is.finite(f_vals[-1L])
  sign_change <- ok & (f_vals[-length(f_vals)] * f_vals[-1L] < 0)
  bracket_idx <- which(sign_change)

  bracket_roots <- numeric(0)
  bracket_left <- numeric(0)
  bracket_right <- numeric(0)

  if (length(bracket_idx) > 0L) {
    for (idx in bracket_idx) {
      a <- grid[idx]
      b <- grid[idx + 1L]

      root_obj <- tryCatch(
        stats::uniroot(
          f = npv_fun,
          interval = c(a, b),
          tol = tol,
          maxiter = maxiter
        ),
        error = function(e) NULL
      )

      if (!is.null(root_obj)) {
        bracket_roots <- c(bracket_roots, root_obj$root)
        bracket_left <- c(bracket_left, a)
        bracket_right <- c(bracket_right, b)
      }
    }
  }

  # --- Combine roots ---
  all_roots <- c(point_roots, bracket_roots)

  if (length(all_roots) == 0L) {
    return(empty_result)
  }

  # Build temporary results table before deduplication
  point_tbl <- if (length(point_roots) > 0L) {
    tibble::tibble(
      irr = point_roots,
      interval_left = point_roots,
      interval_right = point_roots
    )
  } else {
    tibble::tibble(
      irr = numeric(0),
      interval_left = numeric(0),
      interval_right = numeric(0)
    )
  }

  bracket_tbl <- if (length(bracket_roots) > 0L) {
    tibble::tibble(
      irr = bracket_roots,
      interval_left = bracket_left,
      interval_right = bracket_right
    )
  } else {
    tibble::tibble(
      irr = numeric(0),
      interval_left = numeric(0),
      interval_right = numeric(0)
    )
  }

  raw_tbl <- dplyr::bind_rows(point_tbl, bracket_tbl) |>
    dplyr::arrange(.data$irr)

  # --- Deduplicate numerically close roots ---
  keep <- rep(TRUE, nrow(raw_tbl))

  if (nrow(raw_tbl) > 1L) {
    for (j in 2:nrow(raw_tbl)) {
      if (abs(raw_tbl$irr[j] - raw_tbl$irr[j - 1L]) <= sqrt(tol)) {
        keep[j] <- FALSE
      }
    }
  }

  out_tbl <- raw_tbl[keep, , drop = FALSE]

  # --- Compute derived quantities ---
  out_tbl$root_id <- seq_len(nrow(out_tbl))
  out_tbl$i_effective_annual <- out_tbl$irr
  out_tbl$j_nominal_interest <- m * ((1 + out_tbl$irr)^(1 / m) - 1)
  out_tbl$delta <- log1p(out_tbl$irr)
  out_tbl$npv <- vapply(out_tbl$irr, npv_fun, numeric(1L))
  out_tbl$n_cashflows <- length(cf)
  out_tbl$has_both_signs <- has_both
  out_tbl$n_sign_changes_cashflow <- n_sign_changes

  out_tbl |>
    dplyr::select(
      "root_id",
      "irr",
      "i_effective_annual",
      "j_nominal_interest",
      "delta",
      "npv",
      "interval_left",
      "interval_right",
      "n_cashflows",
      "has_both_signs",
      "n_sign_changes_cashflow"
    )
}
