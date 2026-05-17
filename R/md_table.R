#' Multiple decrement table (annual, discrete ages)
#'
#' Builds a multiple decrement table from cause-specific annual decrement
#' probabilities \eqn{q_x^{(j)}}. This function is annual/discrete: ages must be
#' integer-valued and the input probabilities are interpreted as one-year
#' decrement probabilities for each cause.
#'
#' @details
#' Let the cause columns be \eqn{q_x^{(1)}, \dots, q_x^{(J)}}. The total decrement
#' probability is \eqn{q_x^{(\tau)} = \sum_j q_x^{(j)}} and the total survival
#' probability is \eqn{p_x^{(\tau)} = 1 - q_x^{(\tau)}}. The cohort is generated
#' recursively by \eqn{\ell_{x+1} = \ell_x \, p_x^{(\tau)}} with starting radix
#' \eqn{\ell_{x_0} = \text{radix}}.
#'
#' If \code{close = TRUE}, the last age (omega) must satisfy
#' \eqn{q_{\omega}^{(\tau)} = 1} (within tolerance), so that the table closes
#' naturally.
#'
#' @param qx_df A data.frame/tibble with an age column (default \code{x}) and one
#'   or more cause columns containing annual probabilities in \eqn{[0,1]}.
#'   Recommended naming convention: cause columns start with \code{"q_"}
#'   (e.g., \code{q_death}, \code{q_disability}).
#' @param age_col Character. Name of the age column (default \code{"x"}).
#' @param cause_cols Character vector. Names of the cause columns. If \code{NULL}
#'   (default), all columns other than \code{age_col} are treated as causes.
#' @param radix Numeric. Starting cohort size at the first age (default \code{1e5}).
#' @param close Logical. If \code{TRUE}, requires \eqn{q_{\omega}^{(\tau)} = 1}
#'   at the last age (default \code{TRUE}).
#' @param check Logical. If \code{TRUE}, performs input validation (default \code{TRUE}).
#' @param tol Numeric tolerance used in checks (default \code{1e-10}).
#'
#' @return A tibble with columns:
#' \itemize{
#'   \item \code{x}: integer ages.
#'   \item \code{lx}: cohort \eqn{\ell_x}.
#'   \item \code{q_total}: total decrement probability \eqn{q_x^{(\tau)}}.
#'   \item \code{p_total}: total survival probability \eqn{p_x^{(\tau)}}.
#'   \item \code{d_total}: total decrements \eqn{d_x^{(\tau)} = \ell_x q_x^{(\tau)}}.
#'   \item cause columns (as provided).
#'   \item cause-specific decrements \code{d_*} with \eqn{d_x^{(j)} = \ell_x q_x^{(j)}}.
#' }
#'
#' @examples
#' qx_df <- tibble::tibble(
#'   x = 30:35,
#'   q_death = c(0.001, 0.0012, 0.0014, 0.0017, 0.0020, 1.0000),
#'   q_disability = c(0.002, 0.0021, 0.0022, 0.0023, 0.0024, 0.0000)
#' )
#' md <- md_table(qx_df, radix = 1e5, close = TRUE)
#' md
#'
#' @export
md_table <- function(
    qx_df,
    age_col = "x",
    cause_cols = NULL,
    radix = 1e5,
    close = TRUE,
    check = TRUE,
    tol = 1e-10
) {

  if (missing(qx_df)) stop("`qx_df` is required.")
  if (!is.data.frame(qx_df)) stop("`qx_df` must be a data.frame/tibble.")
  if (!is.character(age_col) || length(age_col) != 1L) stop("`age_col` must be a single character string.")
  if (!(age_col %in% names(qx_df))) stop("`qx_df` must contain the age column specified by `age_col`.")
  if (!is.numeric(radix) || length(radix) != 1L || is.na(radix) || radix <= 0) stop("`radix` must be a single positive number.")

  df <- dplyr::as_tibble(qx_df) |>
    dplyr::arrange(.data[[age_col]])

  if (is.null(cause_cols)) {
    cause_cols <- setdiff(names(df), age_col)
  } else {
    if (!is.character(cause_cols)) stop("`cause_cols` must be a character vector.")
    missing_cols <- setdiff(cause_cols, names(df))
    if (length(missing_cols) > 0L) {
      stop("`cause_cols` contains columns not found in `qx_df`: ", paste(missing_cols, collapse = ", "), ".")
    }
  }

  if (length(cause_cols) == 0L) stop("At least one cause column is required.")

  # --- checks ---
  ages <- as.numeric(df[[age_col]])

  if (check) {
    if (any(!is.finite(ages))) stop("Ages must be finite.")
    if (any(abs(ages - round(ages)) > tol)) stop("Ages must be integer-valued for an annual multiple decrement table.")
    if (anyDuplicated(ages)) stop("Age column contains duplicated ages.")
    if (length(ages) > 1L && any(diff(as.integer(round(ages))) != 1L)) {
      stop("Ages must be consecutive integers (step = 1) for an annual table.")
    }

    for (cc in cause_cols) {
      if (!is.numeric(df[[cc]])) stop("Cause column `", cc, "` must be numeric.")
      if (any(!is.finite(df[[cc]]))) stop("Cause column `", cc, "` must be finite.")
      if (any(df[[cc]] < -tol | df[[cc]] > 1 + tol)) stop("Cause column `", cc, "` must be in [0,1].")
    }
  }

  # --- totals ---
  df <- df |>
    dplyr::mutate(
      q_total = rowSums(dplyr::across(dplyr::all_of(cause_cols)), na.rm = FALSE),
      p_total = 1 - .data$q_total
    )

  if (check) {
    if (any(df$q_total < -tol | df$q_total > 1 + tol)) stop("Row-wise sum of causes (q_total) must be in [0,1].")
    if (close) {
      q_last <- df$q_total[nrow(df)]
      if (abs(q_last - 1) > tol) {
        stop("With close = TRUE, require q_total at the last age (omega) to be 1 (within tol).")
      }
    }
  }

  # --- build lx recursively ---
  n <- nrow(df)
  lx <- numeric(n)
  lx[1] <- radix
  if (n > 1L) {
    for (k in 2:n) {
      lx[k] <- lx[k - 1] * df$p_total[k - 1]
    }
  }

  df <- df |>
    dplyr::mutate(
      lx = lx,
      d_total = .data$lx * .data$q_total
    )

  # --- cause-specific decrements d_x^(j) = lx * q_x^(j) ---
  # Name d columns by stripping a leading "q_" when present.
  d_names <- vapply(cause_cols, function(cc) {
    if (startsWith(cc, "q_")) {
      paste0("d_", substring(cc, 3))
    } else {
      paste0("d_", cc)
    }
  }, character(1))

  # Ensure uniqueness if user supplied tricky names
  d_names <- make.unique(d_names)

  for (j in seq_along(cause_cols)) {
    cc <- cause_cols[j]
    dn <- d_names[j]
    df[[dn]] <- df$lx * df[[cc]]
  }

  # --- output ordering ---
  out <- df |>
    dplyr::select(
      !!age_col,
      lx,
      q_total,
      p_total,
      d_total,
      dplyr::all_of(cause_cols),
      dplyr::all_of(d_names)
    ) |>
    dplyr::rename(x = !!age_col)

  out
}
