#' Total-decrement lifetable from a multiple decrement table: lt_tau
#'
#' @description
#' Builds a single-decrement lifetable for the *total* decrement (any cause),
#' using \eqn{q_x^{(\tau)}} from a multiple decrement table produced by
#' \code{\link{md_table}}. This enables direct re-use of single-life functions
#' (e.g., \code{t_px}, \code{t_qx}, \code{t_Ex}, annuities, insurances) under the
#' total decrement model.
#'
#' @details
#' Given cause-specific decrement probabilities \eqn{q_x^{(j)}}, the total decrement
#' is \eqn{q_x^{(\tau)} = \sum_j q_x^{(j)}}. This function simply passes
#' \code{x = md$x} and \code{qx = md$q_total} to \code{\link{lifetable}}.
#'
#' @param md A multiple decrement table (typically the output of \code{\link{md_table}}),
#'   containing columns \code{x} and \code{q_total}.
#' @param ... Additional arguments passed to \code{\link{lifetable}} (e.g., \code{radix},
#'   \code{omega}, \code{close}, \code{ax}, \code{type}, \code{frac}, \code{check}, \code{tol}).
#'
#' @return A lifetable object as produced by \code{\link{lifetable}}.
#'
#' @examples
#' qx_df <- tibble::tibble(
#'   x = 30:35,
#'   q_death = c(0.001, 0.0012, 0.0014, 0.0017, 0.0020, 1.0000),
#'   q_disability = c(0.002, 0.0021, 0.0022, 0.0023, 0.0024, 0.0000)
#' )
#' md <- md_table(qx_df, radix = 1e5, close = TRUE)
#' lt <- lt_tau(md, radix = 1e5, close = TRUE, frac = "UDD")
#' t_px(lt, x = 30, t = 5)
#'
#' @export
lt_tau <- function(md, ...) {
  if (missing(md)) stop("`md` is required.")
  if (!is.data.frame(md)) stop("`md` must be a data.frame/tibble.", call. = FALSE)
  if (!all(c("x", "q_total") %in% names(md))) {
    stop("`md` must contain columns `x` and `q_total`.", call. = FALSE)
  }

  # Delegate all modeling/validation details to lifetable()
  lifetable(x = md$x, qx = md$q_total, ...)
}
