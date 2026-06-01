#' Compute Monte Carlo net premiums for life contingencies
#'
#' Computes Monte Carlo net premiums from simulated present values of insurance
#' benefits and premium annuities, using compact actuarial notation.
#'
#' This function applies the actuarial equivalence principle to simulated
#' present value random variables. If \eqn{Z} denotes the simulated present
#' value of benefits and \eqn{Y} denotes the simulated present value of the
#' premium annuity, the Monte Carlo net premium is estimated as
#'
#' \deqn{
#'   \hat{P} = \frac{\bar{Z}}{\bar{Y}}.
#' }
#'
#' Equivalently, this estimates
#'
#' \deqn{
#'   P = \frac{E[Z]}{E[Y]}.
#' }
#'
#' @param .data A data frame or tibble containing simulated present values.
#'   Usually this object is obtained after applying \code{\link{mc_insurance}}
#'   and \code{\link{mc_annuity}} to the same simulated lifetime sample.
#' @param col_Z Character string. Name of the column containing the simulated
#'   present value of insurance benefits. Default is \code{"pv_benefit"}.
#' @param col_Y Character string. Name of the column containing the simulated
#'   present value of the premium annuity. Default is \code{"pv_annuity"}.
#' @param col_P Character string. Name of the output column containing the
#'   Monte Carlo net premium. Default is \code{"P"}.
#' @param by Optional character vector with grouping columns. If supplied, the
#'   premium is computed separately within each group. If \code{by = NULL} and
#'   \code{.data} is already grouped with \code{dplyr::group_by()}, the current
#'   grouping structure is used.
#' @param na_rm Logical scalar. Should missing values be removed when computing
#'   simulated means? Default is \code{TRUE}.
#' @param ... Transitional compatibility for older calls using \code{data},
#'   \code{benefit_col}, \code{annuity_col}, and \code{premium_col}.
#'
#' @details
#' This function follows the compact actuarial notation used throughout
#' \code{tidyactuarial}: \code{Z} denotes the present value random variable of
#' benefits, \code{Y} denotes the present value random variable of the premium
#' annuity, and \code{P} denotes the net premium per payment.
#'
#' The function does not simulate lifetimes and does not calculate present
#' values directly. It only computes the Monte Carlo net premium from columns
#' that already contain simulated present values.
#'
#' In a typical workflow, \code{\link{simulate_lifetime}} generates simulated
#' values of \eqn{K_x} and possibly \eqn{T_x}; \code{\link{mc_insurance}}
#' creates the simulated benefit present value \eqn{Z};
#' \code{\link{mc_annuity}} creates the simulated premium annuity present value
#' \eqn{Y}; and \code{mc_premium()} estimates the net level premium.
#'
#' The estimated premium is attached to every row of the input data. This is
#' intentional: it makes it easy to construct the simulated loss random
#' variable with \code{\link{mc_loss}}, for example
#'
#' \deqn{
#'   L = Z - \hat{P}Y.
#' }
#'
#' The function is also valid for premiums payable more than once per year. In
#' that case, the payment frequency is not specified in \code{mc_premium()};
#' it is already embedded in the simulated premium annuity present value
#' supplied through \code{col_Y}.
#'
#' For example, if \code{\link{mc_annuity}} was called with
#' \code{payment = 1 / 12} and \code{k = 12}, then \code{pv_annuity} represents
#' the present value of a monthly premium stream whose total annual amount is
#' 1. The premium estimated by \code{mc_premium()} is then consistent with that
#' monthly payment structure.
#'
#' If \code{\link{mc_annuity}} was called with \code{payment = 1} and
#' \code{k = 12}, then \code{pv_annuity} represents a stream of payments of 1
#' each month, and the resulting premium should be interpreted relative to that
#' payment pattern.
#'
#' This function computes net premiums only. It does not include expenses,
#' safety loadings, profit margins, taxes, surrender charges, commissions, or
#' other practical pricing adjustments.
#'
#' @return A tibble with the original columns and one additional column
#' containing the simulated net premium. The name of this column is controlled
#' by \code{col_P}. For transition, if \code{col_P != "premium"} and the input
#' does not already contain a column named \code{premium}, a legacy column
#' \code{premium} is also added with the same value.
#'
#' @seealso
#' \code{\link{simulate_lifetime}}, \code{\link{simulate_lifetimes}},
#' \code{\link{mc_insurance}}, \code{\link{mc_annuity}},
#' \code{\link{mc_loss}}, \code{\link{mc_reserve}},
#' \code{\link{summary_mc}}
#'
#' @references
#' Bowers, N. L., Gerber, H. U., Hickman, J. C., Jones, D. A.,
#' and Nesbitt, C. J. (1997). \emph{Actuarial Mathematics}. Second Edition.
#' Society of Actuaries.
#'
#' @family monte-carlo
#'
#' @examples
#' # Example 1: direct use with simulated present values
#' sim_values <- tibble::tibble(
#'   sim_id = 1:6,
#'   pv_benefit = c(0.82, 0.74, 0.61, 0.95, 0.70, 0.88),
#'   pv_annuity = c(8.2, 7.5, 6.1, 9.0, 7.2, 8.8)
#' )
#'
#' sim_values |>
#'   mc_premium()
#'
#' # Example 2: grouped premiums by age
#' sim_by_age <- tibble::tibble(
#'   sim_id = rep(1:6, times = 2),
#'   x = rep(c(40, 50), each = 6),
#'   pv_benefit = c(
#'     0.82, 0.74, 0.61, 0.95, 0.70, 0.88,
#'     0.91, 0.86, 0.79, 0.98, 0.83, 0.94
#'   ),
#'   pv_annuity = c(
#'     8.2, 7.5, 6.1, 9.0, 7.2, 8.8,
#'     6.8, 6.4, 5.9, 7.1, 6.2, 6.7
#'   )
#' )
#'
#' sim_by_age |>
#'   mc_premium(by = "x")
#'
#' # Example 3: using dplyr grouping
#' sim_by_age |>
#'   dplyr::group_by(x) |>
#'   mc_premium()
#'
#' # Example 4: annual whole life net premium
#' lt <- tibble::tibble(
#'   x = 40:100,
#'   qx = seq(0.002, 1, length.out = 61)
#' )
#'
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     seed = 123
#'   ) |>
#'   mc_insurance(
#'     i = 0.05,
#'     type = "whole",
#'     benefit = 1
#'   ) |>
#'   mc_annuity(
#'     i = 0.05,
#'     type = "whole",
#'     payment = 1,
#'     k = 1,
#'     timing = "due"
#'   ) |>
#'   mc_premium()
#'
#' # Example 5: monthly whole life net premium
#' lt |>
#'   simulate_lifetime(
#'     x = 40,
#'     n_sim = 25,
#'     frac = "udd",
#'     seed = 123
#'   ) |>
#'   mc_insurance(
#'     i = 0.05,
#'     type = "whole",
#'     benefit = 1
#'   ) |>
#'   mc_annuity(
#'     i = 0.05,
#'     type = "whole",
#'     payment = 1 / 12,
#'     k = 12,
#'     timing = "due"
#'   ) |>
#'   mc_premium()
#'
#' # Example 6: term insurance with monthly premium annuity
#' lt |>
#'   simulate_lifetime(
#'     x = 45,
#'     n_sim = 25,
#'     frac = "udd",
#'     seed = 123
#'   ) |>
#'   mc_insurance(
#'     i = 0.04,
#'     type = "term",
#'     n = 20,
#'     benefit = 100000
#'   ) |>
#'   mc_annuity(
#'     i = 0.04,
#'     type = "temporary",
#'     n = 20,
#'     payment = 1 / 12,
#'     k = 12,
#'     timing = "due"
#'   ) |>
#'   mc_premium()
#'
#' # Example 7: joint-life annual net premium
#' lt |>
#'   simulate_lifetimes(
#'     x = c(60, 58),
#'     n_sim = 25,
#'     seed = 123
#'   ) |>
#'   mc_multilife_status(status = "joint") |>
#'   mc_insurance(
#'     i = 0.04,
#'     type = "whole",
#'     benefit = 100000,
#'     col_K = "K_status",
#'     col_T = "T_status"
#'   ) |>
#'   mc_annuity(
#'     i = 0.04,
#'     type = "whole",
#'     payment = 1,
#'     k = 1,
#'     timing = "due",
#'     col_K = "K_status",
#'     col_T = "T_status"
#'   ) |>
#'   mc_premium()
#'
#' # Transitional compatibility with old column arguments
#' sim_values |>
#'   mc_premium(
#'     benefit_col = "pv_benefit",
#'     annuity_col = "pv_annuity",
#'     premium_col = "premium"
#'   )
#'
#' @export
mc_premium <- function(
    .data = NULL,
    col_Z = "pv_benefit",
    col_Y = "pv_annuity",
    col_P = "P",
    by = NULL,
    na_rm = TRUE,
    ...
) {
  dots <- list(...)

  # -------------------------------------------------------------------------
  # Transitional compatibility with the previous public API
  # -------------------------------------------------------------------------

  allowed_old <- c(
    "data",
    "benefit_col",
    "annuity_col",
    "premium_col"
  )

  bad_dots <- setdiff(names(dots), allowed_old)

  if (length(bad_dots) > 0L) {
    stop(
      "Unused argument(s): ",
      paste(sprintf("`%s`", bad_dots), collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!is.null(dots$data)) {
    if (!is.null(.data)) {
      stop("Provide only one of `.data` or deprecated `data`.", call. = FALSE)
    }

    .data <- dots$data
  }

  if (!is.null(dots$benefit_col)) {
    if (!identical(col_Z, "pv_benefit")) {
      stop("Provide only one of `col_Z` or deprecated `benefit_col`.",
           call. = FALSE)
    }

    col_Z <- dots$benefit_col
  }

  if (!is.null(dots$annuity_col)) {
    if (!identical(col_Y, "pv_annuity")) {
      stop("Provide only one of `col_Y` or deprecated `annuity_col`.",
           call. = FALSE)
    }

    col_Y <- dots$annuity_col
  }

  if (!is.null(dots$premium_col)) {
    if (!identical(col_P, "P")) {
      stop("Provide only one of `col_P` or deprecated `premium_col`.",
           call. = FALSE)
    }

    col_P <- dots$premium_col
  }

  # -------------------------------------------------------------------------
  # Validation
  # -------------------------------------------------------------------------

  if (!is.data.frame(.data)) {
    stop("`.data` must be a data frame or tibble.", call. = FALSE)
  }

  .mc_assert_numeric_column(.data, col_Z, "col_Z")
  .mc_assert_numeric_column(.data, col_Y, "col_Y")
  .mc_assert_character_scalar(col_P, "col_P")

  if (!is.logical(na_rm) || length(na_rm) != 1L || is.na(na_rm)) {
    stop("`na_rm` must be a logical scalar.", call. = FALSE)
  }

  if (!is.null(by)) {
    if (!is.character(by) || anyNA(by)) {
      stop("`by` must be NULL or a character vector.", call. = FALSE)
    }

    if (!all(by %in% names(.data))) {
      stop("All columns supplied in `by` must exist in `.data`.", call. = FALSE)
    }
  }

  active_by <- by

  if (is.null(active_by)) {
    active_by <- dplyr::group_vars(.data)
  }

  compute_P <- function(Z, Y) {
    denominator <- mean(Y, na.rm = na_rm)

    if (is.na(denominator) || denominator <= 0) {
      stop(
        "The mean simulated premium annuity present value must be positive.",
        call. = FALSE
      )
    }

    numerator <- mean(Z, na.rm = na_rm)

    if (is.na(numerator)) {
      stop(
        "The mean simulated benefit present value could not be computed.",
        call. = FALSE
      )
    }

    numerator / denominator
  }

  add_legacy_premium <- function(out) {
    if (!identical(col_P, "premium") && !"premium" %in% names(out)) {
      out <- dplyr::mutate(
        out,
        premium = .data[[col_P]]
      )
    }

    out
  }

  if (length(active_by) == 0L) {
    P_value <- compute_P(
      Z = .data[[col_Z]],
      Y = .data[[col_Y]]
    )

    out <- .data |>
      dplyr::ungroup() |>
      dplyr::mutate(
        "{col_P}" := P_value
      )

    return(add_legacy_premium(out))
  }

  out <- .data |>
    dplyr::ungroup() |>
    dplyr::group_by(dplyr::across(dplyr::all_of(active_by))) |>
    dplyr::mutate(
      "{col_P}" := compute_P(
        Z = .data[[col_Z]],
        Y = .data[[col_Y]]
      )
    ) |>
    dplyr::ungroup()

  add_legacy_premium(out)
}
