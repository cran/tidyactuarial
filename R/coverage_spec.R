#' Specify an insurance coverage contract
#' @md
#'
#' @description
#' Creates a lightweight S3 specification for per-loss insurance coverage,
#' including an ordinary or franchise deductible, a coverage limit,
#' coinsurance, inflation, and optional indexation of contract terms.
#'
#' @param deductible Nonnegative finite numeric scalar. The deductible applied
#'   to the inflated loss, possibly adjusted by `index_deductible`.
#' @param limit Nonnegative numeric scalar, or `Inf`. Maximum covered amount
#'   **after the deductible but before coinsurance**. This is not the maximum
#'   ground-up loss and not necessarily the insurer's final payment limit.
#' @param coinsurance Numeric scalar in `[0, 1]`. Fraction of the covered amount
#'   paid by the insurer, applied after the limit. Zero is allowed to represent
#'   a contract with no insurer participation.
#' @param inflation_factor Positive finite numeric scalar multiplying the
#'   original loss. This is a **factor**, not a rate: use `1.10` for 10 percent
#'   inflation.
#' @param deductible_type Either `"ordinary"` (default) or `"franchise"`.
#'   Under a franchise deductible, a loss exceeding the deductible is payable
#'   from the first monetary unit, subject to the limit and coinsurance.
#' @param index_deductible Logical scalar. If `TRUE`, multiply the deductible
#'   by `inflation_factor` before applying it to the inflated loss.
#' @param index_limit Logical scalar. If `TRUE`, multiply the limit by
#'   `inflation_factor` before applying it. An infinite limit remains infinite.
#' @param label Optional nonempty, non-whitespace character scalar identifying
#'   the contract.
#'
#' @details
#' This function only records contractual terms. It does not accept a severity
#' distribution, evaluate payments, simulate losses, or change a frequency model.
#'
#' Let `X >= 0` denote the original loss, `f` the inflation factor, `d` the
#' deductible, `L` the limit, and `c` the coinsurance fraction. Write
#' `Z = f * X`, `d_f = d * f` if the deductible is indexed (otherwise `d`),
#' and `L_f = L * f` if the limit is indexed (otherwise `L`). The specified
#' processing order is inflation, deductible, limit, and coinsurance.
#'
#' For an ordinary deductible, the intended payment per loss is
#' \deqn{Y_L = c \min\{\max(Z-d_f,0),L_f\}.}
#'
#' For a franchise deductible, the intended payment per loss is
#' \deqn{Y_L = c \min\{Z,L_f\} \mathbf{1}_{\{Z>d_f\}}.}
#' Thus a franchise loss exactly equal to its deductible produces no payment.
#' With a finite limit, the maximum insurer payment is `c * L_f` (or a smaller
#' achievable maximum if the loss distribution is bounded).
#'
#' In standard actuarial notation, the *maximum covered loss* for an ordinary
#' deductible with finite `d_f` and `L_f` is `u = d_f + L_f`; it is distinct
#' from this function's `limit` argument.
#'
#' All monetary amounts must use the same currency and reference monetary
#' basis. Inputs are deliberately scalar: one object represents one contract.
#'
#' @return
#' An object of class `tidyact_coverage_spec`, implemented as a named list
#' containing `deductible`, `limit`, `coinsurance`, `inflation_factor`,
#' `deductible_type`, `index_deductible`, `index_limit`, and `label`.
#' No derived payments or distributional quantities are stored.
#'
#' @references
#' Klugman, S. A., Panjer, H. H., and Willmot, G. E. (2012).
#' *Loss Models: From Data to Decisions*, 4th edition, Chapter 8.
#'
#' @examples
#' # Default: full insurance of each loss, without inflation.
#' coverage_spec()
#'
#' # Ordinary deductible with a capped covered layer and 80% coinsurance.
#' contract <- coverage_spec(
#'   deductible = 500,
#'   limit = 10000,
#'   coinsurance = 0.80,
#'   label = "Standard contract"
#' )
#' contract$limit
#'
#' # A qualifying loss is covered from its first unit under a franchise.
#' coverage_spec(deductible = 1000, deductible_type = "franchise")
#'
#' # Inflation with indexed deductible, but a fixed nominal limit.
#' coverage_spec(
#'   deductible = 500,
#'   limit = 10000,
#'   inflation_factor = 1.12,
#'   index_deductible = TRUE,
#'   index_limit = FALSE
#' )
#'
#' # Completely indexed contract, including an unbounded limit.
#' coverage_spec(
#'   deductible = 500,
#'   limit = Inf,
#'   inflation_factor = 1.05,
#'   index_deductible = TRUE,
#'   index_limit = TRUE
#' )
#'
#' @export
coverage_spec <- function(
  deductible = 0,
  limit = Inf,
  coinsurance = 1,
  inflation_factor = 1,
  deductible_type = c("ordinary", "franchise"),
  index_deductible = FALSE,
  index_limit = FALSE,
  label = NULL
) {
  check_nonnegative <- function(value, name, allow_infinity = FALSE) {
    valid <- is.numeric(value) &&
      length(value) == 1L &&
      !is.na(value) &&
      value >= 0 &&
      (is.finite(value) || (allow_infinity && is.infinite(value) && value > 0))

    if (!valid) {
      stop(
        sprintf("`%s` must be a nonnegative numeric scalar%s.",
                name, if (allow_infinity) " (or Inf)" else " (finite)"),
        call. = FALSE
      )
    }
  }

  check_nonnegative(deductible, "deductible")
  check_nonnegative(limit, "limit", allow_infinity = TRUE)
  check_nonnegative(coinsurance, "coinsurance")
  check_nonnegative(inflation_factor, "inflation_factor")

  if (coinsurance > 1) {
    stop("`coinsurance` must be between 0 and 1.", call. = FALSE)
  }
  if (inflation_factor == 0) {
    stop("`inflation_factor` must be strictly positive.", call. = FALSE)
  }

  # Preserve the exact published default while rejecting partial/ambiguous names.
  if (identical(deductible_type, c("ordinary", "franchise"))) {
    deductible_type <- "ordinary"
  }
  if (!is.character(deductible_type) ||
      length(deductible_type) != 1L ||
      is.na(deductible_type) ||
      !deductible_type %in% c("ordinary", "franchise")) {
    stop("`deductible_type` must be 'ordinary' or 'franchise'.", call. = FALSE)
  }

  is_flag <- function(x) is.logical(x) && length(x) == 1L && !is.na(x)
  if (!is_flag(index_deductible)) {
    stop("`index_deductible` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is_flag(index_limit)) {
    stop("`index_limit` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.null(label) &&
      (!is.character(label) || length(label) != 1L ||
       is.na(label) || !nzchar(trimws(label)))) {
    stop("`label` must be NULL or a nonempty character scalar.", call. = FALSE)
  }

  structure(
    list(
      deductible = as.numeric(deductible),
      limit = as.numeric(limit),
      coinsurance = as.numeric(coinsurance),
      inflation_factor = as.numeric(inflation_factor),
      deductible_type = deductible_type,
      index_deductible = index_deductible,
      index_limit = index_limit,
      label = label
    ),
    class = "tidyact_coverage_spec"
  )
}
