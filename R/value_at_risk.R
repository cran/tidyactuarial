#' Value at Risk for actuarial losses
#'
#' Computes the lower quantile VaR at level p: the smallest x for which
#' F_X(x) >= p. No interpolation across probability atoms; for discrete
#' claims, the quantile is an observed support point, not an interpolated
#' monetary amount.
#'
#' @param x A vector of observed losses, a tibble of simulations, a severity
#'   specification with an exact quantile, or an [aggregate_distribution()]
#'   result. A zero-payment specification is also accepted.
#' @param level Probability levels strictly between zero and one.
#' @param value Optional character name of a loss column in a data frame.
#'   By default, `aggregate_loss`, `claim_amount`, `claim_count`, `loss`, or
#'   `payment` is selected (in that order).
#' @return A tibble with `level`, `value_at_risk`, and `method`.
#' @details The tail probability retained by an aggregate-distribution object
#'   is never renormalized away. Quantiles in its uncomputed tail are rejected.
#'   For sampled losses the definition is empirical, not a model-based claim.
#'   FFT count truncation can omit probability at small losses as well as
#'   large losses. If the omitted count mass could move the selected quantile
#'   to a smaller amount, VaR is rejected; recompute with a smaller `tol`.
#' @examples
#' value_at_risk(c(0, 100, 100, 300), level = c(0.5, 0.9))
#' value_at_risk(severity_spec("exponential", rate = 0.01), level = 0.95)
#' @export
value_at_risk <- function(x, level = 0.95, value = NULL) {
  level <- .risk_validate_levels(level)
  if (inherits(x, "tidyact_severity_spec") &&
      !x$dist %in% c("discrete", "empirical", "degenerate")) {
    if (!is.null(value)) stop("`value` only selects a data-frame column.",
                              call. = FALSE)
    quantile <- .risk_exact_quantile(x, level)
    method <- "exact_quantile"
  } else {
    masses <- .risk_masses(x, value)
    quantile <- .risk_discrete_var(masses, level)
    method <- masses$method
  }
  tibble::tibble(level = level, value_at_risk = quantile,
                 method = rep(method, length(level)))
}
