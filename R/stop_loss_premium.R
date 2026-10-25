#' Net stop-loss premium
#'
#' Computes E[(X-d)+], the expected excess loss over retention d, with
#' no loading or discounting. This is a per-risk expectation for X.
#' Equivalently, the premium is the integral of the survival function
#' from the retention to infinity, provided the expectation is finite.
#'
#' @param x A nonnegative loss sample or simulation tibble, a supported
#'   severity specification, or an aggregate-distribution result.
#' @param retention Finite nonnegative retention amount(s).
#' @param value Optional name of a loss column when `x` is a data frame.
#' @return A tibble with `retention`, `premium`, and `method`.
#' @details Closed forms are available for exponential, gamma, lognormal,
#'   Weibull, and Pareto II severity. Burr XII, log-logistic, and bounded
#'   custom distributions with a survival provider use explicitly identified
#'   numerical integration of the survival function. Finite discrete models and observations
#'   are summed exactly over their mass distributions. A missing aggregate
#'   tail forbids an exact stop-loss expectation; no renormalization or zero
#'   assumption is made. The Pareto II premium is infinite when shape <= 1.
#'   Per-loss exponential severities transformed by [apply_coverage()] use
#'   closed forms for ordinary and franchise deductibles, including unlimited
#'   coverage. This avoids numerical integration over very large limits.
#' @examples
#' stop_loss_premium(c(0, 100, 200, 500), retention = c(0, 150))
#' stop_loss_premium(severity_spec("exponential", rate = 0.01), retention = 50)
#' @export
stop_loss_premium <- function(x, retention = 0, value = NULL) {
  retention <- .risk_validate_retention(retention)
  if (inherits(x, "tidyact_severity_spec") &&
      !x$dist %in% c("discrete", "empirical", "degenerate")) {
    if (!is.null(value)) stop("`value` only selects a data-frame column.",
                              call. = FALSE)
    premium <- vapply(retention, function(d) .risk_stop_loss_exact(x, d),
                      numeric(1))
    method <- .risk_stop_loss_method(x)
  } else {
    masses <- .risk_masses(x, value)
    if (masses$tail > 0) {
      stop("The stop-loss premium is unavailable: uncomputed aggregate tail.",
           call. = FALSE)
    }
    premium <- vapply(retention, function(d) {
      sum(pmax(masses$loss - d, 0) * masses$probability)
    }, numeric(1))
    method <- masses$method
  }
  tibble::tibble(retention = retention, premium = premium,
                 method = rep(method, length(retention)))
}
