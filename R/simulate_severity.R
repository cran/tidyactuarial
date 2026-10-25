#' Simulate individual claim severities
#'
#' Samples independent nonnegative claim amounts from an actuarial severity
#' specification, including atoms at zero when present. A realization is one individual
#' claim amount, not an aggregate or a fitted severity estimate.
#'
#' @param x A [severity_spec()], [empirical_severity()], or compatible payment
#'   severity specification.
#' @param n_sim Positive integer number of simulations.
#' @param seed Optional nonnegative integer seed; preserves the caller's RNG
#'   state when supplied.
#' @return A tibble with `simulation` and `claim_amount` columns.
#' @details A missing RNG provider produces an explicit error. Zero severities
#'   are retained; no implicit conversion from per-loss to per-payment basis
#'   takes place.
#' @examples
#' simulate_severity(severity_spec("exponential", rate = 0.002),
#'                   n_sim = 5, seed = 20)
#' @export
simulate_severity <- function(x, n_sim = 1000L, seed = NULL) {
  .collective_validate_severity(x)
  n_sim <- .risk_check_n_sim(n_sim)
  .risk_with_seed(seed, {
    tibble::tibble(
      simulation = seq_len(n_sim),
      claim_amount = .risk_draw(x, n_sim, "severity")
    )
  })
}
